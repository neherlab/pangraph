import os
import io

from glob   import glob
from Bio    import Entrez, SeqIO
from ftplib import FTP

import pyfaidx as fai

# ------------------------------------------------------------------------
# Global constants

FAPATH   = "data/seq/all.fasta"
PREFIX   = "ftp://"
NCBI_FTP = "ftp.ncbi.nlm.nih.gov"

# ------------------------------------------------------------------------
# Functions

def findid(term, db="nucleotide"):
    handle = Entrez.esearch(db=db, term=term)
    record = Entrez.read(handle)
    return record['IdList']

def getseq(uid, out):
    handle = Entrez.efetch(db="assembly", id=uid[0], rettype="docsum", retmod="xml")
    record = Entrez.read(handle)

    # Get path to relevant FTP site.
    url = record['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_RefSeq']
    assert url.startswith(PREFIX + NCBI_FTP)
    url = url[len(PREFIX) + len(NCBI_FTP) + 1:].split("/")
    fa  = url[-1] + "_genomic.fna.gz"
    url = "/".join(url)

    # Store fasta file in memory.
    ftp = FTP(NCBI_FTP)
    ftp.login()
    ftp.cwd(url)
    ftp.retrbinary(f"RETR {fa}", open(out, 'wb').write)
    ftp.quit()

def main():
    FA = fai.Fasta(FAPATH)
    for isolate in FA.keys():
        out = f"data/genomes/{isolate}.fa.gz"
        if 'carb' in isolate or os.path.exists(out):
            continue

        print(f"Grabbing {isolate}")
        if isolate.startswith("GCA_"):
            terms = isolate.split("_")
            isolate = terms[0] + "_" + terms[1] + "." + terms[2]

            uid = findid(isolate, db="assembly")
            assert len(uid) >= 1

            getseq(uid, out)
        else:
            handle = Entrez.esearch(db="nucleotide", term=isolate)
            record = Entrez.read(handle)
            ID     = record['IdList']
            assert len(ID) >= 1

            handle = Entrez.efetch(db="nucleotide", id=ID[0], rettype="gb", retmod="text")
            gbk    = SeqIO.read(handle, "genbank")

            assembly = None
            for ref in gbk.dbxrefs:
                typ = ref.split(":")
                if typ[0] == "Assembly":
                    assembly = typ[1]

            assert assembly is not None
            uid = findid(assembly, db="assembly")
            assert len(uid) >= 1

            getseq(uid, out)

if __name__ == "__main__":
    Entrez.email = "nnoll523@gmail.com"
    main()
