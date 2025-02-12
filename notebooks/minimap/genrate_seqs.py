# %%
import numpy as np
from Bio import SeqIO, Seq


def generate_seq(L):
    return "".join(np.random.choice(["A", "C", "G", "T"], L))


def mutate(seq, pm, pi, pd):
    seq = np.array(list(seq))
    L = len(seq)
    M_mask = np.random.rand(L) < pm
    D_mask = np.random.rand(L) < pd
    I_mask = np.random.rand(L) < pi

    # perform mutations
    for l in np.where(M_mask)[0]:
        allow = set(["A", "C", "G", "T"]) - set(seq[l])
        seq[l] = np.random.choice(list(allow))

    # perform deletions
    seq[D_mask] = ""

    # perform insertions
    for l in np.where(I_mask)[0]:
        seq[l] = seq[l] + np.random.choice(["A", "C", "G", "T"])

    return "".join(seq)


L = 1000
N = 5
pm, pi, pd = 0.02, 0.002, 0.002
prc = 0.5


np.random.seed(0)
refs = []
qrys = []
for n in range(N):
    ref = generate_seq(L)
    qry = mutate(ref, pm, pi, pd)

    ref = SeqIO.SeqRecord(Seq.Seq(ref), id=f"ref_{n}", description="")
    qry = SeqIO.SeqRecord(Seq.Seq(qry), id=f"qry_{n}", description="")

    if np.random.rand() < prc:
        qry.seq = qry.seq.reverse_complement()

    refs.append(ref)
    qrys.append(qry)

SeqIO.write(refs, f"test_data/ref.fa", "fasta")
SeqIO.write(qrys, f"test_data/qry.fa", "fasta")

# %%

# align with minimap
import subprocess

subprocess.run(
    "minimap2 -c -x asm20 -k 10 test_data/ref.fa test_data/qry.fa > test_data/aln.paf",
    shell=True,
)
