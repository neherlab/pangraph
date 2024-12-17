from dataclasses import dataclass
from .utils import Edit, Block


# see: https://github.com/neherlab/pangraph/blob/98886771cb20cd4bfe7ce33c52dafc2fc33f6faa/packages/pangraph/src/pangraph/edits.rs#L194
# I added the aligned optional argument.
# When true it will not add insertions and have deletions as gaps
# if this is not good practice could also be split into two functions,
# or have an enum for the output type (aligned vs whole sequence)
def apply_edits_to_ref(edits: Edit, ref: str, aligned=False) -> str:
    """
    Apply the edits to the reference to obtain the query sequence.
    If aligned is True, insertions are not added and deletions are replaced by gaps.
    """
    qry = list(ref)
    for S in edits.subs:
        qry[S.pos] = S.alt

    for D in edits.dels:
        for l in range(D.length):
            qry[D.pos + l] = "-" if aligned else ""
    if aligned:
        qry = [q for q in qry if q != "-"]
        for I in edits.ins:
            if I.pos > 0:
                qry[I.pos - 1] += I.ins
            elif I.pos == 0:
                qry[0] = I.ins + qry[0]
    return "".join(qry)


# reimplementation of already existing code in pangraph
# see https://github.com/neherlab/pangraph/blob/rust/packages/pangraph/src/io/fasta.rs
@dataclass
class FastaRecord:
    name: str
    seq: str
    idx: int


def write_to_file(records: list[FastaRecord], file_path: str):
    with open(file_path, "w") as f:
        for record in sorted(records, key=lambda x: x.idx):
            f.write(f">{record.name}\n{record.seq}\n")
