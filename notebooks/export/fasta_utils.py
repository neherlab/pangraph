from dataclasses import dataclass
from .utils import Edit, Pangraph, Block
from Bio import Seq


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
            # add deletions as gaps if aligned
            # otherwise remove the base
            qry[D.pos + l] = "-" if aligned else ""

    # only add insertions if sequences are not aligned
    if not aligned:
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

    def reverse_complement(self):
        """
        Returns a new FastaRecord object with the reverse complement of the sequence.
        """
        return FastaRecord(
            name=self.name, seq=Seq.reverse_complement(self.seq), idx=self.idx
        )


def concatenate_records(records_list: list[list[FastaRecord]]) -> list[FastaRecord]:
    """
    Given a list of lists of FastaRecord objects, concatenates all records
    with the same name in the same sequence, and returns a single list of FastaRecord objects.
    """
    records = {record.name: "" for record in records_list[0]}
    for entry in records_list:
        for record in entry:
            records[record.name] += record.seq
    return [
        FastaRecord(name=name, seq=seq, idx=idx)
        for idx, (name, seq) in enumerate(records.items())
    ]


def write_to_file(records: list[FastaRecord], file_path: str):
    with open(file_path, "w") as f:
        for record in sorted(records, key=lambda x: x.idx):
            f.write(f">{record.name}\n{record.seq}\n")


def read_from_file(fname: str) -> list[FastaRecord]:
    """
    Given a fasta file, reads the records and returns a list of FastaRecord objects.
    """
    records = []
    with open(fname, "r") as f:
        lines = f.readlines()
    name = None
    seq = ""
    for line in lines:
        if line.startswith(">"):
            if name is not None:
                records.append(FastaRecord(name=name, seq=seq, idx=len(records)))
            name = line.strip()[1:]
            seq = ""
        else:
            seq += line.strip()

    if name is not None:
        records.append(FastaRecord(name=name, seq=seq, idx=len(records)))

    return records


def create_node_record_id(graph: Pangraph, node_id: int) -> str:
    """
    Given a pangraph and a node id, returns the record id for the export fasta file.
    This is in the format: 'node_id path_name-block_id [start-end|strand]'
    """
    node = graph.nodes[node_id]
    block_id = node.block_id
    path_name = graph.paths[node.path_id].name
    start, end = node.position
    strand = "+" if node.strand else "-"
    return f"{node_id} {path_name}-{block_id} [{start}-{end}|{strand}]"


# nb: record_naming could be an enum in rust...
def block_to_aln(
    graph: Pangraph, block: Block, aligned: bool, record_naming: str
) -> list[FastaRecord]:
    """
    Given a block, returns a list of FastaRecord objects containing a sequence per node.
    If aligned is True, it returns aligned sequences, with gaps for deletions and no insertions.
    If aligned is False, it returns the full unaligned sequences.
    Record name can be either "node", in which case the record id is the node id (plus some extra info),
    or "path", in which case the record id is the path name.
    """
    records = []
    for idx, (node_id, edits) in enumerate(block.alignments.items()):
        match record_naming:
            case "node":
                record_id = create_node_record_id(graph, node_id)
            case "path":
                path_id = graph.nodes[node_id].path_id
                record_id = graph.paths[path_id].name
            case _:
                raise ValueError(
                    "Invalid record naming scheme, either 'node' or 'path'"
                )
        seq = apply_edits_to_ref(edits, block.consensus, aligned=aligned)
        records.append(FastaRecord(name=record_id, seq=seq, idx=idx))
    return records
