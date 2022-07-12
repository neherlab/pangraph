# Utilities to parse raw data and extract relevant features.

import numpy as np


def initialize_genmap(seqlen, seqdescr):
    """initializes an empty gmap object, which contains information on a genome"""

    gmap = {
        "gid": [],  # gene id
        "gbeg": [],  # gene beginning
        "glen": [],  # gene length
        "strand": [],  # which strand (+- 1)
        "ann": [],  # annotation
        "accnum": [],  # accession number
        "slen": seqlen,  # total strain sequence length
        "sdescr": seqdescr,  # sequence annotation
    }
    return gmap


def gmap_insert(gmap, gene_id, start, length, strand, ann, accnum):
    """insert the new gene in the correct position, based on gene beginning"""
    idx = np.searchsorted(gmap["gbeg"], start)
    gmap["gid"].insert(idx, gene_id)
    gmap["gbeg"].insert(idx, start)
    gmap["glen"].insert(idx, length)
    gmap["strand"].insert(idx, strand)
    gmap["ann"].insert(idx, ann)
    gmap["accnum"].insert(idx, accnum)


def tag_dict(record):
    """extract information from a record. All features of type 'gene' are
    extracted. They are saved in a tag dictionary whose key is the 'locus_tag'
    property, and which contains gene start, gene length, and strand location."""
    ftrs = [ft for ft in record.features if ft.type == "gene"]
    tags = {}
    for ft in ftrs:
        tag = ft.qualifiers["locus_tag"][0]
        tags[tag] = [
            int(ft.location.start),
            len(ft.location),
            ft.strand,
        ]
        assert ft.location.start <= ft.location.end, "start should be smaller than end"

    return tags


def register_record(record, genes_df):
    """given a genebank file record and the genes_df dataframe, this function
    extracts entries from the file according to matches in genes_df.
    It returns a dictionary with the features, and a list of missed genes.
    This list consists of indices of the genes_df dataframe, so that they can
    later be removed"""

    # list of items in genes_df to be dropped due to them being not present
    missed_idxs = []

    # sequence length and description
    L = len(record.seq)
    seqdescr = record.description

    # initialize empty gene map
    gmap = initialize_genmap(L, seqdescr)

    # extract tag dictionary
    tgdict = tag_dict(record)
    assert len(tgdict) > 0, "error: no annotations"

    # extract from the dataframe the list of relevant loci and genes
    strain = record.name
    mask = genes_df["str"] == strain
    sdf = genes_df[mask]

    # for every gene in the list, append gene to the map
    for index, row in sdf.iterrows():
        gid = row["geneId"]
        tag = row["tag"]
        ann = row["ann"]
        if tag in tgdict:
            start, length, strand = tgdict[tag]
            gmap_insert(gmap, gid, start, length, strand, ann, tag)
        else:
            missed_idxs.append(index)

    return gmap, missed_idxs
