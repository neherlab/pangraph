import glob
import json

PX_config = config["panx"]
# list of different species
PX_species = PX_config["species"]
# kernel options
PX_ker_opt = PX_config["kernel-options"]
# different kernel names
PX_ker_names = list(PX_ker_opt.keys())


wildcard_constraints:
    species=f"({'|'.join(PX_species)})",
    kind=f"({'|'.join(PX_ker_names)})",
    s1="[^_-][^-]+[^_-]",
    s2="[^_-][^-]+[^_-]",


# ------------- Compare pairwise graphs vs pairwise projections from a larger graph -----------
# load dictionary of selected strains & selected pairs for analysis (produced by pick_projection_strain_set.py)


with open(PX_config["projection-data"], "r") as f:
    MG_opt = json.load(f)


# Creates a graph with all of the selected strains.
rule MG_projection_full_graph:
    message:
        "Creating projection graph for species {wildcards.species}"
    input:
        lambda w: expand(
            "panx_data/{{species}}/fa/{acc}.fa", acc=MG_opt[w.species]["strains"]
        ),
    output:
        "projections/{species}/full/pangraph_{kind}.json",
    params:
        opt=lambda w: PX_ker_opt[w.kind],
    conda:
        "../conda_envs/pangraph_build_env.yml"
    shell:
        """
        export JULIA_NUM_THREADS=8
        pangraph build --circular {params.opt} {input} > {output}
        """


# Builds a graph for every selected pair of strains.
rule MG_pairwise_graphs:
    message:
        "Building pairwise graph for strains {wildcards.s1} - {wildcards.s2} ({wildcards.species})"
    input:
        "panx_data/{species}/fa/{s1}.fa",
        "panx_data/{species}/fa/{s2}.fa",
    output:
        "projections/{species}/pairwise/pangraph_{kind}__{s1}-{s2}.json",
    params:
        opt=lambda w: PX_ker_opt[w.kind],
    conda:
        "../conda_envs/pangraph_build_env.yml"
    shell:
        """
        pangraph build --circular {params.opt} {input} > {output}
        """


# For every selected pair of strains, marginalizes the global graph on the pair.
rule MG_pairwise_projection:
    message:
        "project graph on strains {wildcards.s1} - {wildcards.s2} ({wildcards.species})"
    input:
        rules.MG_projection_full_graph.output,
    output:
        "projections/{species}/projected/pangraph_{kind}__{s1}-{s2}.json",
    conda:
        "../conda_envs/pangraph_build_env.yml"
    shell:
        """
        pangraph marginalize -s {wildcards.s1},{wildcards.s2} {input} > {output}
        """


rule MG_shared_kmers_pair:
    message:
        "species {wildcards.species} - evaluating the number of shared kmers {wildcards.s1}|{wildcards.s2}"
    input:
        s1="panx_data/{species}/fa/{s1}.fa",
        s2="panx_data/{species}/fa/{s2}.fa",
    output:
        temp("projections/{species}/kmers/{s1}-{s2}.json"),
    params:
        k=PX_config["kmer-size"],
    conda:
        "../conda_envs/bioinfo_env.yml"
    shell:
        """
        python3 workflow_scripts/shared_kmers_pair.py \
            --s1 {input.s1} --s2 {input.s2} --k {params.k} --json {output}
        """


def all_kmer_pairs(w):
    """given species and kind wildcards produces a list containing all pairs
    of json files with pairwise statistics comparisons."""
    pairs = MG_opt[w.species]["pairs"]
    in_files = []
    for s1, s2 in pairs:
        in_files.append(f"projections/{w.species}/kmers/{s1}-{s2}.json")
    return in_files


rule MG_shared_kmers_summary:
    message:
        "evaluating summary of shared k-mers for species {wildcards.species}"
    input:
        all_kmer_pairs,
    output:
        csv="projections/{species}/kmers/shared_kmers.csv",
    run:
        import pandas as pd
        import json

        data = []
        for i in input:
            with open(i, "r") as f:
                data += json.load(f)
        pd.DataFrame(data).to_csv(output.csv, index=False)



rule MG_compare_projection_pairwise:
    # Compares the pairwise and the marginalized graph, checking on what fraction of the genome
    # the two agree. Results are saved in a json file
    message:
        "comparing projected graph to pairwise graph ({wildcards.s1} - {wildcards.s2} ; {wildcards.species})"
    input:
        pw=rules.MG_pairwise_graphs.output,
        pj=rules.MG_pairwise_projection.output,
    output:
        "projections/{species}/comparison/{kind}__{s1}-{s2}.json",
    shell:
        """
        julia -t 1 --project=. workflow_scripts/pairwise_vs_marginalize.jl\
            {input.pw} {input.pj} {output}
        """


def all_pair_comparisons(w):
    """given species and kind wildcards produces a list containing all pairs
    of json files with pairwise statistics comparisons."""
    pairs = MG_opt[w.species]["pairs"]
    in_files = []
    for s1, s2 in pairs:
        in_files.append(f"projections/{w.species}/comparison/{w.kind}__{s1}-{s2}.json")
    return in_files


rule MG_pairwise_projection_benchmark:
    # For every species and alignment kernel collects all the statistics on the comparison
    # between pairwise and projeted graphs, and produces summary plots.
    message:
        "preparing plots for pairwise comparison - {wildcards.species} - {wildcards.kind}"
    input:
        jsons=all_pair_comparisons,
    output:
        pdf="projections/benchmark/{kind}_{species}.pdf",
        csv_summary="projections/benchmark/{kind}_{species}.summary.csv",
        csv_full="projections/benchmark/{kind}_{species}.full.csv",
    conda:
        "../conda_envs/bioinfo_env.yml"
    shell:
        """
        python3 workflow_scripts/pairwise_vs_marginalize_summary.py \
            --jsons {input.jsons} \
            --csv_summary {output.csv_summary} --csv_full {output.csv_full} --pdf {output.pdf} \
            --species {wildcards.species}
        """


rule MG_explore_proj_examples:
    # Produces detailed plots for specific examples of comparison between pairwise
    # and projected graphs.
    message:
        "plot inspecting example {wildcards.species} - {wildcards.kind} - {wildcards.s1} - {wildcards.s2} for pairwise comparison"
    input:
        pw=rules.MG_pairwise_graphs.output,
        pj=rules.MG_pairwise_projection.output,
    output:
        pdf="projections/benchmark/examples/{species}__{kind}__{s1}-{s2}-{title}.pdf",
    conda:
        "../conda_envs/bioinfo_env.yml"
    shell:
        """
        python3 workflow_scripts/pairwise_vs_marginalize_example.py \
            --pair {input.pw} --proj {input.pj} --pdf {output.pdf}
        """


# list of example plots with comparison of pairwise and projected graphs
MG_example_list = [
    ("escherichia_coli", "NZ_CP014316-NZ_CP018983", "max_agreeF_max_disagreeL"),
    ("escherichia_coli", "NZ_CP012635-NZ_CP017251", "min_agreeF"),
    ("escherichia_coli", "NZ_CP018252-NZ_CP018983", "min2_agreeF"),
    ("escherichia_coli", "NZ_CP017249-NZ_CP017440", "max2_disagreeL"),
    ("klebsiella_pneumoniae", "NZ_CP018427-NZ_CP018428", "max_agreeF"),
    ("klebsiella_pneumoniae", "NZ_CP013322-NZ_CP015382", "min12_agreeF"),
    ("klebsiella_pneumoniae", "NZ_CP010361-NZ_CP018428", "max_disagreeL"),
    ("klebsiella_pneumoniae", "NZ_CP010361-NZ_CP011976", "max2_disagreeL"),
]


rule MG_proj_all:
    # Rule to obtain all comparison plots
    input:
        expand(
            "projections/benchmark/{kind}_{species}.pdf",
            kind=["minimap20-std", "minimap20-noenergy"],
            species=PX_species,
        ),
        [
            f"projections/benchmark/examples/{sp}__minimap20-std__{st}-{t}.pdf"
            for sp, st, t in MG_example_list
        ],
