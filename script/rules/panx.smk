# This snakemake workflow contains recipes to analyze the panx dataset with pangraph,
# and benchmark pangraph performances on real data.

import glob
import json

PX_config = config["panx"]
# list of different species
PX_species = PX_config["species"]
# kernel options
PX_ker_opt = PX_config["kernel-options"]
# different kernel names
PX_ker_names = list(PX_ker_opt.keys())


# function that given a folder containing genbank files returns a list of
# their accession numbers, extracted as the file basename.
def acc_nums(fld):
    files = glob.glob(fld + "/*.gbk")
    p = re.compile("/?([^/]+)\.gbk$")
    return [re.search(p, f).groups()[0] for f in files]


# create dictionary of accession numbers for each species
PX_accnums = {s: acc_nums(f"panx_data/{s}/input_GenBank") for s in PX_species}


wildcard_constraints:
    species=f"({'|'.join(PX_species)})",
    kind=f"({'|'.join(PX_ker_names)})",
    s1="[^_-][^-]+[^_-]",
    s2="[^_-][^-]+[^_-]",


# ------------- Test performances on real PanX data -----------


# final benchmark plots
rule PX_benchmark_all:
    input:
        "panx_data/benchmark/benchmark_compression.csv",
        "panx_data/benchmark/benchmark_summary.csv",


# convert genbank files to fasta
rule PX_gbk_to_fa:
    message:
        "converting genbank to fasta for {wildcards.species} - {wildcards.acc} "
    input:
        "panx_data/{species}/input_GenBank/{acc}.gbk",
    output:
        "panx_data/{species}/fa/{acc}.fa",
    conda:
        "../conda_envs/bioinfo_env.yml"
    shell:
        """
        python3 workflow_scripts/gbk_to_fa.py --gbk {input} --fa {output}
        """


# Build full pangraphs for the species of interest and the specified alignment kernel option.
# Performance are measured and saved in a txt file.
# NB: for correct benchmarking needs to have pangraph binary present in the path.
rule PX_build_full_pangraph:
    message:
        "building full pangraph ({wildcards.kind}) for {wildcards.species}"
    input:
        lambda w: expand("panx_data/{{species}}/fa/{acc}.fa", acc=PX_accnums[w.species]),
    output:
        pg="panx_data/{species}/pangraphs/pangraph-{kind}.json",
        bm="panx_data/{species}/pangraphs/benchmark/pangraph-{kind}.txt",
    params:
        opt=lambda w: PX_ker_opt[w.kind],
    conda:
        "../conda_envs/pangraph_build_env.yml"
    shell:
        """
        echo "species = {wildcards.species}" > {output.bm}
        echo "kind = {wildcards.kind}" >> {output.bm}
        export JULIA_NUM_THREADS=8
        /usr/bin/time --verbose -o {output.bm} -a pangraph build --circular {params.opt} {input} > {output.pg}
        """


# summarizes the performance of pangraph (time, memory, cpu..) on real data from the
# log-file produced by the previous rule.
rule PX_summary_performance_benchmark:
    message:
        "Summary of pangraph performances"
    input:
        expand(
            "panx_data/{species}/pangraphs/benchmark/pangraph-{kind}.txt",
            species=PX_species,
            kind=PX_ker_names,
        ),
    output:
        csv="panx_data/benchmark/benchmark_summary.csv",
        pdf="panx_data/benchmark/benchmark_summary.pdf",
    conda:
        "../conda_envs/bioinfo_env.yml"
    shell:
        """
        python3 workflow_scripts/summary_benchmark.py {output.csv} {output.pdf} {input}
        """


# For any given species this rule summarizes the performance of pangraph concerning
# the characteristic of the pangenome graph (size and number of blocks, fraction of
# core genome, sequence compression...) and saves the result in a json file. This
# includes all tested alignment kernels.
rule PX_compression_benchmark:
    message:
        "Compression performances for species {wildcards.species}"
    input:
        pang=expand(
            "panx_data/{{species}}/pangraphs/pangraph-{kind}.json", kind=PX_ker_names
        ),
        fa=lambda w: expand(
            "panx_data/{{species}}/fa/{acc}.fa", acc=PX_accnums[w.species]
        ),
    output:
        json="panx_data/benchmark/{species}/compression.json",
    conda:
        "../conda_envs/bioinfo_env.yml"
    shell:
        """
        python3 workflow_scripts/compression_benchmark.py \
            --fasta {input.fa} --pangraphs {input.pang} --out_json {output.json}
        """


# Takes the json file produced by the previous rule for all the different species
# and produces a summary plot.
rule PX_summary_compression_benchmark:
    message:
        "Summary of compression performances"
    input:
        expand("panx_data/benchmark/{species}/compression.json", species=PX_species),
    output:
        csv="panx_data/benchmark/benchmark_compression.csv",
        pdf="panx_data/benchmark/benchmark_compression.pdf",
    conda:
        "../conda_envs/bioinfo_env.yml"
    shell:
        """
        python3 workflow_scripts/compression_summary.py \
            --jsons {input} --csv {output.csv} --pdf {output.pdf}
        """


# ------------- Compare pairwise graphs vs pairwise projections from a larger graph -----------

# load dictionary of selected strains & selected pairs for analysis (produced by pick_projection_strain_set.py)
with open(PX_config["projection-data"], "r") as f:
    PX_proj = json.load(f)


# Creates a graph with all of the selected strains.
rule PX_projection_full_graph:
    message:
        "Creating projection graph for species {wildcards.species}"
    input:
        lambda w: expand(
            "panx_data/{{species}}/fa/{acc}.fa", acc=PX_proj[w.species]["strains"]
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
rule PX_pairwise_graphs:
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
rule PX_pairwise_projection:
    message:
        "project graph on strains {wildcards.s1} - {wildcards.s2} ({wildcards.species})"
    input:
        rules.PX_projection_full_graph.output,
    output:
        "projections/{species}/projected/pangraph_{kind}__{s1}-{s2}.json",
    conda:
        "../conda_envs/pangraph_build_env.yml"
    shell:
        """
        pangraph marginalize -s {wildcards.s1},{wildcards.s2} {input} > {output}
        """


# Evaluate pairwise distances between strains using mash
rule PX_mash_triangle:
    message:
        "evaluating mash distance matrix for species {wildcards.species}"
    input:
        lambda w: expand(
            "panx_data/{{species}}/fa/{acc}.fa", acc=PX_proj[w.species]["strains"]
        ),
    output:
        txt="projections/{species}/mash/mash_distance.txt",
        csv="projections/{species}/mash/mash_distance.csv",
    params:
        k=PX_config["mash-kmer-size"],
    conda:
        "../conda_envs/bioinfo_env.yml"
    shell:
        """
        mash triangle -k {params.k} {input} > {output.txt}
        python3 workflow_scripts/mash_triangle_to_csv.py \
            --mash_tri {output.txt} --csv {output.csv}
        """


# Compares the pairwise and the marginalized graph, checking on what fraction of the genome
# the two agree. Results are saved in a json file
rule PX_compare_projection_pairwise:
    message:
        "comparing projected graph to pairwise graph ({wildcards.s1} - {wildcards.s2} ; {wildcards.species})"
    input:
        pw=rules.PX_pairwise_graphs.output,
        pj=rules.PX_pairwise_projection.output,
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
    pairs = PX_proj[w.species]["pairs"]
    in_files = []
    for s1, s2 in pairs:
        in_files.append(f"projections/{w.species}/comparison/{w.kind}__{s1}-{s2}.json")
    return in_files


# For every species and alignment kernel collects all the statistics on the comparison
# between pairwise and projeted graphs, and produces summary plots.
rule PX_pairwise_projection_benchmark:
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


# Produces detailed plots for specific examples of comparison between pairwise
# and projected graphs.
rule PX_explore_proj_examples:
    message:
        "plot inspecting example {wildcards.species} - {wildcards.kind} - {wildcards.s1} - {wildcards.s2} for pairwise comparison"
    input:
        pw=rules.PX_pairwise_graphs.output,
        pj=rules.PX_pairwise_projection.output,
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
PX_example_list = [
    ("escherichia_coli", "NZ_CP014316-NZ_CP018983", "max_agreeF_max_disagreeL"),
    ("escherichia_coli", "NZ_CP012635-NZ_CP017251", "min_agreeF"),
    ("escherichia_coli", "NZ_CP018252-NZ_CP018983", "min2_agreeF"),
    ("escherichia_coli", "NZ_CP017249-NZ_CP017440", "max2_disagreeL"),
    ("klebsiella_pneumoniae", "NZ_CP018427-NZ_CP018428", "max_agreeF"),
    ("klebsiella_pneumoniae", "NZ_CP013322-NZ_CP015382", "min12_agreeF"),
    ("klebsiella_pneumoniae", "NZ_CP010361-NZ_CP018428", "max_disagreeL"),
    ("klebsiella_pneumoniae", "NZ_CP010361-NZ_CP011976", "max2_disagreeL"),
]


# Rule to obtain all comparison plots
rule PX_proj_all:
    input:
        expand(
            "projections/benchmark/{kind}_{species}.pdf",
            kind=["minimap20-std", "minimap20-noenergy"],
            species=PX_species,
        ),
        [
            f"projections/benchmark/examples/{sp}__minimap20-std__{st}-{t}.pdf"
            for sp, st, t in PX_example_list
        ],


# ------------- General pangenome graph statistics vs n. isolates -----------

# nested dictionary {size (str) -> n. trial (str) -> [list of strains]}.
with open(PX_config["incremental-size"], "r") as f:
    PX_IS_strains = json.load(f)

# list of considered sizes and n. trials
PX_IS_sizes = list(sorted(PX_IS_strains.keys()))
PX_IS_trials = list(sorted(PX_IS_strains[PX_IS_sizes[0]].keys()))
PX_IS_allstrains = PX_accnums["escherichia_coli"]
PX_IS_Ntot = len(PX_IS_allstrains)

# input function wildcards -> (size,trial) -> set of input fasta files
def IS_select_strains(w):
    s = w.size
    t = w.trial
    strains = PX_IS_strains[s][t]
    return [f"panx_data/escherichia_coli/fa/{s}.fa" for s in strains]


# rule to build a pangenome graph for the selected set of strains, specified by the size + trial dictionary
rule PX_IS_build_graph:
    message:
        "Build pangenome graph for incremental size analysis - size={wildcards.size} trial={wildcards.trial}"
    input:
        IS_select_strains,
    output:
        "incremental_size/escherichia_coli/{size,[0-9]+}/{trial,[0-9]+}/pangraph.json",
    params:
        opt=PX_ker_opt["minimap20-std"],
    conda:
        "../conda_envs/pangraph_build_env.yml"
    shell:
        """
        export JULIA_NUM_THREADS=8
        pangraph build --circular {params.opt} {input} > {output}
        """


rule PX_IS_extract_stats:
    message:
        "Extracting stats for incremental size analysis - size={wildcards.size} trial={wildcards.trial}"
    input:
        pang=rules.PX_IS_build_graph.output,
        fasta=IS_select_strains,
    output:
        "incremental_size/escherichia_coli/{size,[0-9]+}/{trial,[0-9]+}/stats.json",
    conda:
        "../conda_envs/bioinfo_env.yml"
    shell:
        """
        python3 workflow_scripts/incr_size_extract_stats.py \
            --pangraph {input.pang} --fasta {input.fasta} --json {output}
        """


rule PX_IS_extract_stats_full_graph:
    message:
        "Extracting stats for full pangenome graph"
    input:
        pang="panx_data/escherichia_coli/pangraphs/pangraph-minimap20-std.json",
        fasta=expand(
            "panx_data/escherichia_coli/fa/{acc}.fa",
            acc=PX_IS_allstrains,
        ),
    output:
        stats=f"incremental_size/escherichia_coli/{PX_IS_Ntot}/0/stats.json",
        link=f"incremental_size/escherichia_coli/{PX_IS_Ntot}/0/pangraph.json",
    conda:
        "../conda_envs/bioinfo_env.yml"
    shell:
        """
        ln -s ../../../../{input.pang} {output.link}
        python3 workflow_scripts/incr_size_extract_stats.py \
            --pangraph {output.link} --fasta {input.fasta} --json {output.stats}
        """


rule PX_IS_summary_df:
    message:
        "Building summary dataframe for incremental size analaysis"
    input:
        jsons=expand(
            rules.PX_IS_extract_stats.output,
            size=PX_IS_sizes,
            trial=PX_IS_trials,
        ),
        json_full=rules.PX_IS_extract_stats_full_graph.output.stats,
    output:
        "incremental_size/summary/escherichia_coli_IS_analysis.csv",
    conda:
        "../conda_envs/bioinfo_env.yml"
    shell:
        """
        python3 workflow_scripts/incr_size_summary_df.py \
            --jsons {input.jsons} {input.json_full} --csv {output}
        """


rule PX_IS_all:
    input:
        rules.PX_IS_summary_df.output,
