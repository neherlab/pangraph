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


rule PX_diversity:
    message:
        "Evaluating dataset diversity for species {wildcards.species}"
    input:
        gc="panx_data/{species}/vis/geneCluster.json",
        fld="panx_data/{species}/vis/geneCluster",
    output:
        "panx_diversity/{species}.json",
    params:
        ntot=lambda w: len(PX_accnums[w.species]),
    conda:
        "../conda_envs/bioinfo_env.yml"
    shell:
        """
        python3 workflow_scripts/panx_dataset_diversity.py \
            --gc {input.gc} --species {wildcards.species} \
            --ntot {params.ntot} --gcfolder {input.fld} --out {output}
        """


rule PX_diversity_all:
    input:
        expand(rules.PX_diversity.output, species=PX_species),
    output:
        csv="panx_diversity/all.csv",
    run:
        import pandas as pd
        import json

        data = []
        for i in input:
            with open(i, "r") as f:
                data.append(json.load(f))
        pd.DataFrame(data).to_csv(output.csv, index=False)


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
        fasta=expand("panx_data/escherichia_coli/fa/{acc}.fa", acc=PX_IS_allstrains,),
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
            rules.PX_IS_extract_stats.output, size=PX_IS_sizes, trial=PX_IS_trials,
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
