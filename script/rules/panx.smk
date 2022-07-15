# This snakemake workflow contains recipes to analyze the panx dataset with pangraph,
# and benchmark pangraph performances on real data.

import glob


def acc_nums(fld):
    files = glob.glob(fld + "/*.gbk")
    p = re.compile("/?([^/]+)\.gbk$")
    return [re.search(p, f).groups()[0] for f in files]


SPECIES = [
    "klebsiella_pneumoniae",
    "helicobacter_pylori",
    "prochlorococcus_marinus",
    "mycobacterium_tuberculosis",
    "escherichia_coli",
]
ACCNUMS = {s: acc_nums(f"panx_data/{s}/input_GenBank") for s in SPECIES}

ker_opt = {
    "minimap10-full": "-k minimap2 -s 10",
    "minimap20-full": "-k minimap2 -s 20",
    "minimap-noenergy": "-k minimap2 -s 20 -a 0 -b 0",
    "mmseqs-full": "-k mmseqs -b 10",
    "mmseqs-noenergy": "-k mmseqs -a 0 -b 0",
}

KER_KINDS = list(ker_opt.keys())


wildcard_constraints:
    species=f"({'|'.join(SPECIES)})",
    kind=f"({'|'.join(KER_KINDS)})",


rule all:
    input:
        expand(
            "panx_data/{species}/pangraphs/pangraph-{kind}.json",
            species=SPECIES,
            kind=KER_KINDS,
        ),
        "panx_data/benchmark/benchmark_compression.csv",
        "panx_data/benchmark/benchmark_summary.csv",


rule gbk_to_fa:
    message:
        "converting genbank to fasta for {wildcards.species} - {wildcards.acc} "
    input:
        "panx_data/{species}/input_GenBank/{acc}.gbk",
    output:
        "panx_data/{species}/fa/{acc}.fa",
    conda:
        "cluster/bioinfo_env.yml"
    shell:
        """
        python3 workflow_scripts/gbk_to_fa.py --gbk {input} --fa {output}
        """


# ------------- Test on real PanX data -----------


# for correct benchmarking needs to have pangraph binary present in the path.
rule build_full_pangraph:
    message:
        "building full pangraph ({wildcards.kind}) for {wildcards.species}"
    input:
        lambda w: expand("panx_data/{{species}}/fa/{acc}.fa", acc=ACCNUMS[w.species]),
    output:
        pg="panx_data/{species}/pangraphs/pangraph-{kind}.json",
        bm="panx_data/{species}/pangraphs/benchmark/pangraph-{kind}.txt",
    params:
        opt=lambda w: ker_opt[w.kind],
    conda:
        "cluster/pangraph_build_env.yml"
    shell:
        """
        echo "species = {wildcards.species}" > {output.bm}
        echo "kind = {wildcards.kind}" >> {output.bm}
        export JULIA_NUM_THREADS=8
        /usr/bin/time --verbose -o {output.bm} -a pangraph build --circular {params.opt} {input} > {output.pg}
        """


rule summary_performance_benchmark:
    message:
        "Summary of pangraph performances"
    input:
        expand(
            "panx_data/{species}/pangraphs/benchmark/pangraph-{kind}.txt",
            species=SPECIES,
            kind=KER_KINDS,
        ),
    output:
        csv="panx_data/benchmark/benchmark_summary.csv",
        pdf="panx_data/benchmark/benchmark_summary.pdf",
    conda:
        "cluster/bioinfo_env.yml"
    shell:
        """
        python3 workflow_scripts/summary_benchmark.py {output.csv} {output.pdf} {input}
        """


rule compression_benchmark:
    message:
        "Compression performances for species {wildcards.species}"
    input:
        pang=expand(
            "panx_data/{{species}}/pangraphs/pangraph-{kind}.json", kind=KER_KINDS
        ),
        fa=lambda w: expand(
            "panx_data/{{species}}/fa/{acc}.fa", acc=ACCNUMS[w.species]
        ),
    output:
        json="panx_data/benchmark/{species}/compression.json",
    conda:
        "cluster/bioinfo_env.yml"
    shell:
        """
        python3 workflow_scripts/compression_benchmark.py --fasta {input.fa} --pangraphs {input.pang} --out_json {output.json}
        """


rule summary_compression_benchmark:
    message:
        "Summary of compression performances"
    input:
        expand("panx_data/benchmark/{species}/compression.json", species=SPECIES),
    output:
        csv="panx_data/benchmark/benchmark_compression.csv",
        pdf="panx_data/benchmark/benchmark_compression.pdf",
    conda:
        "cluster/bioinfo_env.yml"
    shell:
        """
        python3 workflow_scripts/compression_summary.py --jsons {input} --csv {output.csv} --pdf {output.pdf}
        """


# ------------- Test pairwise graphs vs graph merging -----------


# rule projection_full_graph:
#     message:
#         "Creating projection graph for species {wildcards.species}"
#     input:
#         lambda w: expand("panx_data/{{species}}/fa/{acc}.fa", acc=ACC50[w.species]),
#     output:
#         "projections/{species}/full/pangraph_{kind}.json",
#     params:
#         opt=lambda w: ker_opt[w.kind],
#     conda:
#         "cluster/pangraph_build_env.yml"
#     shell:
#         """
#         export JULIA_NUM_THREADS=8
#         pangraph build --circular {params.opt} {input} > {output}
#         """

# rule pairwise_graphs:
#     message:
#         "Building pairwise graph for strains {wildcards.s1} - {wildcards.s2} ({wildcards.species})"
#     input:
#         lambda w: expand("panx_data/{{species}}/fa/{acc}.fa", acc=ACC50[w.species]),
#     output:
#         "projections/{species}/pairwise/pangraph_{kind}_{s1}|{s2}.json",
#     params:
#         opt=lambda w: ker_opt[w.kind],
#     conda:
#         "cluster/pangraph_build_env.yml"
#     shell:
#         """
#         pangraph build --circular {params.opt} {input} > {output}
#         """

# rule pairwise_projection:
#     message:
#         "project graph on strains {wildcards.s1} - {wildcards.s2} ({wildcards.species})"
#     input:
#         "projections/{species}/full/pangraph_{kind}.json"
#     output:
#         "projections/{species}/projected/pangraph_{kind}_{s1}|{s2}.json"
#     conda:
#         "cluster/pangraph_build_env.yml"
#     shell:
#         """
#         pangraph marginalize -s {wildcards.s1},{wildcards.s2} > {output}
#         """

# rule compare_projection_pairwise:
#     message:
#         "comparing projected graph to pairwise graph ({wildcards.s1} - {wildcards.s2} ; {wildcards.species})"
#     input:
#         pw=""
#         pj=""
#     output:
#         "projections/{species}/comparison/{s1}|{s2}.json"
#     shell:
#         ""
