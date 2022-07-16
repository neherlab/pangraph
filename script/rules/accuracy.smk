import json

# extract config for accuracy rules
AC_config = config["accuracy"]

# simulation parameters
with open(AC_config["sim-params"], "r") as f:
    sim_params = json.load(f)
HGT = sim_params["hgt"]
SNPS = sim_params["snps"]
SNPS_accplot = sim_params["snps-accplot"]

# kernels and number of trials
AC_ker_opt = AC_config["kernel-options"]
AC_ker_names = list(AC_ker_opt.keys())
AC_trials = list(range(1, AC_config["sim-ntrials"] + 1))

# pangraph project folder
pgf = "./.."


rule AC_all:
    input:
        expand("figs/paper-accuracy-{kernel}.png", kernel=AC_ker_names),
        "figs/paper-accuracycomp.pdf",


rule AC_generate_data:
    message:
        "generating pangraph with hgt = {wildcards.hgt}, snps = {wildcards.snps}, n = {wildcards.n}"
    output:
        graph="synthetic_data/generated/{hgt}_{snps}/known_{n}.json",
        seqs="synthetic_data/generated/{hgt}_{snps}/seqs_{n}.fa",
    params:
        N=100,
        T=50,
        L=50000,
        pgf=pgf,
    shell:
        """
        julia -t 1 --project=. make-sequence.jl -N {params.N} -L {params.L} \
        | julia -t 1 --project={params.pgf} {params.pgf}/src/PanGraph.jl generate \
            -m {wildcards.snps} -r {wildcards.hgt} -t {params.T} -i "1e-2" \
            -o {output.graph} > {output.seqs}
        """


ruleorder: AC_guess_pangraph_mmseqs > AC_guess_pangraph


rule AC_guess_pangraph:
    message:
        """
        reconstructing pangraph with kernel {wildcards.kernel}
        hgt = {wildcards.hgt}, snps = {wildcards.snps}, n = {wildcards.n}
        """
    input:
        rules.AC_generate_data.output.seqs,
    output:
        "synthetic_data/{kernel}/{hgt}_{snps}/guess_{n}.json",
    params:
        ker=lambda w: AC_ker_opt[w.kernel],
        pgf=pgf,
    shell:
        """
        julia -t 1 --project={params.pgf} {params.pgf}/src/PanGraph.jl build \
            --circular -a 0 -b 0 {params.ker} {input} > {output}
        """


rule AC_guess_pangraph_mmseqs:
    message:
        """
        reconstructing pangraph with kernel mmseqs
        hgt = {wildcards.hgt}, snps = {wildcards.snps}, n = {wildcards.n}
        """
    input:
        rules.AC_generate_data.output.seqs,
    output:
        "synthetic_data/mmseqs/{hgt}_{snps}/guess_{n}.json",
    params:
        ker=AC_ker_opt["mmseqs"],
        pgf=pgf,
    conda:
        "../conda_envs/pangraph_build_env.yml"
    shell:
        """
        julia -t 8 --project={params.pgf} {params.pgf}/src/PanGraph.jl build \
            --circular -a 0 -b 0 {params.ker} {input} > {output}
        """


rule AC_single_accuracy:
    message:
        """
        generating partial accuracy database for:
        kernel = {wildcards.kernel} hgt = {wildcards.hgt}, snps = {wildcards.snps}
        """
    input:
        known=expand(
            "synthetic_data/generated/{{hgt}}_{{snps}}/known_{n}.json", n=AC_trials
        ),
        guess=expand(
            "synthetic_data/{{kernel}}/{{hgt}}_{{snps}}/guess_{n}.json", n=AC_trials
        ),
    output:
        temp("synthetic_data/{kernel}/{hgt}_{snps}/partial_accuracy.jld2"),
    shell:
        """
        julia -t 1 --project=. make-accuracy.jl {output} {input}
        """


rule AC_accuracy_database:
    message:
        "generating accuracy database for kernel {wildcards.kernel}"
    input:
        expand(
            "synthetic_data/{{kernel}}/{hgt}_{snps}/partial_accuracy.jld2",
            hgt=HGT,
            snps=SNPS,
        ),
    output:
        "synthetic_data/results/accuracy-{kernel}.jld2",
    shell:
        """
        julia -t 1 --project=. concatenate-database.jl {output} {input}
        """


rule AC_accuracy_plots:
    message:
        "generating accuracy plot for kernel {wildcards.kernel}"
    input:
        rules.AC_accuracy_database.output,
    output:
        "figs/cdf-accuracy-{kernel}.png",
        "figs/heatmap-accuracy-{kernel}.png",
        "figs/paper-accuracy-{kernel}.png",
        "figs/paper-accuracy-{kernel}.pdf",
    params:
        snps=SNPS_accplot,
    shell:
        """
        julia -t 1 --project=. plot-accuracy.jl {input} figs {params.snps}
        """


rule AC_accuracy_comparison_plots:
    message:
        "generating accuracy comparison plots"
    input:
        expand("synthetic_data/results/accuracy-{kernel}.jld2", kernel=AC_ker_names),
    output:
        "figs/paper-accuracycomp.pdf",
        "figs/paper-accuracycomp-mutdens.pdf",
        "figs/paper-accuracycomp-scatter.pdf",
    shell:
        """
        julia -t 1 --project=. plot-accuracy-comparison.jl figs {input}
        """
