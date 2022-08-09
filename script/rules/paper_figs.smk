
rule PF_benchmark_synthetic:
    message:
        "Creating paper plots for time/memory performances on synthetic data."
    input:
        rules.SB_summary_dataframe.output,
    output:
        v1="figs/paper/benchmark_synthetic_data_v1.pdf",
        v2="figs/paper/benchmark_synthetic_data_v2.pdf",
    conda:
        "../conda_envs/bioinfo_env.yml"
    shell:
        """
        python3 workflow_scripts/paper_figs/synthetic_data_benchmark.py \
            --csv {input} --pdf1 {output.v1} --pdf2 {output.v2}
        """


rule PF_extract_accuracy:
    message:
        "Extracting JLD2 file to csv table for plotting - kernel = {wildcards.kernel}"
    input:
        rules.AC_accuracy_database.output,
    output:
        "synthetic_data/results/accuracy-{kernel}.json",
    shell:
        """
        julia -t 1 --project=. workflow_scripts/paper_figs/accuracy_to_json.jl \
            {input} {output}
        """


rule PF_accuracy_plots:
    message:
        "Creating accuracy plots"
    input:
        mm10="synthetic_data/results/accuracy-minimap10.json",
        mm20="synthetic_data/results/accuracy-minimap20.json",
        mmsq="synthetic_data/results/accuracy-mmseqs.json",
    output:
        expand("figs/paper/accuracy/accuracy_{kernel}.pdf", kernel=AC_ker_names),
        "figs/paper/accuracy/accuracy_comparison.pdf",
        "figs/paper/accuracy/snps_rate_vs_divergence.pdf"
    params:
        snps=AC_snps_accplot,
        pdf_fld="figs/paper/accuracy"
    conda:
        "../conda_envs/bioinfo_env.yml"
    shell:
        """
        python3 workflow_scripts/paper_figs/accuracy_plots.py \
            --mm10 {input.mm10} --mm20 {input.mm20} --mmsq {input.mmsq} \
            --snps {params.snps} --pdf_fld {params.pdf_fld}
        """


rule PF_all:
    input:
        rules.PF_benchmark_synthetic.output,
        rules.PF_accuracy_plots.output,
