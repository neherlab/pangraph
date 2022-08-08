
rule PF_benchmark_synthetic:
    message:
        "Creating paper plots for time/memory performances on synthetic data."
    input:
        rules.SB_summary_dataframe.output,
    output:
        v1="figs/paper/benchmark_synthetic_data_v1.pdf",
        v2="figs/paper/benchmark_synthetic_data_v2.pdf",
    shell:
        """
        python3 workflow_scripts/paper_figs/synthetic_data_benchmark.py \
            --csv {input} --pdf1 {output.v1} --pdf2 {output.v2}
        """


rule PF_all:
    input:
        rules.PF_benchmark_synthetic.output,
