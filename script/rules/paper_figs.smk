
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

rule PF_accuracy_single_kernel_plot:
    message:
        "performing accuracy plot for kernel {wildcards.kernel}"
    input:
        rules.PF_extract_accuracy.output,
    output:
        "figs/paper/accuracy_{kernel}.pdf",
    params:
        snps=AC_snps_accplot,
        mut_factor=AC_config["mut-factor"],
    conda:
        "../conda_envs/bioinfo_env.yml"
    shell:
        """
        python3 workflow_scripts/paper_figs/accuracy_single_kernel.py \
            --json {input} \
            --ker {wildcards.kernel} --snps {params.snps} --mut_factor {params.mut_factor} \
            --pdf {output}
        """

# rule PF_accuracy_kernel_comparison_plot:
#     message:
#         "performing accuracy plot for kernel {wildcards.kernel}"
#     input:
#         mm10="synthetic_data/results/accuracy-minimap10.csv",
#         mm20="synthetic_data/results/accuracy-minimap20.csv",
#         mmsq="synthetic_data/results/accuracy-mmseqs.csv",
#     output:
#         "figs/paper/accuracy_kernel_comparison.pdf",
#     shell:
#         """
#         python3 workflow_scripts/paper_figs/accuracy_kernel_comparison.py \
#             --mm10 {input.mm10} -mm20 {input.mm20} --mmsq {input.mmsq} \
#             --pdf {output}
#         """

rule PF_all:
    input:
        rules.PF_benchmark_synthetic.output,
        expand("figs/paper/accuracy_{kernel}.pdf", kernel=AC_ker_names),
        # rules.PF_accuracy_kernel_comparison_plot.output,
