
rule PF_benchmark_synthetic:
    message:
        "Creating paper plots for time/memory performances on synthetic data."
    input:
        rules.SB_summary_dataframe.output,
    output:
        main="figs/paper/synth_data_benchmark/main.pdf",
        suppl="figs/paper/synth_data_benchmark/suppl.pdf",
    conda:
        "../conda_envs/bioinfo_env.yml"
    shell:
        """
        python3 workflow_scripts/paper_figs/synthetic_data_benchmark.py \
            --csv {input} --pdf_main {output.main} --pdf_suppl {output.suppl}
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
        suppl_acc="figs/paper/accuracy/accuracy_comparison.pdf",
        suppl_snps="figs/paper/accuracy/snps_rate_vs_divergence.pdf",
        med="figs/paper/accuracy/median_misplacement_vs_divergence.pdf",
        frac="figs/paper/accuracy/misplaced_fraction_vs_divergence.pdf",
    params:
        snps_suppl=AC_sim_params["snps-accplot-suppl"],
        snps_main=AC_sim_params["snps-accplot-main"],
    conda:
        "../conda_envs/bioinfo_env.yml"
    shell:
        """
        python3 workflow_scripts/paper_figs/accuracy_plots.py \
            --mm10 {input.mm10} --mm20 {input.mm20} --mmsq {input.mmsq} \
            --snps_suppl {params.snps_suppl} --snps_main {params.snps_main} \
            --pdf_acc {output.suppl_acc} --pdf_snps {output.suppl_snps} \
            --pdf_med {output.med} --pdf_frac {output.frac}
        """


rule PF_projection_plot_all:
    message:
        "Creating projection plot for paper"
    input:
        expand(
            rules.MG_pairwise_projection_benchmark.output.csv_full,
            kind="minimap20-std",
            species=PX_species,
        ),
    output:
        "figs/paper/projections/proj_fig.pdf",
    conda:
        "../conda_envs/bioinfo_env.yml"
    shell:
        """
        python3 workflow_scripts/paper_figs/projections_plot.py \
            --csv {input} --pdf {output}
        """


rule PF_projection_plot_single:
    message:
        "Creating projection plot for paper"
    input:
        comp="projections/benchmark/minimap20-std_{species}.full.csv",
        kdist=rules.MG_shared_kmers_summary.output,
        pw_div=rules.PX_pairwise_divergence.output,
    output:
        pdf="figs/paper/projections/proj_single_{species}.pdf",
        svg="figs/paper/projections/proj_single_{species}.svg",
    params:
        klen=PX_config["kmer-size"],
    conda:
        "../conda_envs/bioinfo_env.yml"
    shell:
        """
        python3 workflow_scripts/paper_figs/projections_plot_single.py \
            --csv {input.comp} --kmer_dist {input.kdist} --klen {params.klen} \
            --pairwise_div {input.pw_div} --pdf {output.pdf} --svg {output.svg}
        """


rule PF_incremental_size_plot:
    message:
        "Creating plots for incremental graph size analysis"
    input:
        rules.PX_IS_summary_df.output,
    output:
        pdf_all="figs/paper/incr_size/all_stats.pdf",
        pdf_paper="figs/paper/incr_size/summary.pdf",
    conda:
        "../conda_envs/bioinfo_env.yml"
    shell:
        """
        python3 workflow_scripts/paper_figs/incremental_size_plot.py \
            --csv {input} --pdf_all {output.pdf_all} \
            --pdf_paper {output.pdf_paper}
        """


rule PF_panx_compression_plot:
    message:
        "Creating plots for PanX compression performances"
    input:
        comp=rules.PX_summary_compression_benchmark.output.csv,
        summ=rules.PX_summary_performance_benchmark.output.csv,
    output:
        main="figs/paper/panx/main.pdf",
        suppl="figs/paper/panx/suppl.pdf",
    conda:
        "../conda_envs/bioinfo_env.yml"
    shell:
        """
        python3 workflow_scripts/paper_figs/panx_benchmark_plot.py \
            --csv_comp {input.comp} --csv_summ {input.summ} \
            --pdf_main {output.main} --pdf_suppl {output.suppl}
        """


rule PF_all:
    input:
        rules.PF_benchmark_synthetic.output,
        rules.PF_accuracy_plots.output,
        rules.PF_incremental_size_plot.output,
        rules.PF_panx_compression_plot.output,
        rules.PF_projection_plot_all.output,
        expand(rules.PF_projection_plot_single.output, species=PX_species),
