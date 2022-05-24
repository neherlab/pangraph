Example of slurm launch command:

```bash
snakemake all --jobs=20 \
    --cluster "sbatch --time={cluster.time} --mem={cluster.mem} --cpus-per-task={cluster.n} --qos={cluster.qos}" \
    --jobscript cluster/slurm_submit.sh \
    --cluster-config cluster/cluster_config.json \
    --jobname "{rulename}_{jobid}" \
    --latency-wait 60
    --use-conda
```