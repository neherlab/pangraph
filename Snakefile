
cluster_name = "cluster_066_n_10"

rule map:
	input:
		f"inc_clusters/{cluster_name}.fasta"
	output:
		f"mapping/{cluster_name}.paf"
	params:
		"-x asm5 -D -c"
	shell:
		"minimap2 {params} {input} {input} > {output}"


rule cluster_sequences:
	input:
		rule.map.output
	output:
		blocks = f"graphs/blocks_{cluster_name}.fasta",
		strings = f"graphs/sequences_{cluster_name}.json"
	shell:
		"""
			python scripts/cluster_paf.py {input} {output}
		"""
