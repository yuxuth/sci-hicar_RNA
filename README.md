# snakemake_RNA_feature_count

snakemake -p -j 24 --cluster-config cluster.json --cluster "sbatch -J {cluster.job} --mem={cluster.mem} -N 1 -n {threads} -o {cluster.out} -e {cluster.err} " &> log &


add additional reads trim step for the Hi-CAR RNA-seq fastq file
