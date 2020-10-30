#  BUSCO with singularity

Run BUSCO for nucleotide fasta files using already downloaded BUSCO lineages

Install singularity
```
conda install -c conda-forge singularity=3.6.1 -y
```

Download BUSCO image
```
singularity pull docker://ezlabgva/busco:v4.1.4_cv1
```

## Executing the pipeline in the farm

Submit an interactive job. Choose a queue appropriate to the time of execution of your pipeline. Take into account that the default config submits jobs to the 'normal' queue of 12 hours execution and if the job fails it will be resubmitted into the 'long' queue.
```
mbMem=5000; bsub -n 1 -q long -R"span[hosts=1] select[mem>${mbMem}] rusage[mem=${mbMem}]" -M${mbMem} -Is bash
```

Launch the pipeline
```
nextflow -c /lustre/scratch116/tol/projects/tol-nemotodes/sw/nxf_pipelines/busco_nf.config run /lustre/scratch116/tol/projects/tol-nemotodes/sw/nxf_pipelines/busco.nf
		 --busco_downloads /lustre/scratch116/tol/teams/team301/dbs/busco_2020_08/busco_downloads/
		 --genomes '/lustre/scratch116/tol/teams/team301/dbs/fasta_genomes/insects/*fasta.gz'
		 --outdir /lustre/scratch116/tol/teams/team301/users/pg17/busco_insects_renamed/
		 --odb insecta_odb10,endopterygota_odb10
		 -profile farm
```


