
#  BUSCO with singularity

Run BUSCO for nucleotide fasta files using already downloaded BUSCO lineages

##  Dependencies

- singularity
- nextflow 


Install singularity
```
conda install -c conda-forge singularity=3.6.1 -y
```

Download BUSCO image
```
singularity pull docker://ezlabgva/busco:4.1.2_cv1
```

## Configure BUSCO version and predownload dataset

Whichever version you choose, you must specify the full path of the SIF file in the configuration file (busco_nf.config, search for "container = a/path/busco_vx.sif"). You can get images of different BUSCO versions at [this link](https://hub.docker.com/r/ezlabgva/busco/tags).

The pipeline also expects the lineages (groupX_odb10) to be already downloaded. To set them up, I open an instance of the BUSCO image and run it with the desired target taxonomic group (the file used for -i can be anything as all we want is the automatic download and decompression). Stop the execution after the dataset decompression. 

```
singularity shell busco_4.1.2_cv1.sif
# Inside this image:
busco -i anyFile.txt -l groupX_odb10 --out tmp -f -m geno
#ctrl+C
```


##  Parameters
`--busco_downloads` indicates the directory that busco created to download the reference datasets.
`--genomes`  glob path that captures the assemblies you want to assess. The files should have .fasta or .fasta.gz extension.
`--outdir` directory where you want to save the results
`--odb` comma separated string of datasets you want to assess on each of the fasta files.
`-profile` configuration specific to the machine where you are running the pipeline [default or farm].


## Output
In the output directory you will find for each of your fasta files, the full and short tables together with the single and multicopy sequences assessed per reference set.

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


