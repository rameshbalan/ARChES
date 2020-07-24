## **DCP**: Dosage Compensation Analyses Pipeline using Snakemake.

![](https://img.shields.io/maintenance/yes/2020)
[![](https://img.shields.io/github/license/rameshbalan/dcp)](https://github.com/rameshbalan/dcp/blob/master/LICENSE)

Uses:
[![orthofinder](https://img.shields.io/conda/dn/bioconda/orthofinder?label=orthofinder)](https://github.com/davidemms/OrthoFinder)
[![salmon](https://img.shields.io/conda/dn/bioconda/salmon?label=salmon)](https://salmon.readthedocs.io/en/latest/salmon.html)
[![transdecoder](https://img.shields.io/conda/dn/bioconda/transdecoder?label=transdecoder)](https://github.com/TransDecoder/TransDecoder/wiki)
[![cd-hit](https://img.shields.io/conda/dn/bioconda/cd-hit?label=cdhit)](https://github.com/weizhongli/cdhit/wiki)
[![BUSCO](https://img.shields.io/conda/dn/bioconda/busco?label=busco)](https://busco-archive.ezlab.org/v3/)
[![Trinity](https://img.shields.io/docker/v/trinityrnaseq/trinityrnaseq/2.10.0?label=Trinity)](https://hub.docker.com/r/trinityrnaseq/trinityrnaseq/tags)
[![Trinity](https://img.shields.io/docker/image-size/trinityrnaseq/trinityrnaseq/2.10.0?label=Trinity)](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Trinity-in-Docker)
[![deseq2](https://img.shields.io/conda/dn/bioconda/bioconductor-deseq2?label=deseq2)](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
[![edgeR](https://img.shields.io/conda/dn/bioconda/bioconductor-edgeR?label=edgeR)](https://bioconductor.org/packages/release/bioc/html/edgeR.html)

Depends: [![](https://img.shields.io/github/downloads/conda/conda/4.8.3/conda-4.8.3.tar.gz?label=conda)](https://docs.conda.io/en/latest/miniconda.html)
[![snakemake](https://img.shields.io/conda/dn/bioconda/snakemake.svg?label=snakemake)](https://bioconda.github.io/recipes/snakemake/README.html)
[![](https://img.shields.io/github/downloads-pre/hpcng/singularity/v3.5.3/singularity-3.5.3.tar.gz?label=Singularity)](https://sylabs.io/guides/3.5/user-guide/quick_start.html#quick-installation-steps)

### Summary

This repository contains a [snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline that handles all the steps needed for a sex chromosome dosage compensation analyses. Mainly, the workflow

1. Constructs a non-redundant transcriptome assembly (for each species) using [Trinity](https://github.com/trinityrnaseq/trinityrnaseq).
2. Adds annotation by identifying orthologs using [OrthoFinder](https://github.com/davidemms/OrthoFinder).
3. Estimates Differential Expression (DE) among males and females within each species using [salmon](https://github.com/COMBINE-lab/salmon).
	- This step identifies Dosage Balance between males and females.
4. Using a (given) dated phylogenetic tree of all the species in the analysis, (weighted) ancestral-X chromosome expression is compared against the neo-X chromosome.

__Rules Graph__  

![](https://raw.githubusercontent.com/rameshbalan/rameshbalan.github.io/e8e4494c6c4f9f8bcc6b9b5479f9b983f8db89e3/files/DCA_rulegraph.svg)

### Dependencies.

- [conda](https://docs.conda.io/en/latest/miniconda.html)
	- `conda` should be available in the `$PATH`
	- The workflow was tested in python `3.7` and conda version `4.8.3`.
- [snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
	- This is a conda env and can be installed as follows.
	```bash
	# Install mamba to solve all snakemake Dependencies
	conda install -c conda-forge mamba
	# create/install a snakemake environment
	mamba create -c conda-forge -c bioconda -n snakemake snakemake
	```
	- The snakemake workflow was created and tested in `5.19.3`
- [singularity](https://sylabs.io/guides/3.5/user-guide/quick_start.html#quick-installation-steps)
	- singularity is available as a module in TACC and other HPCs.
	 	```bash
		# This command searches for available singularity version in TACC
		module spider singularity
		```
	- Should be available in the `$PATH`
	- The workflow was tested in singularity version `3.5.3`

### Data Input.

User-specific data and reference files can be configured in the `config.yml` file. For a complete analysis, the following files are needed.

1. RNA-Seq reads
2. Reference Protein file URL
3. Reference GFF file URL
4. Reference Orthodb file URL
5. Dated Phylogenetic Tree

### Usage.

```bash
# Activate the snakemake conda environment
conda activate snakemake
# Do a dry run to verify the workflow and the jobs
snakemake --cores 16 --use-singularity --use-conda -np
# Run the pipeline
snakemake --cores 16 --use-singularity --use-conda
```
> Number of cores can be increased to reflect the compute resource available

### Pipeline Output/ Directory Structure:

```bash
├── config.yml - (analysis specific configuration file)
├── data - (reads, tree and other user input data)
├── envs - (environments for different programs)
├── logs - (log file for each job) - created
├── results - (all the major results) - created
├── report.html - (the final output file with interactive graphs) - created
├── scripts - (scripts for processing)
├── Snakefile - (the driver script for the workflow)
└── tmp_dir - (all temporary files and supporting result files) - created
```

### (Expected) Releases.

Current release is highlighted in bold font below and also is tagged in github.  
- **&alpha;**  
	- Modular workflow  
	- Test on sample dataset  
- &beta;  
	- Parallelize BLAST in 2 bash scripts.
	- Snakemake, Python and R Best Practices  
		- Create Rules directory
- Pre-Release  
	- Test using complete beetle dataset
- Release  
	- Test using some other species dataset (_Drosophila_?)

See [Wiki](https://github.com/rameshbalan/dcp/wiki/Dosage-Compensation-Analyses-Pipeline-(DCP)) for further information.

### Citation:  

Ramesh, Balan and Demuth, Jeff. "A General Framework for Dosage Compensation Analyses using Snakemake" (in prep).2020

### References:

1. [Köster, Johannes and Rahmann, Sven. “Snakemake - A scalable bioinformatics workflow engine”. Bioinformatics 2012.](https://academic.oup.com/bioinformatics/article/28/19/2520/290322)
2. [Grabherr MG, Haas BJ, Yassour M, Levin JZ, Thompson DA, Amit I, Adiconis X, Fan L, Raychowdhury R, Zeng Q, Chen Z, Mauceli E, Hacohen N, Gnirke A, Rhind N, di Palma F, Birren BW, Nusbaum C, Lindblad-Toh K, Friedman N, Regev A. Full-length transcriptome assembly from RNA-seq data without a reference genome.](https://www.nature.com/articles/nbt.1883)
3. [Patro, R., Duggal, G., Love, M. I., Irizarry, R. A., & Kingsford, C. (2017). Salmon provides fast and bias-aware quantification of transcript expression. Nature methods, 14(4), 417-419.](https://www.nature.com/articles/nmeth.4197)
4. [Limin Fu, Beifang Niu, Zhengwei Zhu, Sitao Wu and Weizhong Li, CD-HIT: accelerated for clustering the next generation sequencing data. Bioinformatics, (2012), 28 (23): 3150-3152.](https://academic.oup.com/bioinformatics/article/28/23/3150/192160)
5. [Haas, B., & Papanicolaou, A. J. G. S. (2016). TransDecoder (find coding regions within transcripts).](https://github.com/TransDecoder/TransDecoder/wiki)  
6. [Emms, D. M., & Kelly, S. (2019). OrthoFinder: phylogenetic orthology inference for comparative genomics. Genome biology, 20(1), 1-14.](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1832-y)
7. [Julien, P., Brawand, D., Soumillon, M., Necsulea, A., Liechti, A., Schütz, F., ... & Kaessmann, H. (2012). Mechanisms and evolutionary patterns of mammalian and avian dosage compensation. PLoS Biol, 10(5), e1001328.](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1001328)
8. [Schield, D.R., Card, D.C., Hales, N.R., Perry, B.W., Pasquesi, G.M., Blackmon, H., Adams, R.H., Corbin, A.B., Smith, C.F., Ramesh, B. and Demuth, J.P., 2019. The origins and evolution of chromosomes, dosage compensation, and mechanisms underlying venom regulation in snakes. Genome research, 29(4), pp.590-601.](https://genome.cshlp.org/content/29/4/590.full.pdf)

