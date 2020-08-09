##############################################
# To check the pipeline
# snakemake --use-conda --use-singularity --cores 16 -np
##############################################
# To run the pipeline
# snakemake --use-conda --use-singularity --cores 16
##############################################

configfile: "config.yml"

import os
import glob

basedir = os.path.dirname("Snakefile")
number_of_species = len(config["samples"])

if not os.path.exists('results') == True:
    os.makedirs(os.path.join(basedir, "results"))
if not os.path.exists('results/proteomes') == True:
    os.makedirs(os.path.join(basedir, "results","proteomes"))
if not os.path.exists('tmp_dir') == True:
    os.makedirs(os.path.join(basedir, "tmp_dir"))
if not os.path.exists(os.path.join(basedir, "results", "busco")) == True:
    os.makedirs(os.path.join(basedir, "results", "busco"))
if not os.path.exists('logs') == True:
    os.makedirs(os.path.join(basedir, "logs"))
if not os.path.exists('data/ref_odb') == True:
    os.makedirs(os.path.join(basedir, "data","ref_odb"))


rule all:
    input:
        "report.html"

rule denovo:
    input:
        read1_M1="data/{sample}_M1_R1.fq.gz",
        read1_M2="data/{sample}_M2_R1.fq.gz",
        read1_F1="data/{sample}_F1_R1.fq.gz",
        read1_F2="data/{sample}_F2_R1.fq.gz",
        read2_M1="data/{sample}_M1_R2.fq.gz",
        read2_M2="data/{sample}_M2_R2.fq.gz",
        read2_F1="data/{sample}_F1_R2.fq.gz",
        read2_F2="data/{sample}_F2_R2.fq.gz"
    output:
        "results/{sample}/trinity_out_dir/Trinity.fasta"
    params:
        sample_name="{sample}"
    log:
        "logs/denovo_{sample}_log.txt"
    conda:
        "envs/trinity_v2.11.0.yml"
    shell:
        "(Trinity --seqType fq \
          --left {input.read1_M1},{input.read1_M2},{input.read1_F1},{input.read1_F2}  \
          --right {input.read2_M1},{input.read2_M2},{input.read2_F1},{input.read2_F2} \
          --max_memory 8G --CPU {threads} \
          --output results/{params.sample_name}/trinity_out_dir) 2> {log}"

rule rename_headers:
    input:
        "results/{sample}/trinity_out_dir/Trinity.fasta"
    output:
        "results/{sample}/trinity_out_dir/{sample}.fasta"
    params:
        sample_name="{sample}"
    conda:
        "envs/busco_v3.0.2.yml"
    log:
        "logs/rename_headers_{sample}_log.txt"
    shell:
        "(sed s/TRINITY/{params.sample_name}/g {input} > {output}) 2> {log}"

rule busco_prep:
    input:
        expand("results/{samp}/trinity_out_dir/{samp}.fasta", samp=config["samples"])
    output:
        dir=directory("data/ref_odb"),
        prot="data/ref_protein.fasta",
        gff="data/ref.gff",
        index="data/ref_protein.fasta.phr"
    params:
        orthodb=config["busco_prep"]["ortho_db_url"],
        ref_prot=config["busco_prep"]["ref_prot_url"],
        ref_gff=config["busco_prep"]["ref_gff_url"]
    conda:
        "envs/busco_v3.0.2.yml"
    log:
        "logs/busco_prep_log.txt"
    shell:
        """
        mkdir data/ref_odb
        wget -qO- {params.orthodb} | tar -xvzf - -C data/ref_odb --strip-components 1 2> {log}
        wget -qO- {params.ref_prot} | gunzip - > data/ref_protein.fasta 2> {log}
        wget -qO- {params.ref_gff} | gunzip - > data/ref.gff 2> {log}
        makeblastdb -in data/ref_protein.fasta -dbtype prot 2> {log}
        """

rule trinity_busco:
    input:
        fasta = "results/{sample}/trinity_out_dir/{sample}.fasta",
        busco_lineage = "data/ref_odb"
    output:
        "run_busco_trinity_{sample}/short_summary_busco_trinity_{sample}.txt"
    params:
        sample_name="busco_trinity_{sample}"
    conda:
        "envs/busco_v3.0.2.yml"
    log:
        "logs/trinity_busco_{sample}_log.txt"
    shell:
        """
        (run_BUSCO.py -i {input.fasta} \
                    -o {params.sample_name} \
                    -l {input.busco_lineage} \
                    -m tran \
                    -c {threads} \
                    -f) 2> {log}
        """

rule transdecoder_longorfs:
    input:
        "results/{sample}/trinity_out_dir/{sample}.fasta"
    output:
        "{sample}_transdecoder/longest_orfs.pep"
    params:
        "{sample}_transdecoder"
    conda:
        "envs/transdecoder_v5.5.0.yml"
    log:
        "logs/transdecoder_longorfs_{sample}_log.txt"
    shell:
        "(TransDecoder.LongOrfs -t {input} -O {params}) 2> {log}"

rule AddHomologyBLASTp:
    input:
        pep="{sample}_transdecoder/longest_orfs.pep",
        ref_prot="data/ref_protein.fasta"
    output:
        "results/blastp_{sample}.outfmt6"
    conda:
        "envs/busco_v3.0.2.yml"
    log:
        "logs/AddHomologyBLASTp_{sample}_log.txt"
    shell:
        """
        (blastp -query {input.pep} \
               -db {input.ref_prot} \
               -num_threads {threads} \
               -outfmt 6 \
               -max_target_seqs 1 \
               -out {output}) 2> {log}
        """

rule transdecoder_predict_aa:
    input:
        fasta="results/{sample}/trinity_out_dir/{sample}.fasta",
        blastp="results/blastp_{sample}.outfmt6"
    output:
        "{sample}.fasta.transdecoder.cds",
        "{sample}.fasta.transdecoder.pep"
    params:
        "{sample}_transdecoder"
    conda:
        "envs/transdecoder_v5.5.0.yml"
    log:
        "logs/transdecoder_predict_aa_{sample}_log.txt"
    shell:
        """
        (TransDecoder.Predict -t {input.fasta} \
                             -O {params} \
                             --retain_blastp_hits {input.blastp} \
                             --single_best_only) 2> {log}
        """

rule cluster_cdhitest:
    input:
        "{sample}.fasta.transdecoder.cds"
    output:
        "results/{sample}_nr95.fasta"
    conda:
        "envs/cd-hit_v4.8.1.yml"
    log:
        "logs/cluster_cdhitest_{sample}_log.txt"
    shell:
        "(cd-hit-est -i {input} -o {output} -c 0.95) 2> {log}"

rule nr95_transdecoder_longorfs:
    input:
        "results/{sample}_nr95.fasta"
    output:
        "{sample}_nr95_transdecoder/longest_orfs.pep"
    params:
        "{sample}_nr95_transdecoder"
    conda:
        "envs/transdecoder_v5.5.0.yml"
    log:
        "logs/nr95_transdecoder_longorfs_{sample}_log.txt"
    shell:
        "(TransDecoder.LongOrfs -t {input} -O {params}) 2> {log}"

rule nr95_transdecoder_predict_aa:
    input:
        fasta="results/{sample}_nr95.fasta",
        pep="{sample}_nr95_transdecoder/longest_orfs.pep"
    output:
        "{sample}_nr95.fasta.transdecoder.cds",
        "{sample}_nr95.fasta.transdecoder.pep"
    params:
        "{sample}_nr95_transdecoder"
    conda:
        "envs/transdecoder_v5.5.0.yml"
    log:
        "logs/nr95_transdecoder_predict_aa_{sample}_log.txt"
    shell:
        """
        (TransDecoder.Predict -t {input.fasta} \
                              -O {params} \
                              --single_best_only) 2> {log}
        """

rule nr95_clean_transdecoder_predict_aa:
    input:
        cds="{sample}_nr95.fasta.transdecoder.cds",
        pep="{sample}_nr95.fasta.transdecoder.pep"
    output:
        "results/{sample}_nr95.fasta.transdecoder.cds",
        "results/proteomes/{sample}_nr95_pep.fasta"
    params:
        "results/proteomes/{sample}_nr95_pep.fasta"
    conda:
        "envs/transdecoder_v5.5.0.yml"
    log:
        "logs/nr95_clean_transdecoder_predict_aa_{sample}_log.txt"
    shell:
        """
        (cp {input.cds} results) 2> {log}
        (cp {input.pep} {params}) 2> {log}
        """

rule transdecoder_busco:
    input:
        trinity=expand("run_busco_trinity_{species}/short_summary_busco_trinity_{species}.txt",species=config["samples"]),
        fasta="results/{sample}_nr95.fasta.transdecoder.cds",
        busco_lineage = "data/ref_odb"
    output:
        "run_transdecoder_{sample}_nr95/short_summary_transdecoder_{sample}_nr95.txt"
    params:
        "transdecoder_{sample}_nr95"
    conda:
        "envs/busco_v3.0.2.yml"
    log:
        "logs/transdecoder_busco_{sample}_log.txt"
    shell:
        """
        (run_BUSCO.py -i {input.fasta} \
                    -o {params} \
                    -l {input.busco_lineage}\
                    -m tran \
                    -c {threads} \
                    -f) 2> {log}
        """

rule transcriptome_index:
    input:
        fasta="results/{sample}_nr95.fasta.transdecoder.cds"
    output:
        dir = directory("results/{sample}_nr95")
    conda:
        "envs/salmon_v1.2.1.yml"
    log:
        "logs/transcriptome_index_{sample}_log.txt"
    shell:
        "(salmon index -t {input.fasta} -i {output.dir}) 2> {log}"

rule quantify_salmon:
    input:
        index=rules.transcriptome_index.output.dir,
        read1="data/{sample}_{sex}{rep}_R1.fq.gz",
        read2="data/{sample}_{sex}{rep}_R2.fq.gz"
    output:
        "results/quants/{sample}_{sex}{rep}/quant.sf"
    params:
        "{sample}_{sex}{rep}"
    conda:
        "envs/salmon_v1.2.1.yml"
    log:
        "logs/quantify_salmon_{sample}_{sex}{rep}_log.txt"
    shell:
        "(salmon quant -i {input.index} -l IU \
         -1 {input.read1} \
         -2 {input.read2} \
         -p {threads} \
         --validateMappings \
         -o results/quants/{params}) 2> {log}"

rule orthologs_orthofinder:
    input:
        expand("results/proteomes/{species}_nr95_pep.fasta",species=config["samples"])
    output:
        dir=directory("results/proteomes/OrthoFinder")
    params:
        config["orthologs_orthofinder"]["path"]
    conda:
        "envs/orthofinder_v2.4.0.yml"
    log:
        "logs/orthologs_orthofinder_log.txt"
    shell:
        "(orthofinder -f {params} -t {threads}) 2> {log}"


rule annotate_orthofinder:
    input:
        ortho_dir=rules.orthologs_orthofinder.output.dir,
        ref_gff="data/ref.gff",
        ref_prot="data/ref_protein.fasta"
    params:
        species_count=str(number_of_species)
        # dir=glob.glob("results/proteomes/OrthoFinder/Results_*")[0]
    output:
        "results/annotated_blast_file.txt"
    conda:
        "envs/orthofinder_v2.4.0.yml"
    log:
        "logs/annotate_orthofinder_log.txt"
    shell:
        """
        (python3 scripts/pull_common_orthogroups.py {params.species_count} {input.ortho_dir}) 2> {log}
        (bash scripts/blast_common_orthogroups.sh {input.ref_prot} {input.ortho_dir} {threads}) 2> {log}
        (bash scripts/annotate_all_OG.sh {input.ref_gff} {output}) 2> {log}
        """

rule TMM_and_annotate:
    input:
        sample="data/{sample}_samples.txt",
        quant_files=expand("results/quants/{sample}_{sex}{rep}/quant.sf",sample=config["samples"],sex=config["sex"], rep=config["rep"]),
        blast="results/annotated_blast_file.txt"
    output:
        TMM="results/{sample}_TMM.txt",
        ann_TMM="results/{sample}_TMM_ann.txt"
    conda:
        "envs/deseq2_v1.28_edgeR_v3.30.yml"
    log:
        "logs/TMM_and_annotate_{sample}_log.txt"
    shell:
        """
        (Rscript scripts/rawcounts_and_TMM.R -f {input.sample} -o {output.TMM}) 2> {log}
        (python3 scripts/add_chrom_OG_to_quant.py {input.blast} {output.TMM} {output.ann_TMM}) 2> {log}
        """

rule plot_busco:
    input:
        transdecoder_busco=expand("run_transdecoder_{sample}_nr95/short_summary_transdecoder_{sample}_nr95.txt",sample=config["samples"])
    params:
        "results/busco"
    output:
        "results/busco/busco_figure.png"
    conda:
        "envs/busco_v3.0.2.yml"
    log:
        "logs/plot_busco_log.txt"
    shell:
        """
        cp run_*/*.txt results/busco
        (generate_plot.py -wd {params}) 2> {log}
        """

rule report:
    input:
        ann_TMM=expand("results/{sample}_TMM_ann.txt",sample=config["samples"]),
        busco_fig='results/busco/busco_figure.png'
    params:
        tree=os.path.abspath(config["report"]["tree"]),
        busco_fig=os.path.abspath("results/busco/busco_figure.png"),
        neo_species=config["report"]["neo_species"]
    output:
        "report.html"
    conda:
        "envs/deseq2_v1.28_edgeR_v3.30.yml"
    log:
        "logs/report_log.txt"
    shell:
        """
        mv *transdecoder* tmp_dir
        mv *.cmds tmp_dir
        mv run* tmp_dir
        (Rscript -e \
        "rmarkdown::render('scripts/generate_report.Rmd',\
        params = list(neo_species='{params.neo_species}', \
        busco_fig='{params.busco_fig}',\
        tree='{params.tree}'), \
        output_file='../report.html')") 2> {log}
        """

