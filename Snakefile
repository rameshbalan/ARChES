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
number_of_species = len(config["species"])

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
        # "results/annotated_blast_file.txt"
        # expand("run_transdecoder_{species}_nr95/short_summary_transdecoder_{species}_nr95.txt", species=config["species"]),

rule denovo:
    input:
        "data/{species}_samples_trinity.txt"
    output:
        "results/{species}/trinity_out_dir/Trinity.fasta"
    params:
        sample_name="{species}"
    log:
        "logs/denovo_{species}_log.txt"
    conda:
        "envs/trinity_v2.11.0.yml"
    shell:
        "(Trinity --seqType fq \
          --samples_file {input} \
          --max_memory 8G --CPU {threads} \
          --output results/{params.sample_name}/trinity_out_dir) 2> {log}"

rule rename_headers:
    input:
        "results/{species}/trinity_out_dir/Trinity.fasta"
    output:
        "results/{species}/trinity_out_dir/{species}.fasta"
    params:
        sample_name="{species}"
    conda:
        "envs/busco_v3.0.2.yml"
    log:
        "logs/rename_headers_{species}_log.txt"
    shell:
        """
        (sed s/TRINITY/{params.sample_name}/g {input} > {output}) 2> {log}
        """

rule ref_prep:
    input:
        expand("results/{species}/trinity_out_dir/{species}.fasta", species=config["species"])
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
        fasta = "results/{species}/trinity_out_dir/{species}.fasta",
        busco_lineage = "data/ref_odb"
    output:
        "run_busco_trinity_{species}/short_summary_busco_trinity_{species}.txt"
    params:
        sample_name="busco_trinity_{species}"
    conda:
        "envs/busco_v3.0.2.yml"
    log:
        "logs/trinity_busco_{species}_log.txt"
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
        "results/{species}/trinity_out_dir/{species}.fasta"
    output:
        "{species}_transdecoder/longest_orfs.pep"
    params:
        "{species}_transdecoder"
    conda:
        "envs/transdecoder_v5.5.0.yml"
    log:
        "logs/transdecoder_longorfs_{species}_log.txt"
    shell:
        "(TransDecoder.LongOrfs -t {input} -O {params}) 2> {log}"

rule AddHomologyBLASTp:
    input:
        pep="{species}_transdecoder/longest_orfs.pep",
        ref_prot="data/ref_protein.fasta"
    output:
        "results/blastp_{species}.outfmt6"
    conda:
        "envs/busco_v3.0.2.yml"
    log:
        "logs/AddHomologyBLASTp_{species}_log.txt"
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
        fasta="results/{species}/trinity_out_dir/{species}.fasta",
        blastp="results/blastp_{species}.outfmt6"
    output:
        "{species}.fasta.transdecoder.cds",
        "{species}.fasta.transdecoder.pep"
    params:
        "{species}_transdecoder"
    conda:
        "envs/transdecoder_v5.5.0.yml"
    log:
        "logs/transdecoder_predict_aa_{species}_log.txt"
    shell:
        """
        (TransDecoder.Predict -t {input.fasta} \
                             -O {params} \
                             --retain_blastp_hits {input.blastp} \
                             --single_best_only) 2> {log}
        """

rule cluster_cdhitest:
    input:
        "{species}.fasta.transdecoder.cds"
    output:
        "results/{species}_nr95.fasta"
    params:
        config["cd_hit_threshold"]
    conda:
        "envs/cd-hit_v4.8.1.yml"
    log:
        "logs/cluster_cdhitest_{species}_log.txt"
    shell:
        "(cd-hit-est -i {input} -o {output} -c {params}) 2> {log}"

rule nr95_transdecoder_longorfs:
    input:
        "results/{species}_nr95.fasta"
    output:
        "{species}_nr95_transdecoder/longest_orfs.pep"
    params:
        "{species}_nr95_transdecoder"
    conda:
        "envs/transdecoder_v5.5.0.yml"
    log:
        "logs/nr95_transdecoder_longorfs_{species}_log.txt"
    shell:
        "(TransDecoder.LongOrfs -t {input} -O {params}) 2> {log}"

rule nr95_transdecoder_predict_aa:
    input:
        fasta="results/{species}_nr95.fasta",
        pep="{species}_nr95_transdecoder/longest_orfs.pep"
    output:
        "{species}_nr95.fasta.transdecoder.cds",
        "{species}_nr95.fasta.transdecoder.pep"
    params:
        "{species}_nr95_transdecoder"
    conda:
        "envs/transdecoder_v5.5.0.yml"
    log:
        "logs/nr95_transdecoder_predict_aa_{species}_log.txt"
    shell:
        """
        (TransDecoder.Predict -t {input.fasta} \
                              -O {params} \
                              --single_best_only) 2> {log}
        """

rule nr95_clean_transdecoder_predict_aa:
    input:
        cds="{species}_nr95.fasta.transdecoder.cds",
        pep="{species}_nr95.fasta.transdecoder.pep"
    output:
        "results/{species}_nr95.fasta.transdecoder.cds",
        "results/proteomes/{species}_nr95_pep.fasta"
    params:
        "results/proteomes/{species}_nr95_pep.fasta"
    conda:
        "envs/transdecoder_v5.5.0.yml"
    log:
        "logs/nr95_clean_transdecoder_predict_aa_{species}_log.txt"
    shell:
        """
        (cp {input.cds} results) 2> {log}
        (cp {input.pep} {params}) 2> {log}
        """

rule transdecoder_busco:
    input:
        trinity=expand("run_busco_trinity_{species}/short_summary_busco_trinity_{species}.txt",species=config["species"]),
        fasta="results/{species}_nr95.fasta.transdecoder.cds",
        busco_lineage = "data/ref_odb"
    output:
        "run_transdecoder_{species}_nr95/short_summary_transdecoder_{species}_nr95.txt"
    params:
        "transdecoder_{species}_nr95"
    conda:
        "envs/busco_v3.0.2.yml"
    log:
        "logs/transdecoder_busco_{species}_log.txt"
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
        fasta="results/{species}_nr95.fasta.transdecoder.cds"
    output:
        "results/{species}_nr95_index/versionInfo.json"
    params:
        "results/{species}_nr95_index"
    conda:
        "envs/salmon_v1.2.1.yml"
    log:
        "logs/transcriptome_index_{species}_log.txt"
    shell:
        "(salmon index -t {input.fasta} -i {params}) 2> {log}"

rule quantify_salmon:
    input:
        index_check=expand("results/{species}_nr95_index/versionInfo.json",species=config["species"]),
        read1="data/{sample}_R1.fq.gz",
        read2="data/{sample}_R2.fq.gz"
    output:
        "results/quants/{sample}/quant.sf"
    params:
        samplename="{sample}",
        index=lambda wildcards: \
                    "results/"+config["samples_species"][wildcards.sample]+"_nr95_index"
    conda:
        "envs/salmon_v1.2.1.yml"
    log:
        "logs/quantify_salmon_{sample}_log.txt"
    shell:
        "(salmon quant -i {params.index} -l IU \
         -1 {input.read1} \
         -2 {input.read2} \
         -p {threads} \
         --validateMappings \
         -o results/quants/{params.samplename}) 2> {log}"

rule orthologs_orthofinder:
    input:
        expand("results/proteomes/{species}_nr95_pep.fasta",species=config["species"])
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
        species_count=str(number_of_species),
        neox=config["annotate_orthofinder"]["neox"],
        x=config["annotate_orthofinder"]["x"],
        un=config["annotate_orthofinder"]["un"]
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
        (bash scripts/annotate_all_OG.sh {input.ref_gff} {output} {params.neox} {params.x} {params.un}) 2> {log}
        """


rule TMM_and_annotate:
    input:
        sample="data/{species}_samples.txt",
        quant_files=expand("results/quants/{sample}/quant.sf",sample=config["samples"]),
        blast="results/annotated_blast_file.txt"
    output:
        TMM="results/{species}_TMM.txt",
        ann_TMM="results/{species}_TMM_ann.txt"
    conda:
        "envs/deseq2_v1.28_edgeR_v3.30.yml"
    log:
        "logs/TMM_and_annotate_{species}_log.txt"
    shell:
        """
        (Rscript scripts/rawcounts_and_TMM.R -f {input.sample} -o {output.TMM}) 2> {log}
        (python3 scripts/add_chrom_OG_to_quant.py {input.blast} {output.TMM} {output.ann_TMM}) 2> {log}
        """

rule plot_busco:
    input:
        transdecoder_busco=expand("run_transdecoder_{species}_nr95/short_summary_transdecoder_{species}_nr95.txt",species=config["species"])
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

rule summary_table:
    input:
        expand("{species}_nr95.fasta.transdecoder.cds",species=config["species"])
    params:
        species_name = expand("{species}",species=config["species"])
    output:
        summary_table="results/summary_table.csv",
        species=expand("{species}.txt",species=config["species"])
    log:
        "logs/summary_table_log.txt"
    shell:
        """
        bash scripts/summary_table.sh {params.species_name}
        """

rule report:
    input:
        ann_TMM=expand("results/{species}_TMM_ann.txt",species=config["species"]),
        busco_fig='results/busco/busco_figure.png',
        summary_table=os.path.abspath("results/summary_table.csv")
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
        # mv *transdecoder* tmp_dir
        # mv *.cmds tmp_dir
        # mv run* tmp_dir
        # mv *.txt tmp_dir
        (Rscript -e \
        "rmarkdown::render('scripts/generate_report.Rmd',\
        params = list(neo_species='{params.neo_species}', \
        busco_fig='{params.busco_fig}',\
        summary_table='{input.summary_table}',\
        tree='{params.tree}'), \
        output_file='../report.html')") 2> {log}
        """
