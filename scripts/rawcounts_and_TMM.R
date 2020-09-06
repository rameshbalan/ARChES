library("DESeq2")
library("tximport")
library("readr")
library("edgeR")
library("reshape2")
library("optparse")
 
option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="species samples file", metavar="character"),
	make_option(c("-o", "--out"), type="character", default=NULL, 
              help="output file name [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$file) | is.null(opt$out)){
  print_help(opt_parser)
  stop("Input and Output files must be specified", call.=FALSE)
}

samples <- read.table(opt$file, header=TRUE)
files <- file.path("results/quants",samples$sample, "quant.sf")
txi <- tximport(files, type="salmon",txOut=TRUE)
ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ gender)
dds <- DESeq(ddsTxi)
species_rawcount <- as.data.frame(counts(dds),optional = TRUE)
colnames(species_rawcount)<- dds@colData@listData$sample
species_rawcount$id <- row.names(species_rawcount)
row.names(species_rawcount) <- NULL
colnames(species_rawcount) <- c(samples[,1],"id")

one_file_for_length <- read.table(files[1], header = T, sep = "\t")
species <- DGEList(counts = species_rawcount[,1:4], group = samples$gender)
species <- calcNormFactors(species, method = "TMM")
species_TMM <- as.data.frame(rpkm(species, gene.length = one_file_for_length$Length, prior.count=1))
species_TMM$id <- species_rawcount$id
melted_species_TMM <- melt(species_TMM, id=c("id"))
colnames(melted_species_TMM) <- c("id", "sample","log2(TMM)")
melted_species_TMM$gender <- samples$gender[match(melted_species_TMM$sample,samples$sample)]
melted_species_TMM$species <- samples$species[match(melted_species_TMM$sample,samples$sample)]
write.table(melted_species_TMM, file = opt$out, sep = "\t", row.names = FALSE)










# python3 gcorn_annotate.py gcorn_TMM.txt
# my_file <- read.table("annotated_TMM.txt")
#     colnames(my_file) <- c("gene","sample","log2(TMM)","chrom_num","chrom_name","orthogroup","gender","species")
#     my_file_neox_males <- my_file[my_file$chrom_name == "Neo-X" & my_file$gender == "Male",]
#     my_file_neox_males$category <- ifelse(my_file_neox_males$species == "T_conf", "Neo-X", "Proto-X")

#     ## Scaling Neo-X in males for each species and combining the dataframes
#     x <- my_file[my_file$chrom_name != "Un" & my_file$chrom_name != "Unknown" & my_file$species == "G_corn" & my_file$gender == "Male",]
#     x_gcorn <- my_file[my_file$chrom_name != "Un" & my_file$chrom_name != "Unknown" & my_file$species == "G_corn" & my_file$gender == "Male" & my_file$chrom_name == "Neo-X",]
#     x_gcorn$scaled <- x[x$chrom_name == "Neo-X",3]/median(x[x$chrom_name == "Autosome",3])

#     x <- my_file[my_file$chrom_name != "Un" & my_file$chrom_name != "Unknown" & my_file$species == "T_cast" & my_file$gender == "Male",]
#     x_tcast <- my_file[my_file$chrom_name != "Un" & my_file$chrom_name != "Unknown" & my_file$species == "T_cast" & my_file$gender == "Male" & my_file$chrom_name == "Neo-X",]
#     x_tcast$scaled <- x[x$chrom_name == "Neo-X",3]/median(x[x$chrom_name == "Autosome",3])

#     x <- my_file[my_file$chrom_name != "Un" & my_file$chrom_name != "Unknown" & my_file$species == "T_conf" & my_file$gender == "Male",]
#     x_tconf <- my_file[my_file$chrom_name != "Un" & my_file$chrom_name != "Unknown" & my_file$species == "T_conf" & my_file$gender == "Male" & my_file$chrom_name == "Neo-X",]
#     x_tconf$scaled <- x[x$chrom_name == "Neo-X",3]/median(x[x$chrom_name == "Autosome",3])

#     x <- my_file[my_file$chrom_name != "Un" & my_file$chrom_name != "Unknown" & my_file$species == "T_brev" & my_file$gender == "Male",]
#     x_tbrev <- my_file[my_file$chrom_name != "Un" & my_file$chrom_name != "Unknown" & my_file$species == "T_brev" & my_file$gender == "Male" & my_file$chrom_name == "Neo-X",]
#     x_tbrev$scaled <- x[x$chrom_name == "Neo-X",3]/median(x[x$chrom_name == "Autosome",3])

#     x <- my_file[my_file$chrom_name != "Un" & my_file$chrom_name != "Unknown" & my_file$species == "T_frem" & my_file$gender == "Male",]
#     x_tfrem <- my_file[my_file$chrom_name != "Un" & my_file$chrom_name != "Unknown" & my_file$species == "T_frem" & my_file$gender == "Male" & my_file$chrom_name == "Neo-X",]
#     x_tfrem$scaled <- x[x$chrom_name == "Neo-X",3]/median(x[x$chrom_name == "Autosome",3])

#     scaled_Neo_X <- rbind(x_gcorn,x_tbrev,x_tconf,x_tcast,x_tfrem)
