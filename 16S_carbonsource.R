##### 16S rRNA sequencing results #####
#Analyzed based on https://github.ugent.be/LabMETNGS/OTUgenRator/blob/master/NoMeta.R


#location
setwd("/media/projects2/CristinaGT/MicrobialProtein_Myrsini/16Sseq/fastq/")

#libs
library(data.table)
library(plyr)
library(tidyr)
library(stringr)
library(dplyr)
library(readxl) 
library(ggplot2)
library(vegan)
library(scales)
library(grid)
library(reshape2)
library(phyloseq)
library(cowplot)

##### generate OTU table #####

sharedfile = "stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared"
taxfile = "stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.cons.taxonomy"
repfile = "stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.rep.fasta"

#generate merged table shared and taxonomy
otu.trans <- t(fread(sharedfile, header = T)) 
tax.table <- fread(taxfile, header = T)
tax.table.split <- tax.table %>%
  separate(Taxonomy, c("regnum", "phylum", "classis", "ordo", "familia", "genus"), ";", extra = "drop")
otu.trans.crop <- otu.trans[-1, ] # remove label
otu.trans.crop <- otu.trans.crop[-2, ] # remove OTU count
colnames(otu.trans.crop) = otu.trans.crop[1, ] # the first row will be the header
otu.trans.crop <- otu.trans.crop[-1, ] # remove sample row
otu.df = data.frame(otu.trans.crop)
setDT(otu.df, keep.rownames = TRUE)[]
otu.df <- plyr::rename(otu.df, c("rn"="OTU"))
merged.table <- merge(tax.table.split, otu.df, by="OTU",all=TRUE)

#filter out concencus sequences
rep.table <- read.table(repfile, header = FALSE, sep = "\t", col.names = paste0("V",seq_len(2)), fill = TRUE)
rep.table.crop.col1 <- subset(rep.table, select = c(V1) )
rep.table.crop.col2 <- subset(rep.table, select = c(V2) )
rep.table.crop.seq <- as.data.frame(rep.table.crop.col1[seq(0, nrow(rep.table.crop.col1), 2), ])
rep.table.crop.id <- as.data.frame(rep.table.crop.col1[seq(1, nrow(rep.table.crop.col1), 2), ])
rep.table.crop.ori <- as.data.frame(rep.table.crop.col2[seq(1, nrow(rep.table.crop.col2), 2), ])

#columns have shitty names
colnames(rep.table.crop.seq)[1] <- "Consensus"

#merge concencus into otu table
merged.rep <- cbind(merged.table,rep.table.crop.seq)

#spit out otu table
dir.create("./output")
write.csv2(file="./output/otu_table.csv",merged.rep)


##### create phyloseq objects and simple graphs #####

### load metadata ###

# map <- as.data.frame(read.csv2("20160617biostat sample list.csv"))

# convert map to phyloseq object

# map <- sample_data(map)

# Assign rownames to be Sample ID's

# rownames(map) <- map$ID

# Import mothur data
mothur_data <- import_mothur(mothur_shared_file = sharedfile,
                             mothur_constaxonomy_file = taxfile)
sample_names(mothur_data)

# Merge mothurdata object with sample metadata
moth_merge <- merge_phyloseq(mothur_data)

# Rename the comumns to reflect tax
colnames(tax_table(moth_merge)) <- c("Kingdom", "Phylum", "Class", 
                                     "Order", "Family", "Genus")

# We can do substampling here

# moth_sub <- moth_merge %>%
#  subset_samples(Type == "sample") %>%
#  prune_taxa(taxa_sums(.) > 0, .)

# Filter out unusefull lineages

# phylobj <- moth_sub %>%
#   subset_taxa(
#     Kingdom == "Bacteria" &
#       Family  != "mitochondria" &
#       Class   != "Chloroplast"
#   )

# phylobj

# bypass last 2 steps

phylobj <- moth_merge

# Make a data frame with a column for the read counts of each sample
sample_sum_df <- data.frame(sum = sample_sums(phylobj))


df_moth_merge <- data.frame(phylobj@otu_table)

sample_names(moth_merge)
# Histogram of sample read counts

ggplot(sample_sum_df, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 2500) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank())

# mean, max and min of sample read counts
smin <- min(sample_sums(phylobj))
smean <- mean(sample_sums(phylobj))
smax <- max(sample_sums(phylobj))

### make some plots ###

#Select top OTUs
TopNOTUs <- names(sort(taxa_sums(phylobj), TRUE)[1:15])
physeqobj15  <- prune_taxa(TopNOTUs, phylobj)
relabsphyseq <- transform_sample_counts(phylobj, function(x) x/sum(x))
TopNOTUs.relab <- names(sort(taxa_sums(relabsphyseq), TRUE)[1:15])
physeqobj15.relab  <- prune_taxa(TopNOTUs.relab, relabsphyseq)


#barplots of top 25 most abundant OTUs

# plot_bar(physeqobj12, x="days", fill = "Family")

nsamp <- ((as.numeric(ncol(df_moth_merge)))*30)+300
dev.off
AbsAbunGen = plot_bar(physeqobj15, fill="Genus", title="Absolute abundance top 15 OTU's")
AbsAbunGen + guides(fill=guide_legend(ncol=1)) + 
  theme_classic()+
  theme(text = element_text(size=18))+
  scale_x_discrete(labels=c("AA lag","AA log","AA stat", "FA lag","FA log","FA stat"))

dev.print(png, width=nsamp, height=600, './output/AbsAbunGen.png') 


dev.off

RelAbunGen = plot_bar(physeqobj15.relab, fill="Genus", title="Relative abundance top 15 OTU's")
RelAbunGen + guides(fill=guide_legend(ncol=1)) + theme(axis.text.x=element_text(angle=90,hjust=0.95,vjust=0.3))
dev.print(png, width=nsamp, height=600, './output/RelAbunGen.png')
dev.off

RelAbunFam = plot_bar(physeqobj15.relab, fill = "Family", title="Relative abundance top 15 OTU's")
RelAbunFam + guides(fill=guide_legend(ncol=1)) + 
  theme_classic()+
  theme(text = element_text(size=18))+
  scale_x_discrete(labels=c("AA lag","AA log","AA stat", "FA lag","FA log","FA stat"))


dev.print(png, width=nsamp, height=600, './output/RelAbunFam.png')
dev.off

