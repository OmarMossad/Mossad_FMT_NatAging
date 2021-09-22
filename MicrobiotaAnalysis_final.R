#set working directory to where the files are
# Load libraries
library(ggplot2)
library(phyloseq)
library(microbiome)
library(microbiomeutilities)
library(knitr)
library(taxa)
library(metacoder)
library(gridExtra)
library(vegan)
library(phangorn)
library(genefilter)
library(dplyr)
library(plyr)
library(metagMisc)
library(ggpubr)
library(DESeq2)
library(dada2)
library(msa)
library(corrplot)
library(ggcorrplot)
library(RColorBrewer)
library(agricolae)
library(Hmisc)
#####
#Fecal pellet data

#Import
Meta_all <- read.csv(file="Metadata_All.csv", header = TRUE, row.names=1) 
meta_fp <- Meta_all %>% filter(Tissue=="FP")
write.csv(meta_fp, file="Meta_FP.csv")

df <- read.csv(file="DataframeForDiff_all.csv", header = TRUE, row.names=1)
df1 <- df[colnames(df) %in% rownames(meta_fp)]
write.csv(df1, file="DataframeForDiff_FP.csv")

dim(df1)
df1_filtered <- df1[rowSums(df1) > 0, ]  # filtered row read count above 0
dim(df1_filtered)
write.csv(df1_filtered, file="DataframeForDiff_FP_filtered.csv")

pseq <- read_phyloseq(otu.file = "DataframeForDiff_FP_filtered.csv", taxonomy.file = "BacteriaAnnotation_FP.csv", metadata.file = "Meta_FP.csv", type = "simple")

# Show available ranks in the dataset
rank_names(pseq)

# Create table, number of features for each phyla
table(tax_table(pseq)[, "Phylum"], exclude = NULL)

#filtering and normalization
#pseq_fil <- prune_taxa(taxa_sums(pseq) > 4, pseq)
wh0 = genefilter_sample(pseq, filterfun_sample(function(x) x > 0), A=0.1*nsamples(pseq))
pseq_fil = prune_taxa(wh0, pseq)
pseq_fil <- subset_taxa(pseq_fil, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
pseq_norm = transform_sample_counts(pseq_fil, function(x) 1E4 * x/sum(x))

# The sample counts stored within the phyloseq pseq_normect
colSums(otu_table(pseq_norm))

# phyloseq holds all information within one R pseq_normect
str(pseq_norm)

#export df from pseq_norm
normalizedDF <- phyloseq_to_df(pseq_norm)
write.csv(normalizedDF, file="DataframeForDiff_FP_normalized.csv")

# How many genera would be present after filtering?
length(get_taxa_unique(pseq, taxonomic.rank = "Genus"))
length(get_taxa_unique(pseq_fil, taxonomic.rank = "Genus"))
length(get_taxa_unique(pseq_norm, taxonomic.rank = "Genus"))
ntaxa(pseq)
ntaxa(pseq_fil)
ntaxa(pseq_norm)

#beta diversity
ordinate(pseq_norm, "NMDS", "bray") %>% 
  plot_ordination(pseq_norm, ., color='Group', shape = "Timepoint", title = "Bray-Curtis") + geom_point(size=5) + stat_ellipse()

#####
#Cecum data
#Import
Meta_all <- read.csv(file="Metadata_All.csv", header = TRUE, row.names=1) 
meta_cc <- Meta_all %>% filter(Tissue=="cc")
write.csv(meta_fp, file="Meta_cecum.csv")

df <- read.csv(file="DataframeForDiff_all.csv", header = TRUE, row.names=1)
df_cc <- df[colnames(df) %in% rownames(meta_cc)]
write.csv(df_cc, file="DataframeForDiff_cecum.csv")
dim(df_cc)
df_filtered <- df_cc[rowSums(df_cc) > 0, ]  # filtered row read count above 0
dim(df_filtered)
write.csv(df_filtered, file="DataframeForDiff_cecum_filtered.csv")

pseq <- read_phyloseq(otu.file = "DataframeForDiff_cecum.csv", taxonomy.file = "BacteriaAnnotation.csv", metadata.file = "Metadata_cecum.csv", type = "simple")

# Show available ranks in the dataset
rank_names(pseq)

# Create table, number of features for each phyla
table(tax_table(pseq)[, "Phylum"], exclude = NULL)

#filtering and normalization
#pseq_fil <- prune_taxa(taxa_sums(pseq) > 4, pseq)
wh0 = genefilter_sample(pseq, filterfun_sample(function(x) x > 0), A=0.1*nsamples(pseq))
pseq_fil = prune_taxa(wh0, pseq)
pseq_fil <- subset_taxa(pseq_fil, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
pseq_norm = transform_sample_counts(pseq_fil, function(x) 1E4 * x/sum(x))

# The sample counts stored within the phyloseq pseq_normect
colSums(otu_table(pseq_norm))

# phyloseq holds all information within one R pseq_normect
str(pseq_norm)

#export df from pseq_norm
normalizedDF <- phyloseq_to_df(pseq_norm)
write.csv(normalizedDF, file="DataframeForDiff_cecum_normalized.csv")

# How many genera would be present after filtering?
length(get_taxa_unique(pseq, taxonomic.rank = "Genus"))
length(get_taxa_unique(pseq_fil, taxonomic.rank = "Genus"))
length(get_taxa_unique(pseq_norm, taxonomic.rank = "Genus"))
ntaxa(pseq)
ntaxa(pseq_fil)
ntaxa(pseq_norm)

# Plotting a stacked bar chart of taxon abundance
p <- plot_bar(pseq_norm, fill="Phylum")
pd <- p$data
ggplot(pd, aes(x = Sample, y = Abundance, fill = Phylum)) + scale_fill_brewer(palette = "Set3")+
  geom_bar(stat = "identity") +
  scale_y_continuous(
    limits =   c(0,10000),
    expand = expansion(mult = c(0,0.01)))

# Plotting a stacked bar chart of taxon abundance
plot_bar(pseq_norm, fill="Genus")
#making colors look better
# Define the number of colors you want
nb.cols <- length(get_taxa_unique(pseq_norm, taxonomic.rank = "Genus")) 
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
# Create a ggplot with 18 colors 
p <- plot_bar(pseq_norm, fill="Genus")
pd <- p$data
ggplot(pd, aes(x = Sample, y = Abundance, fill = Genus)) + scale_fill_manual(values = mycolors, na.value="gray67")+
  geom_bar(stat = "identity")+ 
  scale_y_continuous(
    limits =   c(0,10000),
    expand = expansion(mult = c(0,0.01))
  )

# Plotting a heatmap of taxon abundance
plot_heatmap(pseq_norm)
plot_heatmap(pseq_norm, taxa.label='Genus')

#Global indicators
tab <-microbiome::alpha(pseq_norm, index = "all")
kable(head(tab))
write.csv(tab, file="GlobalIndicators.csv")

#Alpha diversity
tab <-microbiome::alpha(pseq_norm, index = "all")
kable(head(tab))
write.csv(tab, file="AlphaDiversity.csv")

#Richness
tab <- richness(pseq_norm)
kable(head(tab))
write.csv(tab, file="Richness.csv")

# Absolute abundances for the single most abundant taxa in each sample
tab <- dominance(pseq_norm, index = "all")
kable(head(tab))
write.csv(tab, file="Dominance.csv")

#function to list the dominating (most abundant) taxa in each sample.
dominant(pseq_norm)

#Rarity and low abundance
tab <- rarity(pseq_norm, index = "all")
kable(head(tab))
write.csv(tab, file="rare.csv")


#Testing differences in alpha diversity
# Construct the data
d <- meta(pseq_norm)
d$Group <- factor(d$Group)
d$shannon <- alpha(pseq_norm, "shannon")$diversity_shannon
d$simpson <- alpha(pseq_norm, "gini_simpson")$diversity_gini_simpson

hist(d$shannon)
hist(d$simpson)

# create a list of pairwise comaprisons
Group <- levels(d$Group) # get the variables

# make a pairwise list that we want to compare.
Group.pairs <- combn(seq_along(Group), 2, simplify = FALSE, FUN = function(i)Group[i])

print(Group.pairs)
#Violin
#library(ggpubr)
#p1 <- ggviolin(d, x = "Group", y = "simpson", add = "boxplot", add.params = list(fill = "white"), fill = "Group") 
#print(p1)
#Statistics
#p1 <- p1 + stat_compare_means(comparisons = Group.pairs) 
#print(p1)

#Box plot
p1 <- ggboxplot(d, x = "Group", y = "shannon",
                color = "Group",
                add = "dotplot")+
  ggtitle("Shannon") +
  xlab("Group") +
  ylab("Alpha diversity index")
#Statistics
p1 <- p1 + p1 <- p1 + stat_compare_means(method = "anova", aes(label=..p.adj..), comparisons = Group.pairs)

#Box plot
p2 <- ggboxplot(d, x = "Group", y = "simpson",
                color = "Group",
                add = "dotplot")+
  ggtitle("Simpson") +
  xlab("Group") +
  ylab("Alpha diversity index")
#Statistics
p2 <- p2 + stat_compare_means(method = "anova", aes(label=..p.adj..), comparisons = Group.pairs)

# group plots together
grid.arrange(nrow = 1, p1, p2)

#Alpha diversity from phyloseq
plot_richness(pseq, x = "Group", color='Group', measures=c("Shannon", "Simpson"))

#beta diversity
ordinate(pseq_norm, "NMDS", "bray") %>% 
  plot_ordination(pseq_norm, ., color='Group', title = "Bray-Curtis") + geom_point(size=5) + stat_ellipse()

#Comparing taxon abundance with 2 groups
sample_data(pseq_norm)$Group <- as.factor(sample_data(pseq_norm)$Group)
diagdds = phyloseq_to_deseq2(pseq_norm, ~ Group)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
#choose groups basic comparisons
resultsNames(diagdds)

#compare aged vs young
res = results(diagdds, cooksCutoff = FALSE, name = "Group_b_Old_vs_a_Young")
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(pseq_norm)[rownames(sigtab), ], "matrix"))
head(sigtab)
write.csv(sigtab, file="Diff_OTUs_old_vs_Young.csv")
#plotting DeOTUs
theme_set(theme_bw())
sigtabgen = subset(sigtab, !is.na(Genus))
# Phylum order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels=names(x))
# Genus order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels=names(x))
ggplot(sigtabgen, aes(y=Genus, x=log2FoldChange, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=5) + 
  ggtitle("Aged vs Young") +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

#compare oFMT vs yFMT
res = results(diagdds, cooksCutoff = FALSE, contrast = c("Group","d_oFMT", "c_yFMT"))
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(pseq_norm)[rownames(sigtab), ], "matrix"))
head(sigtab)
write.csv(sigtab, file="Diff_OTUs_oFMT_vs_yFMT.csv")
#plotting DeOTUs
theme_set(theme_bw())
sigtabgen = subset(sigtab, !is.na(Genus))
# Phylum order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels=names(x))
# Genus order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels=names(x))
ggplot(sigtabgen, aes(y=Genus, x=log2FoldChange, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=5) + 
  ggtitle("oFMT vs yFMT") +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

#intersect of both df
df_oldvsyoung <- read.csv("Diff_OTUs_Old_vs_Young.csv", row.names=1)
df_oldvsyoung_up <- subset(df_oldvsyoung, df_oldvsyoung$log2FoldChange>0)
write.csv(df_oldvsyoung_up, file="df_oldvsyoung_up.csv")
df_oldvsyoung_down <- subset(df_oldvsyoung, df_oldvsyoung$log2FoldChange<0)
write.csv(df_oldvsyoung_down, file="df_oldvsyoung_down.csv")

df_oFMTvsyFMT <- read.csv("Diff_OTUs_oFMT_vs_yFMT.csv", row.names=1)
df_oFMTvsyFMT_up <- subset(df_oFMTvsyFMT, df_oFMTvsyFMT$log2FoldChange>0)
write.csv(df_oFMTvsyFMT_up, file="df_oFMTvsyFMT_up.csv")
df_oFMTvsyFMT_down <- subset(df_oFMTvsyFMT, df_oFMTvsyFMT$log2FoldChange<0)
write.csv(df_oFMTvsyFMT_down, file="df_oFMTvsyFMT_down.csv")

#get list of overlapping OTUs
df_oldvsyoung$OTU <- row.names(df_oldvsyoung)
df_oFMTvsyFMT$OTU <- row.names(df_oFMTvsyFMT)
interesting <- intersect(df_oldvsyoung$OTU,df_oFMTvsyFMT$OTU)

df_oldvsyoung_fil <- subset(df_oldvsyoung, rownames(df_oldvsyoung) %in% interesting)
df_oldvsyoung_fil$Comparison <- "Old_vs_Young"

df_oFMTvsyFMT_fil <- subset(df_oFMTvsyFMT, rownames(df_oFMTvsyFMT) %in% interesting)
df_oFMTvsyFMT_fil$Comparison <- "oFMT_vs_yFMT"

new_df <- rbind(df_oldvsyoung_fil, df_oFMTvsyFMT_fil)
#plotting Diff. OTUs
theme_set(theme_bw())
sigtabgen = new_df #subset(new_df, !is.na(Genus))
# Phylum order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels=names(x))
# Genus order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels=names(x))
ggplot(sigtabgen, aes(y=Genus, x=log2FoldChange, color=Phylum, shape = Comparison)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=5) + 
  ggtitle("Diff. Abundant Taxa") +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

#####
#Construct the phylogenetic tree
seqtab <- read.csv(file="SeqTable.csv")
df_filtered$OTU <- rownames(df_filtered)
seqtab <- subset(seqtab, OTU %in% df_filtered$OTU)
rownames(seqtab) <- seqtab[,1]
seqtab <- subset(seqtab, select=-c(1))
seqs <- getSequences(seqtab$sequence)
names(seqs) <- seqs # This propagates to the tip labels of the tree

mult <- msa(seqs, method="ClustalW", type="dna", order="input", verbose = TRUE)

phang.align <- as.phyDat(mult, type="DNA", names=seqtab$sequence)
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)

## negative edges length changed to 0!

fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
detach("package:phangorn", unload=TRUE)


filtered_df <- phyloseq_to_df(pseq_fil)

samdf <- read.csv(file="Metadata_cecum.csv")
samdf <- merge(seqtab, samdf)
rownames(samdf) <- samdf[,1]
samdf <- subset(samdf, select=-c(1))
write.csv(samdf, file="samdf.csv")


taxtab <- read.csv(file="BacteriaAnnotation.csv")
taxtab <- subset(taxtab, OTU %in% normalizedDF$OTU)
rownames(taxtab) <- taxtab[,1]
taxtab <- subset(taxtab, select=-c(1))
taxtab <- merge(seqtab, taxtab, by=0)
taxtab <- subset(taxtab, select=-c(1))
write.csv(taxtab, file="taxtab.csv")


otutab <- subset(filtered_df, select=-c(2,3,4,5,6,7,8))
rownames(otutab) <- otutab[,1]
otutab <- subset(otutab, select=-c(1))
otutab <- merge(seqtab, otutab, by=0)
otutab <- subset(otutab, select=-c(1))
write.csv(otutab, file="otutab_fil.csv")


taxtab <- subset(taxtab, sequence %in% otutab$sequence)

tree <- fitGTR$tree

pseq_tree = phyloseq(otu_table(as.matrix(otutab), taxa_are_rows=TRUE), 
                     sample_data(samdf), 
                     tax_table(as.matrix(taxtab)), phy_tree(tree, errorIfNULL = T))

random_tree = rtree(ntaxa(pseq_tree), rooted=TRUE, tip.label=taxa_names(pseq_tree))
plot(random_tree)

physeq1 = merge_phyloseq(pseq_tree, samdf, tree)
physeq1
dna <- Biostrings::DNAStringSet(taxa_names(physeq1))
names(dna) <- taxa_names(physeq1)
physeq1 <- merge_phyloseq(physeq1, dna)
taxa_names(physeq1) <- paste0("OTU_", seq(ntaxa(physeq1)))
physeq1

taxa_names(physeq1) <- paste0(seq(ntaxa(physeq1)))

# Convert to taxmap
test <- parse_phyloseq(physeq1)

test$filter_taxa(taxon_names == "Bacteria", subtaxa = TRUE)

# Convert counts to proportions
test$data$otu_table <- calc_obs_props(test,
                                      data = "otu_table",
                                      cols = test$data$sample_data$sample_id)

# Calculate per-taxon proportions 
test$data$tax_table <- calc_taxon_abund(test, 
                                        data = "otu_table", 
                                        cols = test$data$sample_data$sample_id)

print(test$data$tax_table)
# Calculate difference between treatments
test$data$diff_table <- compare_groups(test,
                                       data = "tax_table",
                                       cols = test$data$sample_data$sample_id,
                                       groups = test$data$sample_data$Group2,
)
print(test$data$diff_table)


test <- mutate_obs(test, "diff_table",
                   wilcox_p_value = p.adjust(wilcox_p_value, method = "fdr"))

test$data$diff_table$log2_median_ratio[test$data$diff_table$wilcox_p_value > 0.05] <- 0

#general plot
color_interval <- c(-5, 5) # The range of values (log 2 ratio of median proportion) to display
test %>%
  metacoder::filter_taxa(taxon_names == "Bacteria", subtaxa = TRUE) %>%
  heat_tree(node_size_axis_label = "Number of OTUs",
            node_size = n_obs,
            node_color_axis_label = "Log 2 ratio of median proportions",
            node_color = log2_median_ratio,
            node_color_range = diverging_palette(),
            node_color_trans = "linear",
            node_color_interval = color_interval,
            edge_color_interval = color_interval,
            node_label = taxon_names,
            node_label_max = 150)

#Separate plots
test %>%
  metacoder::filter_taxa(taxon_names == "Bacteria", subtaxa = TRUE) %>%
  heat_tree_matrix(data = "diff_table",
                   node_label = taxon_names,
                   node_size = n_obs, # number of OTUs
                   node_color = log2_median_ratio, # difference between groups
                   node_color_trans = "linear",
                   node_color_interval = c(-5, 5), # symmetric interval
                   edge_color_interval = c(-5, 5), # symmetric interval
                   node_color_range = diverging_palette(), # diverging colors
                   node_size_axis_label = "OTU count",
                   node_color_axis_label = "Log 2 ratio of median counts",
                   layout = "da", initial_layout = "re",
                   key_size = 0.67,
                   seed = 2)

#####
#Stats for correlation
microbe <- subset(normalizedDF, OTU %in% interesting)
row.names(microbe) <- microbe$OTU
microbe <- subset(microbe, select= -c(1,2,3,4,5,6,7,8))
head(microbe)
microbe <- t(microbe)

metab <- read.csv(file= "Metab.csv", row.names = 1)
metab <- log(1 + metab, base = 10)
head(metab)
metab <- t(metab)

cormat <- cbind(microbe, metab)
res <- cor(cormat, use = "complete.obs")
round(res, 2)

# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

res2 <- rcorr(as.matrix(cormat))
res2
#plot correlations

# Insignificant correlations are leaved blank
corrplot(res2$r, type="lower", order="hclust", 
         p.mat = res2$P, tl.col="black", insig = "label_sig",
         sig.level = c(.001, .01, .05), pch.cex = 1, pch.col = "white", col=colorRampPalette(c("darkblue","white","darkred"))(200)) #col=brewer.pal(n=8, name="RdYlBu")

ggcorrplot(res2$r, type="upper", hc.order = TRUE, method = "circle", 
           p.mat = res2$P, insig = "blank",
           sig.level = c(.05), pch.cex = 1.5, pch.col = "white") #col=brewer.pal(n=8, name="RdYlBu")

#extract significantly correlated OTUs
res3 <- flattenCorrMatrix(res2$r, res2$P)
res3 <- filter(res3, column == "Valero")
write.csv(res3, file="correlation_microbes_valero.csv")
res4 <- filter(res3, p<0.05)
write.csv(res4, file="correlation_sig_microbes_valero.csv")

#get the data
BacteriaAnnotation <- read.csv(file="BacteriaAnnotation.csv")
names(res4)[names(res4) == "row"] <- "OTU"
df <- merge(res4, BacteriaAnnotation, by="OTU")
write.csv(df, file="OTUsforValero.csv")

#plotting Diff. OTUs
theme_set(theme_bw())
sigtabgen = subset(new_df, OTU %in% res4$OTU)
# Phylum order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels=names(x))
# Genus order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels=names(x))
ggplot(sigtabgen, aes(y=Genus, x=log2FoldChange, color=Phylum, shape = Comparison)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=5) + 
  ggtitle("Diff. correlated Taxa") +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

