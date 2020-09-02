library(tidyverse)
library(dplyr)
library(limma)
library(edgeR)
library(org.Hs.eg.db)
library(ggrepel)
library(readr)
library(plot3D)
library(goseq)
library(fgsea)

#read in counts table.  There is one column that's the gene, and then one column for each Sample
pA <-
  read_csv("pA_Raw.csv")

#read in the Sample sheet.
sample_names <- colnames(pA)[2:7]
condition <-
  c("infected", "infected", "infected", "mock", "mock", "mock")

Sample_info <- tibble(Sample = sample_names, Condition = condition)
#everything is good to go when all the colnames for the counts table match the sample name column in the Sample_info table
rm(sample_names)
rm(condition)

names(pA)[2:7] == Sample_info$Sample

#This will create an object with the normalized (for read depth) expression levels for each gene in each sample
##
myCPM <- cpm(pA[2:7])
row.names(myCPM) <- pA$X1
pA$X1 <- NULL
rownames(pA) <- row.names(myCPM)

dge <- DGEList(counts = pA)

dim(dge)

#Now, I'm going to start building my variables from the samples.  In this case, we care about condition
condition <- Sample_info$Condition
#group <- paste(sample, condition, sep = "_") <- use group function if you have more than 1 variable

#I use the filterByExprs function in the edgeR package to filter out lowly expressed genes.  Joe did it by hand.  Either works
keep <-
  filterByExpr(dge, condition = condition)  #This takes into account the group that the sample was in

#removes the lowly expressed from the dataset
y <- dge[keep, , keep.lib.sizes = FALSE]
dim(y) #15526 genes kept

# Get log2 counts per million
logcounts <- cpm(y, log = TRUE)

#RAW MDS plot
#plotMDS(y, top = 500) #What is an MDS Plot

#normalization
y <-
  calcNormFactors(dge, method = "TMM") 

design <- model.matrix( ~ 0 + condition) #???

colnames(design) = sub("condition", "", colnames(design)) #???

v <- voom(dge, design, plot = TRUE)
v

fit = lmFit(v, design)

cont.matrix <- makeContrasts(infected - mock,
                             levels = design)

fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)

#
summa.fit <- decideTests(fit.cont)
summary(summa.fit)

results.table.1 <- topTable(fit.cont, coef = 1, number = Inf) %>%
  as_tibble(rownames = "Gene")

results.table.1 <- results.table.1 %>%
  separate(Gene,
           sep = "_",
           into = c("Gene", "Symbol", "extra1", "extra2")) %>%
  dplyr::select(-extra1,-extra2)

results.table.2 <- results.table.1 %>%
  filter(adj.P.Val <= 0.05)

results.table.2$logPVal <- -log10(results.table.2$adj.P.Val)

###############

results.table.2 <- results.table.2 %>%
  filter(is.na(Symbol) == FALSE,
         duplicated(Symbol) == FALSE)

#So, if you want the best list of genes that are significantly differentially expressed:
#We'll also get rid of genes without gene symbols
sig.genes <- results.table.2 %>%
  filter(adj.P.Val < .05)

#If that's a zillion genes, you can add in a logFC cut-off

#THIS IS OPTIONAL
sig.genes <- results.table.2 %>%
  filter(adj.P.Val < .05 &
           abs(logFC) > 2)  #This takes genes whose absolute value of logFC is greater than 1

#to run GO analysis using "changed genes" vs "expressed genes"
#this pulls out all of the expressed genes
all.genes <- results.table.2 %>% pull(Symbol)

#this puts a 1 for genes that are significantly changed and a 0 for those that are not
gene.vector <- as.integer(all.genes %in% sig.genes$Symbol)

#look at it, "1" should be the number of significantly changed genes
table(gene.vector)

#we need to name the values
names(gene.vector) <- all.genes

#look at it
# head(gene.vector)

#run the goseq analysis
pwf = nullp(gene.vector, "hg19", "geneSymbol")
# head(pwf)
GO.wall <- goseq(pwf, "hg19", "geneSymbol")
# head(GO.wall)

#calculate FDR
GO.wall$padj <-
  p.adjust(GO.wall$over_represented_pvalue, method = "BH")

GO.sig <- GO.wall %>%
  filter(over_represented_pvalue < .001) %>%
  mutate(rank = rank(over_represented_pvalue))

ggplot(data = GO.sig %>%
         filter(rank < 50 &
                  #you can adjust this if it doesn't give very many categories
                  numInCat < 1000),
       #this removes really large categories, which can be vague, you can adjust this
       aes(
         x = reorder(term, -over_represented_pvalue),
         y = -log10(over_represented_pvalue)
       )) +
  geom_bar(stat = "identity", fill = "lightgray", color = "black") +
  coord_flip() +
  theme(text = element_text(size = 12), axis.text = element_text(size = 12)) +
  theme_classic()


#another option, this is "GSEA" analysis
#In this version, we rank all of the expressed genes by how much they're changing
#I like the t-statistic from limma results

#For this one, we need to input the pathways
Pathways <- gmtPathways("Pathway_Raw.xls")

#Build the input.  It looks a lot like the GO, except that instead of 1s and 0s, we
#have a value for ranking
gsea.input <- results.table.2$t
names(gsea.input) <- results.table.2$Symbol

#look at it
#head(gsea.input)

#run the analysis.  This may take a while
gsea_res <-
  fgsea(
    Pathways,
    gsea.input,
    minSize = 15,
    maxSize = 500,
    eps = 0
  )

#sort it by NES, which is the enrichment score that indicates how much a given
#pathways is enriched in the top (positive) or bottom (negative) of the list.
#genes at the top are going up, genes at the bottom are going down
gsea_res <- gsea_res[order(gsea_res$NES, decreasing = TRUE), ]

#take the significant ones
gsea_res.sig <- gsea_res %>%
  filter(padj < .05) %>%
  mutate(rank = rank(padj))

#For GSEA, we can plot the individual Pathways, which are cute.
#plotEnrichment(Pathways[["GO_RESPONSE_TO_CYTOKINE"]], gsea.input)
#The ES reflects the degree to which the genes in a gene set are
#over-represented at the top or bottom of the entire ranked list of genes.

gsea_res.top <- gsea_res.sig %>%  top_n(n = 10, wt = NES)
gsea_res.top <-
  rbind(gsea_res.top, top_n(gsea_res.sig, n = 10, wt = -NES))

##########

rm(keep)
rm(gene.vector)
rm(design)
rm(cont.matrix)
rm(condition)
rm(all.genes)
rm(dge)
rm(fit)
rm(fit.cont)
rm(logcounts)
rm(myCPM)
rm(pwf)
rm(Sample_info)
rm(sig.genes)
rm(summa.fit)
rm(y)
rm(v)
rm(GO.wall)
rm(change)
rm(GO.sig)
rm(hg19.geneSymbol.LENGTH)
rm(gsea_res)
rm(results.table.a)

# ggplot(data = gsea_res.top, aes(x=fct_reorder(pathway, abs(NES)), y=abs(NES))) +
#   geom_point(data = gsea_res.top, aes(x=abs(NES),  y=fct_reorder(pathway, -abs(NES)), size = size, color = padj)) +
#   #scale_color_manual(values=c("red", "turquoise3")) +
#   labs(x="Norm. Enrichment Score", y="", colour="Adj. P Val.", size="Count") + theme_minimal() +
#   theme(legend.position = "bottom" )+ theme(axis.text = element_text(size = 5,face='bold'))

# results.table.a <- filter(results.table.2,Symbol %in% Pathways$GO_CONTRACTILE_FIBER)
# change <- results.table.a %>% pull(logFC)  

#this puts a 1 for genes that are significantly changed and a 0 for those that are not
# gene.vector <- change > 0

# with(
#   results.table.a,
#   scatter3D(
#     abs(logFC),
#     logPVal,
#     AveExpr,
#     labels = results.table.a$Symbol,
#     colvar = gene.vector,
#     col = gg.col(2),
#     theta = 50,
#     phi = 10,
#     xlab = "Fold Change",
#     ylab = "Lvl of Significance",
#     zlab = "Avg. Expr. of Gene",
#     main = c("iPSC-Derived CM Increase in Expression (Infected v. Mock)","","","",""),
#     cex = 0.7,
#     bty = "g",
#     ticktype = "detailed",
#     d = 7,
#     clab = c("", "", "Inc.(Black)", "or Dec.(Blue)", "in Expr."),
#     adj = 0.1
#   )
# )
# 
# 
#   scatter2D(
#     results.table.a$logFC,
#     results.table.a$logPVal,
#     labels = results.table.a$Symbol,
#     colvar = results.table.a$AveExpr,
#     col = gg.col(100),
#     theta = 50,
#     phi = 10,
#     xlab = "Fold Change",
#     ylab = "Lvl of Significance",
#     main = c("iPSC-Derived CM Increase in Expression (Infected v. Mock)","","",""),
#     cex = 0.7,
#     bty = "g",
#     ticktype = "detailed",
#     d = 7,
#     clab = c("", "", "Avg. Expr", "of Gene"),
#     adj = 0.1
#   )
ggplot(data = results.table.1, aes(x=logFC, y=adj.P.Val)) +
  geom_point(data = results.table.1, aes(x=logFC, y=-log10(adj.P.Val), color=AveExpr)) +
  scale_color_gradient(low="blue", high="red") +
  ggtitle("iPSC Derived CM Genes")+
  labs(x="log2 Fold Change", y="Stat. Significance", color="Avg. Expr. of Gene") +
  theme_minimal() + theme(axis.text = element_text(size = 10))


ggplot(data=results.table.1, aes(x=logFC, y = -log10(adj.p.Val)))+
  geom_point(data=results.table.1, aes(x=logFC, y = -log10(adj.p.Val), color=AveExpr)) +
  scale_color_gradient(low="blue", high="red") +
  ggtitle("iPSC Derived CM Genes") +
  labs(x="log2 Fold Change", y="Stat. Significance", color="Avg. Expr. of Gene") +
  theme_minimal() + theme(axis.text = element_text(size = 10))

 genes_a <- c(gsea_res.top$leadingEdge[2])
# 
results.table.1 %>%
  #filter(Symbol %in% genes_a) %>%
  ggplot(aes(x=logFC, y=-log10(adj.P.Val)))+
  geom_point(color="grey") +
  geom_point(data=subset(results.table.1, Symbol %in% genes_a[[1]] & adj.P.Val < 0.05), color="blue") +
  #scale_color_gradient(low="blue", high="red") +
  ggtitle("iPSC Derived CM Genes") +
  labs(x="log2 Fold Change", y="Stat. Significance", color="Avg. Expr. of Gene") +
  theme_minimal() + theme(axis.text = element_text(size = 10))
# 
# 
 c<-results.table.1$Symbol %in% genes_a


