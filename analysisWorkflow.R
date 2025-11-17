#required packages, if not installed, single name can be replaced with a vector of lib names e.g. c("lib1","lib2")
#BiocManager::install("clusterProfiler")
#BiocManager::install("pathview")
#BiocManager::install("enrichplot")
#BiocManager::install("func2vis")
#BiocManager::install(organism, character.only = TRUE)
#BiocManager::install('EnhancedVolcano')

#load libraries
library(SPIA)
library(stringr)
library(readxl)
library(pathview)
library(clusterProfiler)
library(msigdbr)
library(EnhancedVolcano)

#organism
organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)

# Keep the starting directory, just in case.
original_direcotry <- getwd()
# Variable name for the analysis. This was meant to analyze TCGA data, however, it can used to analyse other data sets as well.
cName <- "KIRP" # the name of disease/type, or the data to be analyzed

# First load the workspace where the analysis has been stored
# This is meant to be independent. No other variable with same name should be in the environment.
# It is lot better to clean the evn beforehand with rm(list = ls())
load(paste0("TCGA-",cName,".RData", sep =""))

#######################################
# Create directory for Analysis
#######################################

# Create and set the parent dir for analysis; it will create the analysis dir in the current working dir
# A time stampe should be included in the name, so that the re-analysis (if needed) should not overwrite this. Plus, 
# it can be useful to track the analysis date.

parentDir <- paste0(getwd(), "/", cName, "_Analysis", sep = "")
dir.create(parentDir)
setwd(parentDir)

# function for subsequent dir creation and setting working dir
# This will create other directories within "parentDir", thus, only "dirName"
# should be used

create_and_setwd <- function(dirName) {
	if(file.exists(dirName))
		{setwd(file.path(parentDir, dirName)) }
	else { 
		dir.create(file.path(parentDir, dirName))
		setwd(file.path(parentDir, dirName))  
		}
					}

######################################
# Get gene list for the analysis
######################################

# it is often seen/recommended that the reseachers used all the identified proteins as background, instead of global background
background <- row.names(dataFilt) 
dataDEG2 <- dataDEG2[which(dataDEG2$gene_type == "protein_coding"), ]

# list of differentially expressed genes for EnrichGO analysis
genes_only_DEG2 <- as.vector(dataDEG2$gene_name)

# up regulated in cancers
up_in_cancer <- dataDEG2[which(dataDEG2$logFC > 0  ),] # dataDEG2 refers log2FC > abs(2)
up_in_cancer <- as.vector(up_in_cancer$gene_name)

# up regulated in normal
up_in_normal <- dataDEG2[which(dataDEG2$logFC < 0  ),]
up_in_normal <- as.vector(up_in_normal$gene_name)

geneList <- dataDEG2[, c("logFC", "FDR")]

############################################
# Convert geneIDs from ENSEM to ENTREIDs
############################################

# While converting the names, BiTr usually returns a df with two columns, "fromType" & "toType", and it is oftenly difficult to match them
# to the original dataset. This is a simple function around that concept that BiTr should bind to converted IDs back to the original df 
# for easier downstream analysis

bitr2 <- function(df) {
	stopifnot(class(df)  == "data.frame")
	message("Input is not dataframe")
	items <- row.names(df)
	itemsID <- clusterProfiler::bitr(items, fromType="ENSEMBL", 
					toType="ENTREZID", OrgDb=organism, drop=TRUE)
	df <- merge(df, itemsID, by.x = 0, by.y="ENSEMBL", all.x = TRUE)
	df <- transform(df, logFC = as.numeric(logFC), ENTREZID = as.numeric(ENTREZID))
	df <- df[complete.cases(df),]
	df = df[!duplicated(df$ENTREZID),]
	message("after removing NA: ", dim(df))
	df
			}

geneList <- bitr2(geneList)

#background gene list
backgroundID <- clusterProfiler::bitr(background, fromType="ENSEMBL", 
						toType="ENTREZID", OrgDb=organism, drop=TRUE)
background <- backgroundID$ENTREZID
background <- as.numeric(background)

#
# Sorted list of Differentially Expressed Genes; required for GSEA
#

geneList_DegEntrz <- geneList$logFC
names(geneList_DegEntrz) <- as.vector(geneList$ENTREZID)
geneList_DegEntrz <- sort(geneList_DegEntrz, decreasing = TRUE)

##################################
# Volcano Plot
##################################

create_and_setwd("01_VolcanoPlot")

DEGsP <- dataDEGs[which(dataDEGs$gene_type == "protein_coding"), ]

DEGsP$logFC <- as.numeric(as.character(DEGsP$logFC))
DEGsP$FDR <- as.numeric(as.character(DEGsP$FDR))

EnhancedVolcano(DEGsP, lab = DEGsP$gene_name,
			x = 'logFC', FCcutoff = 2,
			y = 'FDR', pCutoff = 10e-6,
			pointSize = 2.0, labSize = 6.0,
			xlab = bquote(~Log[2]~ 'fold change'),
			colAlpha = 1,
    			legendPosition = 'top',
    			legendLabSize = 12,
    			legendIconSize = 4.0,
    			drawConnectors = FALSE,
			widthConnectors = 0.75)

# ggsave will save the last display graph
ggsave("1.volcanoPlot.pdf", width = 10, height = 10, units = "in", dpi = 300)
dev.off()  

#################################
# 1. Gene Ontology Analysis (OverRepresentation analysis)
#################################

# It is not encouraged to do this, unless you are doing it for clustered proteins; the drawbacks in this method can be found here:
# https://pnnl-comp-mass-spec.github.io/proteomics-data-analysis-tutorial/ora.html#ora-drawbacks

create_and_setwd("02_GO_Analysis")

require(clusterProfiler)

geneList_DEG <- geneList$ENTREZID
cpGO <- enrichGO(geneList_DEG,
	OrgDb = organism,
	keyType = "ENTREZID",
	ont = "BP",		# it can be "CC", "MF" and "ALL"
	pvalueCutoff = 0.05,
	pAdjustMethod = "BH",
	universe = as.character(background),
	qvalueCutoff = 0.2,
	minGSSize = 10,
	maxGSSize = 500,
	readable = TRUE,	# will convert EntrezIDs to GeneNames
	pool = FALSE)

write.table(cpGO, file = "1.GOBP_clusterPro.txt", sep = "\t")
cpGO_df <- as.data.frame(cpGO)

#plot a barplot for top 10 enriched terms ordered by q-values
ggplot(cpGO_df[1:10, ], aes(x = -log10(qvalue[1:10]), y =  reorder(Description[1:10], -log10(qvalue[1:10]) )) + 
	geom_bar(stat = "identity") + 
	theme_classic()

ggsave("1.enrichGO_plot.pdf", width = 10, height = 10, units = "in", dpi = 300)
dev.off()

#################################
# 2. SPIA (Ontology based Pathway analysis)
#################################

create_and_setwd("03_SPIA")

#To create SPIAdata for SPIA analysis; since mouse is not provided, it has to be created
#makeSPIAdata(kgml.path="../mouse_KEGG_pathways/", organism="mmu", out.path='.')
#load(file = paste(system.file("extdata/hsaSPIA.RData",package ="SPIA")))

#pathway analysis based on combined evidence; # use nB=2000 or more for more accurate result

#pGFdr < 0.05 shuld be the cutoff as author suggested

res <- spia(de=geneList_DegEntrz,
		all=as.character(background), 
		organism="hsa",
		nB=3000,
		plots=F,
		beta=NULL)

#if it is desired to fit the result for better visualization; description terms will be truncated
#res$Name <- substr(res$Name,1,10)
#show all pathways with pGDFR < 0.05, omit KEGG links

res <- res[res$pGFdr < 0.05, -12]

write.table(res, file="03.SPIA_pGFDR_0.05.txt", sep='\t', row.names=F)

# It will loop over and create PathView annotated plots
for (i in res$ID){ 
		pathview(gene.data=geneList_DegEntrz,
			pathway.id=i,
			species="hsa")
		}

dev.off()
################################
# GSEA using ClusterProfiler	 
################################

create_and_setwd("04_GSEA")


#change the type of columns to numeric, otherwise unary (aristhmatic comparison) doesnot work
geneList <- transform(geneList, logFC = as.numeric(logFC), ENTREZID = as.numeric(ENTREZID))

gse <- gseGO(geneList = geneList_DegEntrz,
		ont = "BP",
		keyType = "ENTREZID",
		minGSSize = 5,
		maxGSSize = 800, 
            	pvalueCutoff = 0.05,
		verbose = FALSE,
		OrgDb = organism,
		pAdjustMethod = "BH",
		eps = 0)

write.table(as.data.frame(gse), file="4.gseGO_padjust_0.05.txt", sep='\t', row.names=F)

# There are few plot that can be created
  
require(DOSE)
dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)
ggsave("4.gseGO_dotplot.pdf", width = 10, height = 10, units = "in", dpi = 300)

#enrichment map
emapplot(gse, showCategory = 10)
# categorySize can be either 'pvalue' or 'geneNum'
cnetplot(gse, categorySize="pvalue", foldChange=geneList_DegEntrz, showCategory = 3)
ggsave("4.gseGO_cnetplot.pdf", width = 10, height = 10, units = "in", dpi = 300)

# Grouped by gene set, density plots are generated by using the frequency of fold change values 
# per gene within each set. Helpful to interpret up/down-regulated pathways.
ridgeplot(gse) + labs(x = "enrichment distribution")
ggsave("4.gseGO_ridgeplot.pdf", width = 10, height = 10, units = "in", dpi = 300)

# Use the `Gene Set` param for the index in the title, and as the value for geneSetId

# loop over to create and save plots
for (i in 1:30) {
	gseaplot(gse, by = "all", title = gse$Description[i], geneSetID = i)
	ggsave(paste(gse$Description[i],'_',gse$Description[i],'_',i,"_gsePlot.pdf" ,sep=''), 
			 width = 10, height = 10, units = "in", dpi = 300)
 		}

dev.off()

###################################################
# gseKEGG; gene set enrichment using KEGG data set
###################################################

create_and_setwd("05_GSEA_with_KEGG")

gsekegg <- gseKEGG(
		geneList = geneList_DegEntrz,
		organism = "hsa",
		keyType = "ncbi-geneid",
		exponent = 1,	
		minGSSize = 10,
		maxGSSize = 500,
		eps = 0,
		pvalueCutoff = 0.05,
		pAdjustMethod = "BH",
		verbose = TRUE,
		use_internal_data = FALSE,
		seed = FALSE,
		by = "fgsea")

write.table(as.data.frame(gsekegg), file="5.gseKEGG_padjust_0.05.txt", sep='\t', row.names=F)

require(DOSE)
dotplot(gsekegg, showCategory=11, split=".sign") + facet_grid(.~.sign)

ggsave("5.gsekegg_dotplot.pdf", width = 10, height = 10, units = "in", dpi = 300)

#Grouped by gene set, density plots are generated by using the frequency of fold change values #per gene within each set. Helpful to interpret up/down-regulated pathways.
ridgeplot(gsekegg) + labs(x = "enrichment distribution")
ggsave("5.gseKEGG_ridgeplot.pdf", width = 10, height = 10, units = "in", dpi = 300)


for (i in 1:length(gsekegg$ID)) {
	gseaplot(gsekegg, by = "all", title = gsekegg$Description[i], geneSetID = i)
	ggsave(paste("5.",gsekegg$Description[i],'_',gsekegg$ID[i],'_',i,"_gseKEGGPlot.pdf" ,sep=''), 
			 width = 10, height = 10, units = "in", dpi = 300)}
dev.off()
############################################
# Gene set enrichment with predefined genes
############################################

create_and_setwd("06_GSEA_ImmuneSigDB")

# The follwing cab be used to see relevant information
# msigdbr_show_species()
# print(msigdbr_collections(), n=30)
# head(m_t2g)

#----------------- 1 ------------------------#

m_t2g <- msigdbr(species = "Homo sapiens", category = "C7", subcategory = "IMMUNESIGDB") %>%
         dplyr::select(gs_name, entrez_gene)

immuneSig_gse <- GSEA(geneList_DegEntrz,
			TERM2GENE = m_t2g,
			pvalueCutoff = 0.05,
			eps=0,
			verbose = FALSE)

write.table(as.data.frame(immuneSig_gse), file="6.immunesigdb_gse.txt", row.names=F, quote=F)

for (i in 1:length(immuneSig_gse$ID)) {
	gseaplot(immuneSig_gse, by = "all", title = immuneSig_gse$Description[i], geneSetID = i)
	ggsave(paste("6.",immuneSig_gse$Description[i],'_',immuneSig_gse$ID[i],'_',i,".pdf" ,sep=''), 
			 width = 10, height = 10, units = "in", dpi = 300)}

dev.off()
#----------------- 2 ------------------------#

create_and_setwd("06_GSEA_HallMarks")

h_t2g <- msigdbr(species = "Homo sapiens", category = "H", subcategory = "") %>%
         dplyr::select(gs_name, entrez_gene)

hallmark_gse <- GSEA(geneList_DegEntrz, TERM2GENE = h_t2g,
				pvalueCutoff = 0.05, eps=0, verbose = FALSE)

write.table(as.data.frame(hallmark_gse), file="hallmark_gse.txt", row.names=F, quote=F)

for (i in 1:length(hallmark_gse$ID)) {
	gseaplot(hallmark_gse, by = "all", title = hallmark_gse$Description[i], geneSetID = i)
	ggsave(paste(hallmark_gse$Description[i],'_',hallmark_gse$ID[i],'_',i,".pdf" ,sep=''), 
			 width = 10, height = 10, units = "in", dpi = 300)}
dev.off()
#----------------- 3 ------------------------#

create_and_setwd("06_GSEA_C3_TFT-gtrd")

c3GTRD <- msigdbr(species = "Homo sapiens", category = "C3", subcategory = "TFT:GTRD") %>%
         dplyr::select(gs_name, entrez_gene)

c3GTRD_gse <- GSEA(geneList_DegEntrz, TERM2GENE = c3GTRD,
				pvalueCutoff = 0.05, eps=0, verbose = FALSE)

write.table(as.data.frame(c3GTRD_gse), file="C3_TFT-GTRD_gse.txt", row.names=F, quote=F)

for (i in 1:length(c3GTRD_gse$ID)) {
	gseaplot(c3GTRD_gse, by = "all", title = c3GTRD_gse$Description[i], geneSetID = i)
	ggsave(paste(c3GTRD_gse$Description[i],'_',c3GTRD_gse$ID[i],'_',i,".pdf" ,sep=''), 
			 width = 10, height = 10, units = "in", dpi = 300)}
dev.off()
#----------------- 4 ------------------------#

create_and_setwd("06_GSEA_HPO")

c5HPO <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "HPO") %>%
         dplyr::select(gs_name, entrez_gene)

c5HPO_gse <- GSEA(geneList_DegEntrz, TERM2GENE = c5HPO,
				pvalueCutoff = 0.05, eps=0, verbose = FALSE)

write.table(as.data.frame(c5HPO_gse), file="C5HPO_gse.txt", row.names=F, quote=F)

for (i in 1:length(c5HPO_gse$ID)) {
	gseaplot(c5HPO_gse, by = "all", title = c5HPO_gse$Description[i], geneSetID = i)
	ggsave(paste(c5HPO_gse$Description[i],'_', c5HPO_gse$ID[i],'_',i,".pdf" ,sep=''), 
			 width = 10, height = 10, units = "in", dpi = 300)}
dev.off()
#----------------- 5 ------------------------#
# C6 -> Oncogenic signature gene set

create_and_setwd("06_GSEA_Oncogenic_C6")

c6 <- msigdbr(species = "Homo sapiens", category = "C6", subcategory = "") %>%
         dplyr::select(gs_name, entrez_gene)

c6_gse <- GSEA(geneList_DegEntrz, TERM2GENE = c6,
				pvalueCutoff = 0.05, eps=0, verbose = FALSE)

write.table(as.data.frame(c6_gse), file="c6_oncogenicSig_gse.txt", row.names=F, quote=F)

for (i in 1:length(c6_gse$ID)) {
	gseaplot(c6_gse, by = "all", title = c6_gse$Description[i], geneSetID = i)
	ggsave(paste("8.",c6_gse$Description[i],'_', c6_gse$ID[i],'_',i,".pdf" ,sep=''), 
			 width = 10, height = 10, units = "in", dpi = 300)}

dev.off()
#############################################################
# EnrichR
#############################################################

create_and_setwd("07_EnrichR_GO")

#install.packages("enrichR")
library(enrichR)

#by default, species == human 

#Then find the list of all available databases from Enrichr.
#dbs <- listEnrichrDbs()
#listEnrichrSites()

enrichr_GO <- function(geneNames, fName) {

	websiteLive <- TRUE

	dbs <- c("GO_Molecular_Function_2021",	# 1
		   "GO_Cellular_Component_2021",	# 2
		   "GO_Biological_Process_2021",	# 3 
		   "KEGG_2021_Human", 			# 4
		   "Metabolomics_Workbench_Metabolites_2022", 	#5
		   "ChEA_2022", 						# 6	
		   "KEA_2015"						# 7
			)

if (websiteLive) { enriched <- enrichr(geneNames, dbs) }

for (i in 1:length(dbs)){
	write.table(enriched[[i]][which(enriched[[i]]$Adjusted.P.value < 0.05),], paste0(fName,"_",dbs[i], ".txt", sep=""), sep = "\t")
				}
# For plot the GO enrichments
#if (websiteLive) plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
#enriched

}

#GOMF_2021   <- write.table(enriched[[1]][which(enriched[[1]]$Adjusted.P.value < 0.05),], "GOMF_2021.txt", sep = "\t")

enrichr_GO(up_in_cancer, "cancerUp")
enrichr_GO(up_in_normal, "normalUp")
enrichr_GO(genes_only_DEG2, "DEGonly")


############################################################
# Remove all the variables to keep the envireonment clean
############################################################

remove_all <- function() {
	message("
			At the end of this analysis, it is desireable to delete all the variables.
			This way, the next session will not be effected by any varible created 
			in this session
						")

	response = readline(prompt = "Yes | No: ")
	if(tolower(response) == tolower("Yes")) {rm(list = ls())}
	
	else { message("Please carefull, as no variables has been removed")}
	}


remove_all()
Yes

rm(list = ls())
setwd(original_direcotry)

