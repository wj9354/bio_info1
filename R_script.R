setwd("C:/Users/kwj93/bio_info1")


if (!requireNamespace("biomaRt", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("biomaRt")
}

library(biomaRt)


ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

enst_list <- readLines("PABPC4_binding_genes.txt")
enst_list_noversion <- sub("\\..*", "", enst_list) # 버전 제거


result <- getBM(
  attributes = c("ensembl_transcript_id", "ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_transcript_id",
  values = enst_list_noversion,
  mart = ensembl
)


symbols <- unique(result$hgnc_symbol[result$hgnc_symbol != ""])
write(symbols, file = "binding_genes_symbol.txt")


if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
  BiocManager::install("clusterProfiler")
}

if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  BiocManager::install("org.Hs.eg.db")
}

library(clusterProfiler)
library(org.Hs.eg.db)


ego_BP <- enrichGO(
  gene         = symbols,
  OrgDb        = org.Hs.eg.db,
  keyType      = "SYMBOL",
  ont          = "BP",   # "BP": Biological Process, "MF", "CC"도 가능
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE
)

head(ego_BP)


ego_MF <- enrichGO(
  gene         = symbols,
  OrgDb        = org.Hs.eg.db,
  keyType      = "SYMBOL",
  ont          = "MF",   # "BP": Biological Process, "MF", "CC"도 가능
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE
)

head(ego_MF)


ego_CC <- enrichGO(
  gene         = symbols,
  OrgDb        = org.Hs.eg.db,
  keyType      = "SYMBOL",
  ont          = "CC",   # "BP": Biological Process, "MF", "CC"도 가능
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE
)

head(ego_CC)


entrez_ids <- bitr(symbols, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
kegg <- enrichKEGG(
  gene = entrez_ids$ENTREZID,
  organism = "hsa"
)
barplot(kegg, showCategory = 10, title = "KEGG enrichment")


df_ego_BP = (as.data.frame(ego_BP))
df_ego_MF = (as.data.frame(ego_MF))
df_ego_CC = (as.data.frame(ego_CC))
df_kegg = (as.data.frame(kegg))

if (!requireNamespace("enrichR", quietly = TRUE)) {
  install.packages("enrichR")
}

library(enrichR)

dbs <- enrichR::listEnrichrDbs()
dbs[grep("Panther", dbs$libraryName), ]

enrich_results <- enrichr(symbols, databases = c("Panther_2016"))
df_panther = as.data.frame(enrich_results)

install.packages("gprofiler2")
library(gprofiler2)


gostres <- gost(symbols, organism="hsapiens", correction_method="fdr")
head(gostres$result)

gostplot(gostres, capped = TRUE, interactive = FALSE)

df_gost = (as.data.frame(gostres))



####
library(gprofiler2)
library(clusterProfiler)
library(msigdbr)
library(enrichplot)

gost <- gost(symbols,
             organism = "hsapiens",
             correction_method="fdr")

myGOplot <- gostplot(gost,
                     interactive = T,
                     capped = T)

table(gost$result$source)

go.bp <- gost$result %>% filter(source == "GO:BP")
go.cc <- gost$result %>% filter(source == "GO:CC")
go.mf <- gost$result %>% filter(source == "GO:MF")
go.kegg <- gost$result %>% filter(source == "KEGG")
go.hp <- gost$result %>% filter(source == "HP")
go.hpa <- gost$result %>% filter(source == "HPA")
go.mirna <- gost$result %>% filter(source == "MIRNA")
go.reac <- gost$result %>% filter(source == "REAC")
go.tf <- gost$result %>% filter(source == "TF")
go.wp <- gost$result %>% filter(source == "WP")




table(gost$result$source)
