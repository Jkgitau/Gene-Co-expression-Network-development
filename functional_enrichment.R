#########################################################################################
library(biomaRt)
source(here::here("scripts","analysis","libraries.R"))
source("/home/kwanjau/tbrucei_gcn/scripts/analysis/libraries.R")
source("/home/kwanjau/tbrucei_gcn/scripts/utils/enrichment_analysis.R")
source("/home/kwanjau/tbrucei_gcn/scripts/utils/annotations.R")
source("/home/kwanjau/tbrucei_gcn/scripts/utils/wgcna.R")
source("/home/kwanjau/tbrucei_gcn/scripts/utils/util.R")





countsNA <- read.csv('countsNA.txt', sep = "\t")
# try to see what data is in bioMart for Bos taurus
# http://bioconductor.riken.jp/packages/3.4/bioc/vignettes/biomaRt/inst/doc/biomaRt.html
ensembl=useMart("ensembl")
ensembl = useDataset("btaurus_gene_ensembl",mart=ensembl)

# see a list of possible values for filters and attributes
filters = listFilters(ensembl)
ensembl_attributes <- listAttributes(ensembl)

# Attributes and filters to be used in the bioMart query
btaurus_attributes <- c("ensembl_gene_id","hgnc_symbol","description","external_gene_name")
btaurus_filters <- "ensembl_gene_id"

# get the gene ids to be obtained from the query
# btaurus_values <- annotation.ensembl.symbol$ensembl_gene_id
btaurus_values <- countsNA$ENSEMBL_GeneID

# Query
btaurus_ensembl_anno <- getBM(attributes = btaurus_attributes, 
                              filters = btaurus_filters, 
                              values = btaurus_values, 
                              mart = ensembl)


# Query for GO annotation---------------------

btaurus_attributes_GO <- c("ensembl_gene_id","go_id","namespace_1003","name_1006")
btaurus_filters_GO <- "ensembl_gene_id"

# get the gene ids to be obtained from the query
# btaurus_values_GO <- annotation.ensembl.symbol$ensembl_gene_id
btaurus_values_GO <- countsNA$ENSEMBL_GeneID


# Query for GO annotations
btaurus_ensembl_anno_GO <- getBM(attributes = btaurus_attributes_GO, 
                                 filters = btaurus_filters_GO, 
                                 values = btaurus_values_GO, 
                                 mart = ensembl)

#-----------
# Query for transcript lengths

btaurus_attributes_transcript_len <- c("ensembl_gene_id","transcript_length")
btaurus_filters_transcript_len <- "ensembl_gene_id"

# get the gene ids to be obtained from the query
# btaurus_values_transcript_len <- annotation.ensembl.symbol$ensembl_gene_id
btaurus_values_transcript_len <- countsNA$ENSEMBL_GeneID

# Query for GO annotations
btaurus_ensembl_transcript_len <- getBM(attributes = btaurus_attributes_transcript_len, 
                                        filters = btaurus_filters_transcript_len, 
                                        values = btaurus_values_transcript_len, 
                                        mart = ensembl)

# go_terms <- readRDS(file = here::here("/home/kwanjau/tbrucei_gcn/data/intermediate/go_terms.RDS"))


#########################################################################################
# load GO terms associated with each parasite gene
#########################################################################################

# # load go terms from annotation package
# go_terms <- load_go_terms(orgdb, rownames(logcpm.norm.counts.combat),
#                           keytype='GID')
# 
# # this take time to run, so save it to avoid re-running.
# saveRDS(go_terms, file = here::here("data","intermediate","go_terms.RDS"))
go_terms <- readRDS(file = here::here("data","intermediate","go_terms.RDS"))

go_terms <- btaurus_ensembl_anno_GO

# Exclude genes not found in count table --not run--
#go_terms <- go_terms[go_terms$GID %in% rownames(logcpm.norm.counts.combat),]

# gene / go term mapping
gene_go_mapping <- as.data.frame(unique(go_terms %>% select(ensembl_gene_id, go_id, namespace_1003)))
colnames(gene_go_mapping) <- c('gene', 'category', 'ontology')

# go id / term mapping
go_term_id_mapping <- as.data.frame(unique(go_terms[c('go_id', 'namespace_1003', 'name_1006')]))
colnames(go_term_id_mapping) <- c("category", "term", "ontology")

#########################################################################################
# # Load KEGG annotations
# #########################################################################################
# 
# gene_kegg_mapping <- load_kegg_mapping(orgdb, rownames(logcpm.norm.counts.combat),
#                                        keytype="GID")
# 
# kegg_pathways <- load_kegg_pathways(orgdb, rownames(logcpm.norm.counts.combat),
#                                     keytype="GID")
# 
# 
# # Rename gene/KEGG mapping columns to be consistent with GO mapping
# colnames(gene_kegg_mapping) <- c('gene', 'category')
# colnames(kegg_pathways)     <- c('category', 'name', 'class', 'description')
# 
# kegg_pathways <- unique(kegg_pathways)
# 

################################################
# GO Enrichment
################################################

module_colours <- readRDS("module_colours.RDS")

# get the number of modules
num_modules <- length(unique(module_colours))

# Create gene lengths vector
gene_lengths <- btaurus_ensembl_transcript_len$transcript_length
names(gene_lengths) <- btaurus_ensembl_transcript_len$gene_id

# Gene IDs
log_counts <- readRDS("log_counts.RDS")
gene_ids <- rownames(log_counts)

# save the module sizes 
# Data frame of module sizes
module_counts <- c()
for (color in unique(module_colours)) {
  module_counts <- append(module_counts, sum(module_colours == color))
}

# create a mapping from module id to number of genes for later use
module_sizes <- data.frame(module_id=unique(module_colours),
                           num_genes=module_counts)

# Initialize parallelization
cl <- makeCluster(max(1, min(10, detectCores() - 2, na.rm = TRUE)))
registerDoParallel(cl)
message("Performing GO enrichment")

# Check each module for enrichment in GO terms and save result in a list
module_go_enrichment <- foreach(color=unique(module_colours), .packages=c('goseq')) %dopar% {
  set.seed(1)
  # Measure GO enrichment for module
  enriched <- tryCatch({
    # module gene ids
    in_module_geneids <- gene_ids[module_colours == color]
    message(sprintf("[GO enrichment] %s", color))
    
    # T. brucei GO enrichment
    enriched <- test_gene_enrichment(in_module_geneids, gene_ids,
                                     gene_go_mapping, gene_lengths)
    
    # Add descriptions
    enriched <- merge(enriched, go_term_id_mapping, by='category')
  }, error=function(e) {
    # goseq fails in some cases; have not been able to track down cause yet
    # to avoid errors we will just return an empty result set
    warning(sprintf("GO enrichment failed for module %s", color))
    cbind(
      get_enrichment_placeholder(),
      term=numeric(0),
      ontology=numeric(0)
    )
  })
  enriched
}
names(module_go_enrichment) <- unique(module_colours)

# remove any null/empty entries from the results
module_go_enrichment <- module_go_enrichment[!sapply(module_go_enrichment, is.null)]

# unregister cpus
stopCluster(cl)

# save the GO enrichment results
saveRDS(module_go_enrichment, file = here::here("data","intermediate","module_go_enrichment.RDS"))


#------------------------------------
# Print GO enrichment results
#------------------------------------
# temporarily repeat the gene / go term mapping to add 'term' column
gene_go_mapping_tmp <- as.data.frame(unique(go_terms %>% select(GID, GO, TERM, ONTOLOGY)))
colnames(gene_go_mapping_tmp) <- c('gene', 'category', 'term', 'ontology')

gene_info_tmp <- gene_info %>% select(-chromosome, -strand,
                                      - type, -transcript_length)
colnames(gene_info_tmp) <- c("gene","description","transcript_id")

tmp <- cbind(gene_info_tmp, color=module_colours)

#tmp <- cbind(gene=gene_ids, color=module_colours)
gene_mapping <- merge(gene_go_mapping_tmp, tmp, by='gene')
cat(sprintf('- Total enriched modules: %d\n', 
            sum(sapply(module_go_enrichment, nrow) > 0)))

# create tables of the results in this document
print_enrichment_results(module_go_enrichment, module_sizes, 'GO terms',
                         NULL, gene_mapping, 
                         output_dir=here::here("results","tables"),
                         enrichment_type='go',
                         include_gene_lists=FALSE)

enriched_colors_go <- get_enriched_modules(module_go_enrichment)

# Module enrichment status (used in dendrogram plots)
go_enrichment_status   <- as.numeric(module_colours %in% enriched_colors_go)


saveRDS(gene_mapping, file = here::here("data","intermediate","gene_mapping.RDS"))











































#########################################################################################

## wormbase for M. incognita
listMarts(host = "parasite.wormbase.org")

wormbase = useMart(biomart = "parasite_mart", host = "parasite.wormbase.org")
listDatasets(wormbase)

wormbase <- useDataset(mart = wormbase, dataset = "wbps_gene")
head(listFilters(wormbase))

wormbase_filters <- listFilters(wormbase)
wormbase_attributes <- listAttributes(wormbase)



reads_count_nehemiah <- read.csv("NA_new_nehemiah.txt", sep = "\t")


# Attributes and filters to be used in the bioMart query
mincgonita_attributes <- c("wbps_gene_id","wbps_transcript_id")
mincgonita_filters <- "wbps_gene_id"

# get the gene ids to be obtained from the query
mincgonita_values <- reads_count_nehemiah$ENSEMBL_GeneID

# Query
mincgonita_ensembl_anno <- getBM(attributes = mincgonita_attributes, 
                                 filters = mincgonita_filters, 
                                 values = mincgonita_values, 
                                 mart = wormbase)


