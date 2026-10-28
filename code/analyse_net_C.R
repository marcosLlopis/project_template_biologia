#!/usr/bin/env Rscript
library(optparse)
library(igraph)

option_list <- list(
  make_option(c("-i", "--input_file"), type = "character", default = NULL,
              help = "Input file with graph edges"),
  make_option(c("-g", "--graph_file"), type = "character", default = NULL,
              help = "Output file with network graph"),
  make_option(c("-d", "--dict_file"), type = "character", default = NULL,
              help = "File to translate gene symbol to ensembl protein"),
  make_option(c("-G", "--gene_subset"), type = "character", default = NULL,
              help = "File with gene list associated to HPO"),
  make_option(c("-s", "--summary_file"), type = "character", default = NULL,
              help = "Optional output file with the numbers for the report"),
  make_option(c("-t", "--threshold"), type = "character", default = NA,
              help = "Optional score threshold (only written to the summary)"),
  make_option(c("-H", "--hpo"), type = "character", default = NA,
              help = "Optional HPO id (only written to the summary)")
)

opt <- parse_args(OptionParser(option_list = option_list))

net_table <- read.table(opt$input_file, sep = "\t", header = TRUE)
gene_list <- read.table(opt$gene_subset, sep = "\t", header = TRUE)
dict_gene <- read.table(opt$dict_file, sep = "\t", header = FALSE)

n_hpo_genes <- length(unique(gene_list$name))

match_vector <- match(gene_list$name, dict_gene$V2)
gene_list <- dict_gene$V1[match_vector]
gene_list <- gene_list[!is.na(gene_list)]

net <- graph_from_data_frame(d = net_table)
nodes_in_graph <- match(gene_list, V(net)$name)
nodes_in_graph <- nodes_in_graph[!is.na(nodes_in_graph)]
gene_list <- V(net)$name[nodes_in_graph]

net <- induced_subgraph(net, gene_list)
print(net)  # primera línea: número de nodos y de aristas
trans <- transitivity(net)
print(trans)

pdf(file = opt$graph_file)
plot(net, vertex.color = "orange", vertex.size = 5, vertex.label = NA,
     edge.width = 0.5, main = "Nicely Layout")
invisible(dev.off())

if (!is.null(opt$summary_file)) {
  lines <- c(
    paste0("HPO: ", opt$hpo),
    paste0("X (genes del HPO): ", n_hpo_genes),
    paste0("Y (umbral de corte): ", opt$threshold),
    paste0("W (genes en la red): ", vcount(net)),
    paste0("Z (conexiones): ", ecount(net)),
    paste0("J (transitivity): ", round(trans, 4))
  )
  writeLines(lines, opt$summary_file)
}
