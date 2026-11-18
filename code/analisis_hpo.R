#!/usr/bin/env Rscript
# ============================================================================
# Flujo de análisis mínimo para un término HPO (Biología de Sistemas, UMA)
# Adaptado de la guía "Flujo de análisis mínimo para un término HPO".
#
# Uso (desde la carpeta del repositorio):
#   Rscript code/analisis_hpo.R                 # usa HPO_ID por defecto
#   Rscript code/analisis_hpo.R HP:0000098      # otro término HPO
#
# Salidas (carpeta results/analisis_<HPO>/):
#   genes_hpo.csv, red_aristas.csv, metricas.txt, resumen_clusters.csv,
#   enriquecimiento_global.csv, enriquecimiento_cluster_<n>.csv,
#   red_basica.png, red_grado.png, red_clusters.png,
#   dotplot_global.png, dotplot_cluster_<n>.png
# ============================================================================

# ------------------------- Parámetros -------------------------------------
args    <- commandArgs(trailingOnly = TRUE)
HPO_ID  <- if (length(args) >= 1) args[1] else "HP:0000098"
SCORE   <- 800      # umbral de puntuación de STRING (0-1000)
MIN_SZ  <- 5        # tamaño mínimo de cluster para hacer enriquecimiento
STRING_VERSION <- "12.5"
SEED    <- 123

if (!grepl("^HP:[0-9]{7}$", HPO_ID)) stop("HPO no válido: ", HPO_ID)

# La salida va a results/ del repositorio (se detecta desde code/ o desde la raíz)
base <- if (basename(getwd()) == "code") ".." else "."
out_dir <- file.path(base, "results", paste0("analisis_", sub(":", "", HPO_ID)))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
string_dir <- file.path(base, "results", "string_data")   # caché de STRING
dir.create(string_dir, recursive = TRUE, showWarnings = FALSE)
out <- function(f) file.path(out_dir, f)

# ------------------------- Paquetes ---------------------------------------
instalar_cran <- function(pkgs) {
  for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
}
instalar_bioc <- function(pkgs) {
  instalar_cran("BiocManager")
  for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, ask = FALSE, update = FALSE)
}
instalar_cran(c("httr", "jsonlite", "igraph"))
instalar_bioc(c("STRINGdb", "clusterProfiler", "org.Hs.eg.db", "enrichplot"))

suppressPackageStartupMessages({
  library(httr); library(jsonlite)
  library(STRINGdb); library(clusterProfiler); library(org.Hs.eg.db); library(enrichplot)
  library(igraph)
})

png_plot <- function(file, expr, w = 1600, h = 1400) {
  png(file, width = w, height = h, res = 150); on.exit(dev.off()); force(expr)
}

# ------------------------- Parte 1: genes de la HPO -----------------------
message("== Parte 1: genes de ", HPO_ID)
resp <- httr::GET(paste0("https://ontology.jax.org/api/network/annotation/", HPO_ID))
httr::stop_for_status(resp)
parsed <- jsonlite::fromJSON(httr::content(resp, as = "text", encoding = "UTF-8"))
genes  <- parsed$genes
if (is.null(genes) || nrow(genes) == 0) stop("La API no devolvió genes para ", HPO_ID)
genes <- unique(genes[, c("id", "name")])
nombre_hpo <- if (length(args) >= 2) args[2] else NA_character_
write.csv(genes, out("genes_hpo.csv"), row.names = FALSE)
cat("Genes en la HPO:", nrow(genes), "\n")

# ------------------------- Parte 2: red de STRING -------------------------
message("== Parte 2: STRING (versión ", STRING_VERSION, ", umbral ", SCORE, ")")
string_db  <- STRINGdb$new(version = STRING_VERSION, species = 9606,
                           score_threshold = SCORE, input_directory = string_dir)
string_ids <- string_db$map(genes, "name", removeUnmappedRows = TRUE)
no_mapeados <- setdiff(genes$name, string_ids$name)
cat("Genes mapeados a STRING:", nrow(string_ids), "\n")
cat("Genes no mapeados:", if (length(no_mapeados)) paste(no_mapeados, collapse = " ") else "(ninguno)", "\n")

inter <- string_db$get_interactions(string_ids$STRING_id)
inter <- unique(inter[, 1:2])
inter$from <- string_ids$name[match(inter$from, string_ids$STRING_id)]
inter$to   <- string_ids$name[match(inter$to,   string_ids$STRING_id)]
write.csv(inter, out("red_aristas.csv"), row.names = FALSE)

# ------------------------- Parte 3: igraph y propiedades ------------------
message("== Parte 3: igraph")
g <- graph_from_data_frame(inter, directed = FALSE)
g <- igraph::simplify(g)                       # quita bucles y aristas duplicadas (A-B / B-A)
sin_interacciones <- setdiff(string_ids$name, V(g)$name)

n_comp    <- igraph::components(g)$no
gcc       <- largest_component(g)
armonica  <- harmonic_centrality(g, normalized = TRUE)
top_arm   <- head(sort(armonica, decreasing = TRUE), 10)

metricas <- c(
  paste("HPO:", HPO_ID, if (!is.na(nombre_hpo)) paste0("(", nombre_hpo, ")") else ""),
  paste("Genes en la HPO:", nrow(genes)),
  paste("Genes mapeados a STRING:", nrow(string_ids)),
  paste("Genes sin interacciones (fuera de la red):", length(sin_interacciones)),
  paste("Umbral STRING:", SCORE),
  paste("Nodos:", gorder(g)),
  paste("Aristas:", gsize(g)),
  paste("Grado medio:", round(mean(degree(g)), 3)),
  paste("Densidad:", round(edge_density(g), 5)),
  paste("Transitivity (global):", round(transitivity(g, type = "global"), 4)),
  paste("Componentes conexas:", n_comp),
  paste("Nodos en la componente mayor:", gorder(gcc)),
  "Top 10 centralidad armónica (normalizada):",
  paste0("  ", names(top_arm), ": ", round(top_arm, 3))
)
writeLines(metricas, out("metricas.txt"))
cat(paste(metricas, collapse = "\n"), "\n")

# ------------------------- Parte 4: visualización -------------------------
message("== Parte 4: figuras")
set.seed(SEED); lay <- layout_with_fr(g)
png_plot(out("red_basica.png"),
         plot(g, layout = lay, vertex.color = "orange", vertex.size = 5,
              vertex.label = NA, edge.width = 0.5, main = paste("Red STRING -", HPO_ID)))
png_plot(out("red_grado.png"),
         plot(g, layout = lay, vertex.color = "skyblue", vertex.size = 3 + degree(g) / 2,
              vertex.label = NA, edge.width = 0.5, main = "Tamaño del nodo según el grado"))

# ------------------------- Parte 5: clustering ----------------------------
message("== Parte 5: clustering (Louvain)")
set.seed(SEED)
comunidades <- cluster_louvain(g)
cat("Número de clusters:", length(comunidades), "\n")
print(sizes(comunidades))
png_plot(out("red_clusters.png"),
         plot(comunidades, g, layout = lay, vertex.size = 5, vertex.label = NA,
              main = "Clustering de Louvain"))

# ------------------------- Parte 6: enriquecimiento -----------------------
message("== Parte 6: enriquecimiento GO (BP)")
# Global: genes de la HPO frente a todo el genoma anotado
enr_global <- enrichGO(gene = genes$name, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                       ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
write.csv2(as.data.frame(enr_global), out("enriquecimiento_global.csv"), row.names = FALSE)
if (nrow(as.data.frame(enr_global)) > 0) {
  png_plot(out("dotplot_global.png"),
           print(dotplot(enr_global, showCategory = 10, title = "GO BP (HPO completa)")),
           w = 1800, h = 1400)
}

# Por cluster: fondo = genes de la red
genes_clusters <- split(V(g)$name, membership(comunidades))
genes_clusters <- genes_clusters[lengths(genes_clusters) >= MIN_SZ]
cat("Clusters con al menos", MIN_SZ, "genes:", length(genes_clusters), "\n")

enr_clusters <- list()
for (i in names(genes_clusters)) {
  enr_clusters[[i]] <- enrichGO(gene = genes_clusters[[i]], universe = V(g)$name,
                                OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
                                pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
  df <- as.data.frame(enr_clusters[[i]])
  write.csv2(df, out(paste0("enriquecimiento_cluster_", i, ".csv")), row.names = FALSE)
  if (nrow(df) > 0) {
    png_plot(out(paste0("dotplot_cluster_", i, ".png")),
             print(dotplot(enr_clusters[[i]], showCategory = 10,
                           title = paste("GO BP, cluster", i))), w = 1800, h = 1400)
  }
}

resumen <- data.frame(
  cluster  = names(genes_clusters),
  genes    = lengths(genes_clusters),
  top_term = vapply(enr_clusters, function(x) {
    r <- as.data.frame(x); if (nrow(r) == 0) "(ninguno)" else r$Description[1]
  }, character(1)),
  row.names = NULL
)
write.csv2(resumen, out("resumen_clusters.csv"), row.names = FALSE)
print(resumen)

message("\nListo. Resultados en: ", normalizePath(out_dir))
