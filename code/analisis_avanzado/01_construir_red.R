
#!/usr/bin/env Rscript

# ==================================================
# ANALISIS AVANZADO DE REDES - BIOLOGIA DE SISTEMAS
# HPO: HP:0012347
# Parte 1: Construccion de la red
# ==================================================

library(STRINGdb)
library(igraph)

# 1. Parametros

HPO_ID <- "HP:0012347"
SCORE <- 800
STRING_VERSION <- "12.5"

set.seed(123)
options(timeout = 1200)

out_dir <- "results/analisis_avanzado"
string_dir <- "results/string_data"

dir.create(out_dir, recursive = TRUE,
           showWarnings = FALSE)

dir.create(string_dir, recursive = TRUE,
           showWarnings = FALSE)

# 2. Obtener genes de HPO

url <- paste0(
  "https://ontology.jax.org/api/network/annotation/",
  HPO_ID,
  "/download/gene"
)

hpo_file <- file.path(out_dir, "genes_hpo.tsv")

download.file(url, hpo_file, mode = "wb")

genes <- read.delim(hpo_file,
                    stringsAsFactors = FALSE)

genes$name <- trimws(genes$name)

genes <- unique(genes[
  !is.na(genes$name) & nzchar(genes$name),
  c("id", "name")
])

cat("Genes asociados al HPO:", nrow(genes), "\n")

# 3. Inicializar STRING

string_db <- STRINGdb$new(
  version = STRING_VERSION,
  species = 9606,
  score_threshold = SCORE,
  input_directory = string_dir
)

# 4. Mapear genes a STRING

string_ids <- string_db$map(
  genes,
  "name",
  removeUnmappedRows = TRUE
)

cat("Genes mapeados:",
    length(unique(string_ids$name)), "\n")

# 5. Obtener interacciones

inter <- string_db$get_interactions(
  unique(string_ids$STRING_id)
)

inter <- unique(inter[, 1:2])

inter$from <- string_ids$name[
  match(inter$from, string_ids$STRING_id)
]

inter$to <- string_ids$name[
  match(inter$to, string_ids$STRING_id)
]

stopifnot(!anyNA(inter))

# 6. Construir la red

g <- graph_from_data_frame(
  inter,
  directed = FALSE
)

g <- simplify(g)

# 7. Resultados iniciales

cat("\nRESULTADOS DE LA RED\n")
cat("Nodos:", vcount(g), "\n")
cat("Conexiones:", ecount(g), "\n")
cat("Grado medio:", mean(degree(g)), "\n")
cat("Densidad:", edge_density(g), "\n")
cat("Transitividad:", transitivity(g), "\n")
cat("Componentes conexas:", components(g)$no, "\n")

# 8. Guardar datos y red

write.csv(
  genes,
  file.path(out_dir, "genes_hpo.csv"),
  row.names = FALSE
)

write.csv(
  string_ids,
  file.path(out_dir, "genes_mapeados.csv"),
  row.names = FALSE
)

write.csv(
  as_data_frame(g, what = "edges"),
  file.path(out_dir, "red_aristas.csv"),
  row.names = FALSE
)

saveRDS(
  g,
  file.path(out_dir, "red_hpo.rds")
)

# 9. Visualizar red

pdf(file.path(out_dir, "red_original.pdf"),
    width = 10, height = 10)

plot(
  g,
  layout = layout_with_fr(g),
  vertex.color = "steelblue",
  vertex.size = 7,
  vertex.label.cex = 0.6,
  edge.color = "grey70",
  main = "HP:0012347 - STRING v12.5"
)

dev.off()

cat("\nRed guardada correctamente.\n")
