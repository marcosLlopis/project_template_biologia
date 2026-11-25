#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(STRINGdb)
  library(igraph)
})

HPO_ID <- "HP:0012347"
STRING_VERSION <- "12.5"
SCORES <- c(600, 800, 900)
DATA_DIR <- "results/string_data"
OUT_DIR <- "results/analisis_avanzado"
GENES_FILE <- file.path(OUT_DIR, "genes_hpo.tsv")

if (!file.exists(GENES_FILE)) {
  stop("Falta genes_hpo.tsv. Ejecuta primero 01_construir_red.R desde la raíz.")
}

genes <- read.delim(GENES_FILE, stringsAsFactors = FALSE)
genes$name <- trimws(genes$name)
genes <- unique(genes[
  !is.na(genes$name) & nzchar(genes$name),
  c("id", "name")
])

resumen <- vector("list", length(SCORES))

for (i in seq_along(SCORES)) {
  score <- SCORES[i]
  message("Calculando red con umbral STRING ", score)

  db <- STRINGdb$new(
    version = STRING_VERSION,
    species = 9606,
    score_threshold = score,
    input_directory = DATA_DIR
  )

  mapped <- db$map(genes, "name", removeUnmappedRows = TRUE)
  inter <- db$get_interactions(unique(mapped$STRING_id))

  if (nrow(inter) > 0) {
    inter <- unique(inter[, c("from", "to")])
    inter$from <- mapped$name[match(inter$from, mapped$STRING_id)]
    inter$to <- mapped$name[match(inter$to, mapped$STRING_id)]
    inter <- inter[complete.cases(inter), , drop = FALSE]
  } else {
    inter <- data.frame(from = character(), to = character())
  }

  vertices <- data.frame(name = unique(mapped$name))
  g_completo <- graph_from_data_frame(
    inter,
    vertices = vertices,
    directed = FALSE
  )
  g_completo <- simplify(
    g_completo,
    remove.multiple = TRUE,
    remove.loops = TRUE
  )

  activos <- V(g_completo)[degree(g_completo) > 0]
  g <- induced_subgraph(g_completo, vids = activos)

  if (vcount(g) > 0) {
    componentes <- components(g)$no
    grado_medio <- mean(degree(g))
    densidad <- edge_density(g)
    transitividad <- if (vcount(g) >= 3) {
      transitivity(g, type = "global")
    } else {
      NA_real_
    }

    if (ecount(g) > 0) {
      set.seed(123)
      comunidades <- cluster_louvain(g)
      n_comunidades <- length(comunidades)
      modularidad <- modularity(comunidades)
    } else {
      n_comunidades <- NA_integer_
      modularidad <- NA_real_
    }
  } else {
    componentes <- 0
    grado_medio <- NA_real_
    densidad <- NA_real_
    transitividad <- NA_real_
    n_comunidades <- NA_integer_
    modularidad <- NA_real_
  }

  resumen[[i]] <- data.frame(
    HPO = HPO_ID,
    umbral = score,
    genes_mapeados = length(unique(mapped$name)),
    genes_sin_interaccion = length(unique(mapped$name)) - vcount(g),
    nodos = vcount(g),
    aristas = ecount(g),
    grado_medio = grado_medio,
    densidad = densidad,
    transitividad = transitividad,
    componentes = componentes,
    comunidades = n_comunidades,
    modularidad = modularidad
  )
}

tabla <- do.call(rbind, resumen)
write.csv(
  tabla,
  file.path(OUT_DIR, "sensibilidad_umbral_STRING.csv"),
  row.names = FALSE
)
print(tabla, row.names = FALSE)

base <- tabla[tabla$umbral == 800, ]
if (nrow(base) == 1) {
  cat(
    "\nComprobación del umbral 800: ",
    base$nodos, " nodos y ", base$aristas, " aristas.\n",
    sep = ""
  )
}

pdf(file.path(OUT_DIR, "sensibilidad_umbral_STRING.pdf"), width = 9, height = 7)
par(mfrow = c(2, 2))

plot(tabla$umbral, tabla$nodos, type = "b", pch = 19,
     xlab = "Umbral STRING", ylab = "Nodos con interacción",
     main = "Tamaño de la red")
plot(tabla$umbral, tabla$aristas, type = "b", pch = 19,
     xlab = "Umbral STRING", ylab = "Aristas",
     main = "Conexiones")
plot(tabla$umbral, tabla$transitividad, type = "b", pch = 19,
     xlab = "Umbral STRING", ylab = "Transitividad global",
     main = "Transitividad")
plot(tabla$umbral, tabla$modularidad, type = "b", pch = 19,
     xlab = "Umbral STRING", ylab = "Modularidad",
     main = "Comunidades de Louvain")

dev.off()
cat("\nTabla y figura guardadas en ", OUT_DIR, ".\n", sep = "")
