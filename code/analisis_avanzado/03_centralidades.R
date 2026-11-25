
#!/usr/bin/env Rscript

# ==================================================
# ANALISIS AVANZADO - HP:0012347
# Centralidades y conectividad entre comunidades
# ==================================================

library(igraph)

set.seed(123)

out_dir <- "results/analisis_avanzado"

g <- readRDS(file.path(out_dir, "red_hpo.rds"))

# 1. Centralidades globales

grado <- degree(g)

intermediacion <- betweenness(
  g,
  directed = FALSE,
  normalized = TRUE
)

armonica <- harmonic_centrality(
  g,
  normalized = TRUE
)

pagerank <- page_rank(g)$vector

# 2. Componente conexa mayor

gcc <- largest_component(g)

cercania <- closeness(
  gcc,
  normalized = TRUE
)

autovector <- eigen_centrality(
  gcc,
  directed = FALSE
)$vector

# 3. Comunidades originales de la memoria
# Usamos la particion descrita por el grupo

clusters <- list(
  c("ALG1", "ALG2", "ALG3", "ALG6", "ALG8",
    "ALG9", "ALG11", "ALG12", "ALG13", "PMM2", "RFT1"),

  c("DPAGT1", "DHRSX", "DOLK", "DPM1",
    "DPM2", "DPM3", "MPDU1", "SRD5A3"),

  c("ATP6AP1", "ATP6AP2", "ATP6V0A2",
    "ATP6V1A", "ATP6V1E1", "VMA12"),

  c("COG1", "COG3", "COG4", "COG5",
    "COG6", "COG7", "COG8"),

  c("GET4", "CAMLG"),

  c("MPI", "PGM1", "PGM2L1"),

  c("DDOST", "EDEM3", "MAGT1", "SSR4")
)

genes_clusters <- unlist(clusters)

stopifnot(
  !anyDuplicated(genes_clusters),
  setequal(genes_clusters, V(g)$name)
)

comunidad <- setNames(
  rep(seq_along(clusters), lengths(clusters)),
  genes_clusters
)

# 4. Coeficiente de participacion

participacion <- sapply(V(g)$name, function(gen) {

  vecinos <- as_ids(neighbors(g, gen))

  if (length(vecinos) == 0) return(0)

  frecuencias <- table(comunidad[vecinos])

  1 - sum((frecuencias / length(vecinos))^2)
})

# 5. Conexiones con otras comunidades

conexiones_externas <- sapply(
  V(g)$name,
  function(gen) {

    vecinos <- as_ids(neighbors(g, gen))

    sum(comunidad[vecinos] != comunidad[gen])
  }
)

# 6. Puntos de articulacion

articulaciones <- as_ids(
  articulation_points(g)
)

# 7. Tabla completa

tabla <- data.frame(
  Gen = V(g)$name,
  Cluster = as.integer(comunidad[V(g)$name]),
  Grado = as.numeric(grado),
  Intermediacion = as.numeric(intermediacion),
  Armonica = as.numeric(armonica),
  PageRank = as.numeric(pagerank),
  Cercania = as.numeric(cercania[V(g)$name]),
  Autovector = as.numeric(autovector[V(g)$name]),
  Participacion = as.numeric(participacion),
  Conexiones_externas = as.numeric(conexiones_externas),
  Articulacion = V(g)$name %in% articulaciones
)

tabla <- tabla[
  order(-tabla$Intermediacion),
]

write.csv(
  tabla,
  file.path(out_dir, "centralidades.csv"),
  row.names = FALSE
)

# 8. Rankings

medidas <- c(
  "Grado",
  "Intermediacion",
  "Armonica",
  "PageRank",
  "Cercania",
  "Autovector",
  "Participacion",
  "Conexiones_externas"
)

for (m in medidas) {

  cat("\nTOP 10 -", m, "\n")

  orden <- order(
    tabla[[m]],
    decreasing = TRUE,
    na.last = NA
  )

  print(
    head(tabla[orden, c("Gen", "Cluster", m)], 10),
    row.names = FALSE
  )
}

cat("\nPUNTOS DE ARTICULACION\n")
print(articulaciones)

cat("\nMODULARIDAD DE LA PARTICION\n")
print(modularity(
  g,
  as.integer(comunidad[V(g)$name])
))

# 9. Graficas de centralidad

pdf(
  file.path(out_dir, "ranking_centralidades.pdf"),
  width = 10,
  height = 7
)

par(mar = c(5, 8, 4, 2))

for (m in medidas) {

  orden <- order(
    tabla[[m]],
    decreasing = TRUE,
    na.last = NA
  )

  top <- head(tabla[orden, ], 10)

  barplot(
    rev(top[[m]]),
    names.arg = rev(top$Gen),
    horiz = TRUE,
    las = 1,
    col = "steelblue",
    main = paste("Top 10 -", m),
    xlab = m
  )
}

dev.off()

# 10. Correlaciones de Spearman

correlaciones <- cor(
  tabla[, medidas],
  method = "spearman",
  use = "pairwise.complete.obs"
)

write.csv(
  correlaciones,
  file.path(out_dir, "correlaciones_centralidad.csv")
)

cat("\nCORRELACIONES DE SPEARMAN\n")
print(round(correlaciones, 3))

pdf(
  file.path(out_dir, "correlaciones_centralidad.pdf"),
  width = 10,
  height = 10
)

heatmap(
  correlaciones,
  Rowv = NA,
  Colv = NA,
  scale = "none",
  margins = c(12, 12),
  col = colorRampPalette(
    c("steelblue", "white", "firebrick")
  )(100),
  breaks = seq(-1, 1, length.out = 101),
  main = "Correlaciones de Spearman"
)

dev.off()

# 11. Red coloreada por comunidades
# Tamano proporcional a la intermediacion

colores <- c(
  "steelblue",
  "darkorange",
  "seagreen",
  "orchid",
  "gold",
  "tomato",
  "grey60"
)

tamano <- 5 + 18 * intermediacion /
  max(intermediacion)

pdf(
  file.path(out_dir, "red_intermediacion.pdf"),
  width = 12,
  height = 10
)

plot(
  g,
  layout = layout_with_fr(g),
  vertex.color = colores[comunidad[V(g)$name]],
  vertex.size = tamano,
  vertex.label.cex = 0.65,
  vertex.frame.color = ifelse(
    V(g)$name %in% articulaciones,
    "red", "grey40"
  ),
  edge.color = "grey75",
  main = "Comunidades e intermediacion - HP:0012347"
)

legend(
  "bottomleft",
  legend = paste("Cluster", 1:7),
  col = colores,
  pch = 19,
  bty = "n"
)

dev.off()

# 12. Grado frente a intermediacion

pdf(
  file.path(out_dir, "grado_intermediacion.pdf"),
  width = 10,
  height = 7
)

plot(
  tabla$Grado,
  tabla$Intermediacion,
  pch = 19,
  col = colores[tabla$Cluster],
  xlab = "Grado",
  ylab = "Intermediacion normalizada",
  main = "Grado frente a intermediacion"
)

destacados <- unique(c(
  head(order(-tabla$Intermediacion), 5),
  head(order(-tabla$Grado), 2)
))

text(
  tabla$Grado[destacados],
  tabla$Intermediacion[destacados],
  labels = tabla$Gen[destacados],
  pos = 3,
  cex = 0.8
)

dev.off()

cat("\nAnalisis de centralidades y comunidades completado.\n")
