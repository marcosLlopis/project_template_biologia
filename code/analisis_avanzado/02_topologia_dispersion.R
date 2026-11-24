#!/usr/bin/env Rscript

# ==================================================
# ANALISIS AVANZADO - HP:0012347
# Parte 2: Topologia y dispersion del grado
# ==================================================

library(igraph)

out_dir <- "results/analisis_avanzado"

# 1. Cargar la red previamente construida

g <- readRDS(file.path(out_dir, "red_hpo.rds"))

# 2. Obtener el grado de cada nodo

grados <- degree(g)

# 3. Calcular medidas descriptivas

metricas <- data.frame(
  Medida = c(
    "Numero de nodos",
    "Numero de conexiones",
    "Densidad",
    "Transitividad",
    "Grado medio",
    "Mediana",
    "Varianza",
    "Desviacion tipica",
    "Coeficiente de variacion (%)",
    "Grado minimo",
    "Primer cuartil",
    "Tercer cuartil",
    "Grado maximo",
    "Rango intercuartilico"
  ),
  Valor = c(
    vcount(g),
    ecount(g),
    edge_density(g),
    transitivity(g),
    mean(grados),
    median(grados),
    var(grados),
    sd(grados),
    sd(grados) / mean(grados) * 100,
    min(grados),
    quantile(grados, 0.25),
    quantile(grados, 0.75),
    max(grados),
    IQR(grados)
  )
)

print(metricas, row.names = FALSE)

# 4. Guardar las medidas

write.csv(
  metricas,
  file.path(out_dir, "metricas_dispersion.csv"),
  row.names = FALSE
)

# 5. Crear tabla de genes ordenados por grado

tabla_grados <- data.frame(
  Gen = names(grados),
  Grado = as.numeric(grados)
)

tabla_grados <- tabla_grados[
  order(-tabla_grados$Grado),
]

write.csv(
  tabla_grados,
  file.path(out_dir, "grados_genes.csv"),
  row.names = FALSE
)

cat("\nTOP 10 GENES CON MAYOR GRADO\n")
print(head(tabla_grados, 10), row.names = FALSE)

# 6. Distribucion del grado

pdf(
  file.path(out_dir, "distribucion_grado.pdf"),
  width = 9,
  height = 6
)

hist(
  grados,
  breaks = seq(-0.5, max(grados) + 0.5, by = 1),
  col = "steelblue",
  border = "white",
  main = "Distribucion del grado - HP:0012347",
  xlab = "Numero de conexiones",
  ylab = "Frecuencia"
)

abline(
  v = mean(grados),
  col = "red",
  lwd = 2,
  lty = 2
)

legend(
  "topright",
  legend = "Grado medio",
  col = "red",
  lty = 2,
  lwd = 2
)

dev.off()

# 7. Diagrama de cajas

pdf(
  file.path(out_dir, "boxplot_grado.pdf"),
  width = 7,
  height = 6
)

boxplot(
  grados,
  col = "lightblue",
  main = "Dispersion del grado - HP:0012347",
  ylab = "Numero de conexiones"
)

dev.off()

# 8. Grafico de los 10 genes mas conectados

top10 <- head(tabla_grados, 10)

pdf(
  file.path(out_dir, "top10_grado.pdf"),
  width = 10,
  height = 7
)

barplot(
  rev(top10$Grado),
  names.arg = rev(top10$Gen),
  horiz = TRUE,
  las = 1,
  col = "steelblue",
  xlab = "Numero de conexiones",
  main = "Top 10 genes por grado"
)

dev.off()

cat("\nAnalisis de dispersion completado.\n")