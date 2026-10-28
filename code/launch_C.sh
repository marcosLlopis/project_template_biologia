#!/usr/bin/env bash
# Ejecuta el flujo completo. Todo lo que se descarga o genera va a ../results
set -euo pipefail
cd "$(dirname "$0")"
export R_LIBS_USER="$(pwd)/R_LIBS"

HPO="HP:0000098"
THRESHOLD=900
DATA=../results/data
mkdir -p "$DATA"

# Red de interacciones de STRING (humano) y filtro por puntuación
wget -nv "https://stringdb-downloads.org/download/protein.links.v12.0/9606.protein.links.v12.0.txt.gz" -O "$DATA/string.txt.gz"
gunzip -f "$DATA/string.txt.gz"
awk -F ' ' -v t="$THRESHOLD" 'NR==1 || $3 > t {print $1 "\t" $2 "\t" $3}' "$DATA/string.txt" > "$DATA/string_network.txt"

# Diccionario símbolo HGNC -> proteína STRING
wget -nv "https://stringdb-downloads.org/download/protein.aliases.v12.0/9606.protein.aliases.v12.0.txt.gz" -O "$DATA/aliases.txt.gz"
gunzip -f "$DATA/aliases.txt.gz"
grep 'Ensembl_HGNC_symbol' "$DATA/aliases.txt" | cut -f 1,2 > "$DATA/string_dict.txt"

# Genes asociados al término HPO
wget -nv "https://ontology.jax.org/api/network/annotation/${HPO}/download/gene" -O "$DATA/hpo_query.txt"
cut -f 2 "$DATA/hpo_query.txt" | tr -d ' ' > "$DATA/hpo_gene.txt"

Rscript ./analyse_net_C.R -i "$DATA/string_network.txt" \
  -g ../results/net.pdf \
  -G "$DATA/hpo_gene.txt" \
  -d "$DATA/string_dict.txt" \
  -s ../results/summary.txt \
  -t "$THRESHOLD" \
  -H "$HPO"

echo
echo "== Resumen para el report (results/summary.txt) =="
cat ../results/summary.txt
