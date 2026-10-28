#!/usr/bin/env bash
# Instala las dependencias de R en una carpeta local (sin permisos de administrador)
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p R_LIBS
# Mismo entorno que launch: oculta la biblioteca personal del usuario
export R_LIBS_USER="$(pwd)/R_LIBS"
Rscript -e "install.packages(c('rlang', 'igraph', 'optparse'), lib='./R_LIBS', repos='https://cloud.r-project.org')"
echo "Dependencias instaladas en $(pwd)/R_LIBS"
