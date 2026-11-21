#!/usr/bin/env bash
set -euo pipefail

PASTA_SCRIPT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PASTA_GRAFOS="${1:-$PASTA_SCRIPT/../BRKGA/graph_dataset/splits/final}"
PASTA_SAIDA="${2:-$PASTA_SCRIPT/output/base_final}"
PASTA_GRAFOS="$(realpath "$PASTA_GRAFOS")"
PASTA_SAIDA="$(realpath -m "$PASTA_SAIDA")"
TEMPO="${TEMPO:-300}"
THREADS="${THREADS:-1}"
MAX_BINARIAS="${MAX_BINARIAS:-1000000}"
JULIA_EXECUTAVEL="${JULIA_EXECUTAVEL:-$(command -v julia || true)}"

if [[ -z "$JULIA_EXECUTAVEL" && -x "$HOME/.juliaup/bin/julia" ]]; then
    JULIA_EXECUTAVEL="$HOME/.juliaup/bin/julia"
fi
if [[ -z "$JULIA_EXECUTAVEL" ]]; then
    echo "Julia nao encontrado. Instale com juliaup ou defina JULIA_EXECUTAVEL." >&2
    exit 1
fi

cd "$PASTA_SCRIPT"
"$JULIA_EXECUTAVEL" --project=. burning2_pli_warmstart.jl \
    "$PASTA_GRAFOS" \
    "$PASTA_SAIDA" \
    --tempo "$TEMPO" \
    --threads "$THREADS" \
    --max-binarias "$MAX_BINARIAS"
