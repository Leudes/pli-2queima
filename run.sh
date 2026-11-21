#!/usr/bin/env bash
set -euo pipefail

PASTA_SCRIPT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

cd "$PASTA_SCRIPT"
julia --project=. burning2_pli_warmstart.jl input output --tempo 300 --threads 1
