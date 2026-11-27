#!/usr/bin/env bash
# Reproduce every result, table and figure of the paper.
#
#   ./run_all.sh          regenerate figures and tables from the stored results (about one minute)
#   ./run_all.sh --full   rerun all optimizations from scratch first (about 30 minutes)
set -euo pipefail
cd "$(dirname "$0")"

if [[ "${1:-}" == "--full" ]]; then
  rm -f results/rsm_iae_history.json results/rsm_iae_history_staticinit.json \
        results/bench2_*.json
  python run_rsm_iae.py results/rsm_iae_history.json 200   # SAFEHX, initial point of Table 1
  python run_rsm_iae_static_init.py                        # SAFEHX, static-decoupler start
  python benchmarks2.py Hylte                              # benchmark processes, two-way
  python benchmarks2.py WoodBerry
  python optimize_static.py                                # static inverted decoupler, optimized gains
  python load_disturbance_benchmarks.py                    # load disturbances, benchmarks
fi

python make_figures.py    # figures/ and results/*_table.json, montecarlo.json, rga.json
python gen_macros.py      # results/results_macros.tex (every number quoted in the paper)
echo "Done: see figures/ and results/results_macros.tex"
