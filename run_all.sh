#!/usr/bin/env bash
# Runs every experiment and stores its output under logs/. About an hour and a half on one core,
# two thirds of it in the sections steps and draws of the friction sweep.
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p logs

run() {
    local log="logs/$1.txt"
    shift
    echo "== python -m $*"
    python -u -m "$@" | tee "$log"
}

run table1_feasible_fraction        experiments.table1_feasible_fraction
run entry_register_d2               experiments.entry_register_d2
run table3_table4_conditioning      experiments.table3_table4_conditioning
run refinement_and_dimension        experiments.refinement_and_dimension
run figure1_feasible_slice          experiments.figure1_feasible_slice
run gauge_invariance                experiments.gauge_invariance
run table5_bit_depth                experiments.table5_bit_depth
run table6_friction_sweep           experiments.table6_friction_sweep
run table6_friction_sweep_dense     experiments.table6_friction_sweep dense
run table6_friction_sweep_steps     experiments.table6_friction_sweep steps
run table6_friction_sweep_periodic  experiments.table6_friction_sweep periodic
run table6_friction_sweep_draws     experiments.table6_friction_sweep draws
