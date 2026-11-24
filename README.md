# Attacks on the configuration adaptation rule: simulation code

Code and results for the paper "Атаки на правило адаптації конфігурації платформи: межі часу порушення та нижня оцінка вартості захисту" (R. Tsarov, S. Staikutsa, H. Tsarova).

## Files
- `model.py` – platform model: cost (1), admissibility (2)–(3), optimal configuration x*(λ), admissibility bounds.
- `sim.py` – one-day simulation with 1 s step (Numba): rules R1–R5, PR (Holt + residual margin, optional reserve), HPA-like; reconfiguration time d; periodic pulses, boundary-tracking attacker, bias.
- `trace.py` – loads the per-minute WorldCup98 counts and scales each day to a mean of 1450 req/s.
- `experiments.py` – all main experiments (tables 1–3, 5, period sweep, lower bound) -> `results.json`.
- `extra_runs.py` – PR+reserve and the reserve cost table -> `results_pr5.json`, `results_t4.json`.
- `figs.py` – figures 1–5.
- `validate.py` – comparison with the earlier version of the model.
- `sim2.py` – extended simulator: per-pulse records of initial violations (Proposition 2) and the attacker with a noisy margin estimate (`amode=3`).
- `exp_r4.py` – per-pulse check of Proposition 2, attacker with estimation error, lower margin ρ⁻ -> `results_r4.json`.
- `experiments_full.py` – same as `experiments.py`, stores all raw per-replication values -> `results_full.json`.
- `stats_r4.py` – block bootstrap (blocks of 7 days) for the trace and for paired differences -> `stats_r4.json`; `raw_extra.json`, `stats_extra.json` – the same for PR+reserve and the cost table.
- `verify_p2.py` – check of Proposition 2 on pulses that meet all its assumptions -> `verify_p2.json`.
- `fig2.py`, `figs_r4.py` – current figures 2, 3 and 5.

## Data
Per-minute request counts of the 1998 World Cup Web site trace: https://github.com/nimamahmoudi/worldcup98-dataset (file `invocation_count.csv`). Clone it so that `../trace/wc98/invocation_count.csv` exists relative to this folder, or edit `SRC` in `trace.py`.

## Requirements
Python 3.10+, numpy, scipy, pandas, numba, matplotlib:

    pip install -r requirements.txt

## Run
    python experiments.py      # ~1 min
    python extra_runs.py
    python figs.py
    python experiments_full.py
    python exp_r4.py
    python stats_r4.py
    python verify_p2.py
    python fig2.py
    python figs_r4.py

## Licence
Code: MIT (see `LICENSE`). The WorldCup98 per-minute counts are distributed by their repository under its own MIT licence and are not included here.

## Citation
If you use this code, please cite the paper (reference will be added after publication).
