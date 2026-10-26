# EMFC

**Engagement Maximization via Forwarding Control**

The `./EMFC` program solves EMFC on a graph with node hyperparameters \((\lambda_v, q_v, \tau_v)\): it estimates expected engagement and optimizes forwarding probabilities \(\xi_v \in [\xi_{\mathrm lo}, 1]\) using RR-graph sampling and coordinate-ascent sweeps, then evaluates solutions by Monte Carlo.

Example setup: WC propagation, \(\lambda \sim \mathrm{Uniform}[0,0.5]\) (`lam_unirand_0_0.5`), \(\xi_{\mathrm lo}=0.5\), \(\varepsilon=0.5\).

---

## Build

Requires **C++14**:

```bash
make EMFC
```

All EMFC modes below are invoked as `./EMFC -func=…`.

---

## Inputs

**Graph (WC):** \(p(u\to v)=1/\mathrm{indeg}(v)\).

```bash
./EMFC -func=format -gname=soc-facebook -dir=graphInfo -pdist=wc
```

Creates `graphInfo/{gname}.graph` and `graphInfo/{gname}.nodehyper.vec`.

**Hyperparameters for a run** (does not overwrite the baseline `.nodehyper.vec`):

```bash
./EMFC -func=hyp_derive -gname=soc-facebook -dir=graphInfo \
  -hyp_profile=graphInfo/hyp_profiles/lam_unirand_0_0.5.profile -rand_seed=1
```

Writes `graphInfo/{gname}.nodehyper.lam_unirand_0_0.5.vec`; pass `-nodehyper_suffix=lam_unirand_0_0.5` on solve/eval.

---

## `-func` (EMFC)

| `-func` | Role |
|---------|------|
| `frim_rr_graph` | Solve: RR-graph + FRIM sweeps → save \(\xi\) |
| `frim_rr_graph_prune` | Solve: K/L pruning, then the same sweeps |
| `eval` | MC total benefit for \(\xi\) from a prior solve (or `-eval_xi`) |
| `format` | Prepare graph and baseline node hyperparameters |
| `hyp_derive` | Resample node hyperparameters from a profile |

---

## RR sample count \(R\)

- **`-eps`**: accuracy \(\varepsilon\).
- Without **`-frim_rr`**, RR-graph solves use

\[
R = \left\lceil \frac{(8+2\varepsilon)\,\bigl(\log(1/\delta) + n\log 2 + \log 2\bigr)}{\varepsilon^2} \right\rceil,
\]

\(n=|V|\), \(\delta = (\texttt{-delta}\ \text{or}\ 1/n)\,/\,2\).

- **`-frim_rr=N`**: fixed \(R\); use the same value on matching eval.

Also: **`-delta`**, **`-frim_sweeps`** (default 10), **`-xi_lo`**, **`-xi_init=lo|1`**.

---

## Solve

Outputs under **`-outpath`**: `info/` (times, \(\hat J\)), `xi/` (\(\xi\)).

```bash
OUT=./result/example
mkdir -p "$OUT"

./EMFC -func=frim_rr_graph \
  -gname=soc-facebook -dir=graphInfo -pdist=wc \
  -outpath="$OUT" \
  -rand_seed=1 \
  -eps=0.5 -xi_lo=0.5 \
  -frim_sweeps=10 \
  -nodehyper_suffix=lam_unirand_0_0.5 \
  -rr_graph_gate_sweep_index=1 \
  -rr_graph_store_hit_only=1
```

Pruned solve:

```bash
./EMFC -func=frim_rr_graph_prune \
  -gname=soc-facebook -dir=graphInfo -pdist=wc \
  -outpath="$OUT" \
  -rand_seed=1 \
  -eps=0.5 -xi_lo=0.5 \
  -frim_sweeps=10 \
  -nodehyper_suffix=lam_unirand_0_0.5 \
  -rr_graph_gate_sweep_index=1 \
  -rr_graph_store_hit_only=1
```

| Flag | Meaning |
|------|---------|
| `-gname`, `-dir=graphInfo`, `-pdist=wc` | Graph and WC model |
| `-outpath` | Root for `info/`, `xi/`, `evaluation/` |
| `-rand_seed` | Seed (same on eval) |
| `-nodehyper_suffix` | Which `*.nodehyper.{suffix}.vec` to load |
| `-rr_graph_gate_sweep_index=1` | Skip sweep steps that do not change RR weights |
| `-rr_graph_store_hit_only=1` | Compact RR storage |

Prune logs may show `active=…/num_v`; timings are in `info/` (`time_sample_sec`, `time_total_sec`).

---

## Eval

Run after solve with the **same** graph, `-outpath`, seed, `-eps`, `-frim_sweeps`, `-nodehyper_suffix`, RR flags, and **`-method=`** equal to the solve `-func`.

```bash
./EMFC -func=eval \
  -gname=soc-facebook -dir=graphInfo -pdist=wc \
  -outpath="$OUT" \
  -rand_seed=1 \
  -method=frim_rr_graph \
  -eps=0.5 \
  -frim_sweeps=10 \
  -nodehyper_suffix=lam_unirand_0_0.5 \
  -eval_mc=10000 \
  -rr_graph_gate_sweep_index=1 \
  -rr_graph_store_hit_only=1
```

After `frim_rr_graph_prune` solve, set `-method=frim_rr_graph_prune`.

| Flag | Meaning |
|------|---------|
| `-method` | Solver that wrote \(\xi\) (required if not using `-eval_xi`) |
| `-eval_mc` | MC replications (default 10000) |
| `-eval_xi=1` / `lo` / `half` / `unirand` / number | MC at a fixed or random \(\xi\) without loading solve output |

Results: `{outpath}/evaluation/eval_{gname}_{seed}_…`.

To print the expected `info`/`xi` basename before running:

```bash
./EMFC -func=frim_outname -method=frim_rr_graph -gname=soc-facebook -dir=graphInfo \
  -pdist=wc -rand_seed=1 -eps=0.5 -frim_sweeps=10 -nodehyper_suffix=lam_unirand_0_0.5
```
