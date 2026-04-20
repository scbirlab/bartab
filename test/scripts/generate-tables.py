"""Generate a count table with known fitness ratios for bartab validation."""
from typing import Dict, Mapping, Tuple, Optional

from carabiner import print_err
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.stats import multivariate_hypergeom, nbinom

SEED: int = 42
RNG = np.random.default_rng(SEED)

# --- known ground truth ---
FITNESS = {
    "wt": 1.,
    "spike": 0.,   # non-growing spike-in
    "mut_A": 1.5,
    "mut_B": .5,
    "mut_C": .9,
}

# --- experimental parameters ---
N_TIMEPOINTS: int = 6
MAX_TIME: float = 5.
INOCULUM: int = 1000  # cells per strain at t=0
INOCULUM_VAR: float = .1  # between-strain inoculum variance
CARRYING_CAP: float = 50.      # multiplier of inoculum
SAMPLE_FRAC: float = .1
SEQ_DEPTH: int = 100  # reads per strain per sample
SEQ_VAR: float = .01  # sequencing noise
N_REPLICATES: int = 3


def lotka_volterra(t: float, y: int, w: float, K: int):
    remaining_capacity = np.sum(y) / K
    dy = w * y.flatten() * (1. - remaining_capacity)
    return dy


def grow(t: float, n0: int, w: float, K: Optional[int] = None):
    """Deterministic population size at time t.

    Parameters
    ==========
    n0 : initial cells
    w  : growth rate
    t  : time (arbitrary units)
    K  : carrying capacity  (None → pure exponential)
    """
    if K is None:
        return n0 * np.exp(w[None] * t[:,None])
    # logistic with shared K
    else:
        ode_solution = solve_ivp(
            lotka_volterra,
            t_span=sorted(set([0, max(t)])),
            t_eval=sorted(t),
            y0=n0,
            vectorized=True,
            args=(w, K),
        )
        # print(ode_solution)
        return ode_solution.y


def _nb_np_from_mean_var(mean: int, var: float):
    var = mean + var * np.square(mean)
    p = mean / var
    n = np.square(mean) / (var - mean)
    return n, p


def calculate_growth_curves(
    inoculum: int, 
    fitness: dict, 
    inoculum_var: float = .1, 
    carrying_capacity: Optional[int] = None, 
    n_timepoints: int = 100,
    ref_key: str = "wt"
) -> Tuple[np.ndarray, ...]:
    n, p = _nb_np_from_mean_var(inoculum, inoculum_var)
    w = [v for _, v in fitness.items()]
    n0 = nbinom.rvs(n, p, size=len(w), random_state=SEED)
    t = np.linspace(0., MAX_TIME, num=int(n_timepoints))
    growths = grow(t, n0, w, carrying_capacity if carrying_capacity is None else inoculum * carrying_capacity)
    growths = {key: g for (key, _), g in zip(fitness.items(), growths)}
    n0 = {key: _n0 for (key, _), _n0 in zip(fitness.items(), n0)}
    ref_expansion = growths[ref_key] / n0[ref_key]
    return t, ref_expansion, growths


def reads_sampler(
    population: Mapping[str, np.ndarray], 
    seq_depth: int = 1_000, 
    sample_frac: float = 1., 
    reps: int = 1, 
    variance: float = .1
) -> np.ndarray:
    samples = []
    population = np.stack([p for _, p in population.items()], axis=0)  # strains x time
    for i, timepoint_pop in enumerate(population.T):
        sample_size = np.floor(timepoint_pop.sum() * sample_frac).astype(int)
        samples.append(
            multivariate_hypergeom.rvs(
                m=timepoint_pop.flatten().astype(int), 
                n=sample_size, 
                size=reps, 
                random_state=SEED + i,
            ).T
        )
    samples = np.stack(samples, axis=-1)  # strains x reps x time
    n_strains = population.shape[0]
    read_means = np.floor(seq_depth * n_strains * samples / samples.sum(axis=0, keepdims=True))
    n, p = _nb_np_from_mean_var(read_means, variance)
    return np.stack([
        nbinom.rvs(n[:,i], p[:,i], random_state=SEED + i) 
        for i in range(reps)
    ], axis=-2)  # strains x reps x time


if __name__ == "__main__":
    strains = list(FITNESS.keys())

    t, ref_expansion, growths = calculate_growth_curves(
        inoculum=INOCULUM, 
        inoculum_var=.1, 
        carrying_capacity=CARRYING_CAP, 
        fitness=FITNESS, 
        n_timepoints=N_TIMEPOINTS,
    )
    # print_err(growths)
    
    counts = reads_sampler(
        growths, 
        sample_frac=SAMPLE_FRAC, 
        seq_depth=SEQ_DEPTH, 
        reps=N_REPLICATES, 
        variance=SEQ_VAR,
    )
    # print_err(counts)  # strains x reps x time 

    # --- count table: strains × samples ---
    count_df = []
    meta_df = []
    sample_ids = []
    for strain_name, strain_arr in zip(FITNESS, counts):
        for r, rep_arr in enumerate(strain_arr):
            for _t, time_arr in zip(t, rep_arr):
                sample_id = f"t_{_t}-rep_{r}"
                obs_id = f"{sample_id}-strain_{strain_name}"
                count_df.append({
                    "sample_id": sample_id,
                    "obs_id": obs_id,
                    "strain_id": strain_name,
                    "count": time_arr,
                })
                meta_df.append({
                    "sample_id": sample_id,
                    "timepoint": _t,
                    "replicate": r,
                    "is_t0": _t == 0.,
                })
    
    counts_df = pd.DataFrame(count_df)
    meta_df = pd.DataFrame(meta_df).drop_duplicates()

    strain_meta_df = pd.DataFrame({
        "strain_id": list(FITNESS),
        "true_fitness": [FITNESS[s] for s in FITNESS],
        "is_reference": [s == "wt" for s in FITNESS],
        "is_spike": [s == "spike" for s in FITNESS],
    })

    counts_df.to_csv("test/inputs/test_counts.csv", index=False)
    meta_df.to_csv("test/inputs/test_sample_meta.csv", index=False)
    strain_meta_df.to_csv("test/inputs/test_strain_meta.csv", index=False)

    print_err(counts_df)
    print_err("\nSample metadata:")
    print_err(meta_df)
    print_err("\nStrain metadata (ground truth):")
    print_err(strain_meta_df)
