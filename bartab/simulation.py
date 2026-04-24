"""Tools to simulate count data from strain pools."""
from typing import Optional, Mapping, Tuple

from numpy.typing import ArrayLike
import numpy as np
from scipy.integrate import solve_ivp


def lotka_volterra(
    t: float, 
    y: int, 
    w: float, 
    K: int
):
    remaining_capacity = np.sum(y) / K
    dy = w * y.flatten() * (1. - remaining_capacity)
    return dy


def grow(
    t: float, 
    n0: int, 
    w: float, 
    K: Optional[int] = None
):
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
            # vectorized=True,
            args=(w, K),
        )
        if ode_solution.success:
            return ode_solution.y
        else:
            raise ValueError(f"ODE failed: {t=}, {K=}, {w=}, {n0=}")

def _nb_np_from_mean_var(mean: int, var: float):
    var = mean + var * np.square(mean)
    p = mean / var
    n = np.square(mean) / (var - mean)
    return n, p


def calculate_growth_curves(
    inoculum: int, 
    fitness: Mapping[str, ArrayLike], 
    inoculum_var: float = .1, 
    carrying_capacity: Optional[int] = None, 
    n_timepoints: int = 5,
    max_time: float = 10.,
    ref_key: str = "wt",
    seed: int = 42
) -> Tuple[np.ndarray, ...]:
    """
    Examples
    ========
    >>> from bartab.simulation import calculate_growth_curves
    >>> t, expansion, growths = calculate_growth_curves(
    ...     inoculum=1_000,
    ...     fitness={"wt": 1., "mut": .5},
    ...     n_timepoints=3,
    ...     max_time=5.,
    ...     seed=42,
    ... )
    >>> t.shape
    (3,)
    >>> sorted(growths.keys())
    ['mut', 'wt']
    >>> growths["wt"].shape  # one value per timepoint  # doctest: +SKIP
    (3,)

    """
    from scipy.stats import nbinom

    if carrying_capacity is not None and isinstance(carrying_capacity, (int, float)):
        carrying_capacity *= inoculum * len(fitness)

    n, p = _nb_np_from_mean_var(inoculum, inoculum_var)
    w = np.asarray([v for _, v in fitness.items()])
    t = np.linspace(0., max_time, num=int(n_timepoints))
    # print("t", t)
    # n0 = nbinom.rvs(n, p, size=len(w), random_state=seed)
    n0 = _safe_nbinom_rvs(n, p, size=w.size, random_state=seed)
    n0_named = {key: _n0 for (key, _), _n0 in zip(fitness.items(), n0)}
    
    growths = grow(
        t, n0, w, 
        carrying_capacity,
    )
    # print(f"{growths.shape=}", t.shape, n0.shape, w.shape)
    growths = {key: g for (key, _), g in zip(fitness.items(), growths)}
    
    ref_expansion = growths[ref_key] / n0_named[ref_key]

    return t, ref_expansion, growths


def _safe_nbinom_rvs(
    n,
    p,
    random_state: int = 42,
    **kwargs
):
    from scipy.stats import nbinom
    if np.all(n > 0 & np.isfinite(n)):
        return nbinom.rvs(n, p, random_state=random_state,  **kwargs)
    else:
        out = np.empty_like(n, dtype=np.int64).ravel()
        n_flat = n.ravel()
        p_flat = p.ravel()
        for i, (_n, _p) in enumerate(zip(n_flat, p_flat)):
            if _n > 0 and np.isfinite(_n):
                out[i] = nbinom.rvs(
                    _n, _p,
                    random_state=random_state + i,  
                    **kwargs,
                )
            else:
                out[i] = 0
        return np.reshape(out, n.shape)


def reads_sampler(
    population: Mapping[str, ArrayLike], 
    seq_depth: int = 1_000, 
    sample_frac: float = 1., 
    reps: int = 1, 
    variance: float = .1,
    seed: int = 42
) -> np.ndarray:
    """Simulate sampling sequencing reaad counts from a pooled growth curve.

    Examples
    ========
    >>> import numpy as np; np.set_printoptions(legacy='1.25')
    >>> from bartab.simulation import reads_sampler
    >>> population = {
    ...     "wt":  np.array([1_000, 2_718, 7_389]),
    ...     "mut": np.array([  500,   800, 1_200]),
    ... }
    >>> counts = reads_sampler(population, seq_depth=100, reps=2, seed=42)
    >>> counts.shape  # (n_strains, reps, n_timepoints)
    (2, 2, 3)
    >>> (counts >= 0).all()
    True

    """

    from scipy.stats import multivariate_hypergeom, nbinom

    samples = []
    # print("?", population[list(population)[0]].shape, len(population))
    population = np.stack([p for _, p in population.items()], axis=0)  # strains x time
    for i, timepoint_pop in enumerate(population.T):
        sample_size = np.floor(timepoint_pop.sum() * sample_frac).astype(int)
        samples.append(
            multivariate_hypergeom.rvs(
                m=timepoint_pop.flatten().astype(int), 
                n=sample_size, 
                size=reps, 
                random_state=seed + i,
            ).T
        )
    samples = np.stack(samples, axis=-1)  # strains x reps x time
    n_strains = population.shape[0]
    read_means = np.floor(
        seq_depth * n_strains * samples 
        / samples.sum(axis=0, keepdims=True)
    )
    n, p = _nb_np_from_mean_var(read_means, variance)
    # print("???", population.shape, samples.shape, read_means.shape)
    samples = np.stack([
        _safe_nbinom_rvs(n[:,i], p[:,i], random_state=seed + i)
        for i in range(reps)
    ], axis=-2)  # strains x reps x time
    # neg_samples = samples < 0
    # print("???", n[neg_samples], p[neg_samples])
    return samples
