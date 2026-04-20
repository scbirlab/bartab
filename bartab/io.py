"""Input and output utilities."""
from typing import TYPE_CHECKING, Any, Iterable, Optional, Union
from functools import partial

if TYPE_CHECKING:
    from anndata import AnnData
    from pandas import DataFrame
else:
    AnnData, DataFrame = Any, Any

from carabiner import cast, print_err


def _read_file_or_df(df: Union[str, DataFrame]) -> DataFrame:
    from carabiner import cast, print_err
    from carabiner.pd import read_table
    import pandas as pd
    if isinstance(df, str):
        print_err(f"[INFO] Reading from file {df}")
        return read_table(df)
    elif isinstance(df, pd.DataFrame):
        return df


def _check_cols(df: DataFrame, required: Iterable[str]) -> None:
    missing = set(cast(required, to=list)) - set(df.columns)
    if missing:
        print_err(df.head())
        raise ValueError(f"Required columns are missing: {', '.join(missing)}")
    else:
        return None


def _make_index_value(x: Union[str, Iterable[str]], sep: str = "::"):
    x = cast(x, to=list)
    return sep.join(str(v) for v in x)


def _check_control(
    df: DataFrame, 
    col: Union[str, Iterable[str]], 
    ref,
    name: str,
    check_presence: bool = False,
    sep: str = "::"
):
    new_col = f"__is_{name}__"
    if ref is None:
        df[new_col] = False
    else:
        if col == "__index__":
            ref = _make_index_value(ref, sep=sep)
        if col in df.columns:
            vals_to_check = df[col]
        elif col in df.index.names:
            vals_to_check = df.index.get_level_values(col)
        elif col == "__index__":
            vals_to_check = df.index.values
        else:
            raise KeyError(f"Column {col} not in data")
        df[new_col] = vals_to_check == ref
    n_refs = df[new_col].sum()
    if n_refs == 0 and check_presence:
        raise ValueError(f"No reference samples '{ref}' identified in '{col}'")
    else:
        return df


def _make_index0(df: DataFrame, sep: str = "::"):
    return df.apply(partial(_make_index_value, sep=sep), axis=1)


def _make_index(df: DataFrame, col: Union[str, Iterable[str]], sep: str = "::", name="__index__"):
    col = cast(col, to=list)
    idx = _make_index0(df[col], sep=sep)
    return df.assign(**{name: idx})


def load_anndata(
    counts: Union[str, DataFrame],
    sample_meta: Union[str, DataFrame],
    strain_meta: Union[str, DataFrame],
    reference: str = "wt",
    count_column: str = "count",
    timepoint_column: str = "timepoint",
    t0: Optional[Union[str, float, int]] = None,
    strain_id: Union[str, Iterable[str]] = "barcode_id",
    sample_id: Union[str, Iterable[str]] = "sample_id",
    culture_id: Union[str, Iterable[str]] = "culture_id",
    spike: Optional[str] = None
) -> AnnData:

    from anndata import AnnData
    
    counts = _read_file_or_df(counts)
    sample_meta = _read_file_or_df(sample_meta)
    strain_meta = _read_file_or_df(strain_meta)

    strain_id = cast(strain_id, to=list)
    sample_id = cast(sample_id, to=list)
    culture_id = cast(culture_id, to=list)

    _check_cols(counts, strain_id + sample_id + [count_column])
    _check_cols(sample_meta, sample_id + culture_id + [timepoint_column])
    _check_cols(strain_meta, strain_id)

    sample_meta = _make_index(sample_meta, sample_id)
    sample_meta = _make_index(sample_meta, culture_id, name="__culture_index__")
    strain_meta = _make_index(strain_meta, strain_id)

    # validation
    n_neg = (counts[count_column] < 0).sum()
    if n_neg > 0:
        raise ValueError(f"Counts contains {n_neg} negative values")

    # infer is_t0 if absent
    if t0 is None:
        t0 = sample_meta[timepoint_column].min()

    n_timepoints = sample_meta[timepoint_column].nunique()
    if n_timepoints <= 1:
        raise ValueError(f"Need more than one timepoint, found {n_timepoints}")

    # long → wide: strains × samples
    counts = _make_index(counts, strain_id, name="__row_index__")
    counts = _make_index(counts, sample_id, name="__col_index__")
    counts_wide = (
        counts
        .pivot(
            index="__row_index__", 
            columns="__col_index__", 
            values=count_column,
        )
    )
    # align metadata to matrix axes
    sample_meta = sample_meta.set_index("__index__").loc[counts_wide.columns].copy()
    strain_meta = strain_meta.set_index("__index__").loc[counts_wide.index].copy()

    sample_meta = _check_control(sample_meta, timepoint_column, t0, "t0", check_presence=True)
    strain_meta = _check_control(strain_meta, "__index__", reference, "reference", check_presence=True)
    strain_meta = _check_control(strain_meta, "__index__", spike, "spike")

    adata = AnnData(
        X=counts_wide.values.astype(float),
        obs=strain_meta,
        var=sample_meta,
    )
    adata.uns["reference"] = reference
    adata.uns["spike"] = spike

    return adata
