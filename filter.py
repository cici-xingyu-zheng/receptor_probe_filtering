import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt

def less_than_10_check(probe_df):
    '''
    Function to quickly check how many genes have less than 10 probes designed for them
    '''
    counts = probe_df['gene_name'].value_counts()
    counts_df = counts.rename_axis('gene_name').reset_index(name='count')
    less_than_10_probes_counts_df = counts_df[counts_df['count'] < 10].reset_index(drop=True)

    print(f'Number of genes with less than 10 probes: {len(less_than_10_probes_counts_df)}')

    plt.figure(figsize=(6, 4))
    less_than_10_probes_counts_df['count'].plot(kind='hist')
    plt.title('Less than 10 probes per gene count distribution')
    plt.show()

    return less_than_10_probes_counts_df

def drop_global_offenders(
    P, names, remove_frac=None, sum_threshold=None
):
    """
    Greedy removal for heterodimer fraction matrix P 
    At each step, remove the probe with the largest per-probe sum s_i = sum_j P[i, j],
    then update sums by subtracting column j.

    Stop when:
      - you've removed floor(remove_frac * N0) items (if remove_frac is given), OR
      - the current max sum <= sum_threshold (if sum_threshold is given).
    At least one of (remove_frac, sum_threshold) must be provided.

    Returns:
      keep_mask  : bool array shape (N,)  True = kept
      removed_log: DataFrame with per-iteration info
    """

    N0 = P.shape[0]
    if remove_frac is None and sum_threshold is None:
        raise ValueError("Provide remove_frac and/or sum_threshold")


    # initial sums and active set
    cur_sum = P.sum(axis=1).copy()
    active = np.ones(N0, dtype=bool)

    # compute removal target if remove_frac provided
    target_remove = None
    if remove_frac is not None:
        target_remove = int(remove_frac * N0)

    removed = []
    removed_count = 0

    while True:
        # 1) stop by count
        if target_remove is not None and removed_count >= target_remove:
            break

        # pick current worst among active
        cand = np.where(active)[0]
        if cand.size == 0:
            break
        j = cand[np.argmax(cur_sum[cand])]
        smax = cur_sum[j]

        # 2) stop by threshold
        if sum_threshold is not None and smax <= sum_threshold:
            break

        # log and remove j
        removed.append({
            "iter": removed_count,
            "idx": int(j),
            "name": names[j],
            "sum_before": float(smax),
            "mean_active": float(cur_sum[cand].mean()),
            "std_active": float(cur_sum[cand].std(ddof=0)),
            "n_active": int(cand.size),
        })
        active[j] = False
        removed_count += 1

        # subtract column j from remaining sums
        cur_sum[active] -= P[active, j]
        cur_sum[j] = -np.inf  # set to neg inf so it won't be picked again

    removed_log = pd.DataFrame(removed)
    keep_mask = active

    return keep_mask, removed_log

def collapse_duplicated(
        df, subject_col="subject", verbose=True):
    
    '''
    Function to remove duplicate BLAST hits due to the fact we include both the annotated OR genes and
    the Ensembl background as subjects in the BLAST database.
    Collapse rows in df that are identical in all columns except `subject_col`,
    also ignoring any `sstart`/`send` columns if present
    Returns: 
        Collapsed DataFrame removing duplicates 
    '''
    # drop start/end 
    to_drop = ["sstart", "send"] 
    work = df.drop(columns=to_drop, errors="ignore").copy()
    n_before = len(work)

    # Define identical groups by all not subject columns
    key_cols = [c for c in work.columns if c != subject_col]
    g = work.groupby(key_cols, dropna=False, sort=False)

    # find ambiguous groups (where base is not unique)
    report = g.apply(
        lambda x: (
            x[subject_col].astype("string").str.strip()
              .str.split("|", n=1, expand=True).iloc[:, 0]
              .unique().tolist()
        )
    ).rename("bases").reset_index()

    ambiguous_keys = report[report["bases"].apply(lambda xs: len(xs) > 1)][key_cols]
    unambiguous = report[report["bases"].apply(lambda xs: len(xs) == 1)].copy()
    unambiguous["base"] = unambiguous["bases"].str[0]

    # merge base back, collapse unambiguous groups to a single row with subject=base
    merged = work.merge(unambiguous[key_cols + ["base"]], how="left", on=key_cols, copy=False)

    def _collapse(gr: pd.DataFrame) -> pd.DataFrame:
        if pd.isna(gr["base"].iloc[0]):
            # ambiguous or singleton not listed -> keep as-is
            return gr
        out = gr.iloc[[0]].copy()
        out[subject_col] = gr["base"].iloc[0]
        return out

    collapsed = (
        merged.groupby(key_cols, dropna=False, sort=False, group_keys=False)
              .apply(_collapse)
              .drop(columns=["base"], errors="ignore")
              .reset_index(drop=True)
    )

    if verbose:
        n_after = len(collapsed)
        n_ambiguous = len(ambiguous_keys)

        if to_drop:
            print(f"Dropped columns: {to_drop}")
        print(f"Rows before: {n_before}")
        print(f"Rows after:  {n_after}")

        if n_ambiguous:
            print("\n[Note] Example ambiguous key (bases differ within same key):")
            print(ambiguous_keys.head(1).to_string(index=False))

    return collapsed

def keep_top_per_gene(
    df,
    *,
    score_col="score_total",
    gene_col="gene",
    name_col="padlock_name",   # set to "name" if that's your column
    max_per_gene=10,
    add_overall_ranks=True,
    all_names=None,            # optional: full original name list for full mask
    prior_keep_mask=None,      # optional: bool mask aligned to all_names
    return_full_mask=False,    # set True to also return keep_mask_full
):
    
    """
    Filter DataFrame to keep only top K entries per gene,
    based on ranking by score_col (lower = better).
    Also adds per-gene counts and ranks, and overall ranks/percentiles.
    Returns:
        out: DataFrame with added columns:
                        - n_padlocks_in_gene
                        - rank_in_gene_best
                        - (optional) rank_overall_best
                        - (optional) overall_percentile
                        - keep (bool)
        (optional) keep_mask_full : bool array aligned to all_names, combining prior_keep_mask if given
    """
    out = df.copy()

    # overall ranks / percentile (QC)
    if add_overall_ranks:
        out["rank_overall_best"] = out[score_col].rank(ascending=True, method="dense").astype(int)
        r = out[score_col].rank(ascending=True, method="average")
        out["overall_percentile"] = (r - 0.5) / len(out)

    # per-gene counts
    out["n_padlocks_in_gene"] = out.groupby(gene_col)[name_col].transform("count")

    # per-gene rank (1 = best) – compute via a stable sort, then merge back
    tmp = out[[gene_col, score_col, name_col]].copy()
    tmp = tmp.sort_values([gene_col, score_col, name_col], kind="mergesort")
    tmp["rank_in_gene_best"] = tmp.groupby(gene_col).cumcount() + 1

    out = out.merge(tmp[[name_col, "rank_in_gene_best"]],
                    on=name_col, how="left", validate="one_to_one")

    if "rank_in_gene_best" not in out.columns:
        raise RuntimeError("Failed to create 'rank_in_gene_best' – check name_col uniqueness.")

    # keep top K within genes (>K then drop rest; <= K then keep all)
    mask_large = out["n_padlocks_in_gene"] > max_per_gene
    out["keep"] = True
    out.loc[mask_large & (out["rank_in_gene_best"] > max_per_gene), "keep"] = False

    if not return_full_mask:
        return out

    # build full-length mask aligned to all_names (original ordering)
    if all_names is None:
        raise ValueError("all_names must be provided when return_full_mask=True")

    all_names = np.asarray(all_names)
    # start from prior mask or all True
    keep_full = np.array(prior_keep_mask, dtype=bool) if prior_keep_mask is not None else np.ones(all_names.size, bool)

    # map new decisions
    name_to_keep = dict(zip(out[name_col], out["keep"].astype(bool)))
    for i, nm in enumerate(all_names):
        if nm in name_to_keep:
            keep_full[i] = keep_full[i] and name_to_keep[nm]
        # else leave as-is (remains False if prior mask had dropped it)

    return out,keep_full


# A more tail targeting scoring function than z-score:
def tail_signed_score(x, direction="upper", method="logit", eps=1e-12):
    """
    direction: "upper" for higher-worse, "lower" for lower-worse
    method: 
        "logit" (strong tail emphasis) or 
        "invnorm" (phi^{-1}, smoother)
    returns s with median around 0; larger s = worse
    """
    x = pd.to_numeric(pd.Series(x), errors="coerce").to_numpy()
    # percentile in (0,1): average ranks -> plotting position
    r = stats.rankdata(x, method="average")
    u = (r - 0.5) / len(x)
    u = np.clip(u, eps, 1 - eps)

    if method == "invnorm":
        base = stats.norm.ppf(u)             # ~ N(0,1)
    elif method == "logit":
        base = np.log(u / (1 - u))           # harsher on tails
    else:
        raise ValueError("method must be 'invnorm' or 'logit'")

    # orient so larger = worse
    if direction == "upper":      # punish high values
        s = base
    elif direction == "lower":    # punish low values
        s = -base
    else:
        raise ValueError("direction must be 'upper' or 'lower'")
    return s
