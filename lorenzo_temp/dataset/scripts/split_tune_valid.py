#!/usr/bin/env python3
"""Split a table into TUNE (20 rows, bins 7/8/5 by `rot`) and VALID (the rest),
keeping every `target` value entirely inside one subset.

Usage: python split_tune_valid.py table.csv [seed]
"""
import random
import sys

import pandas as pd

BIN_SIZES = (7, 8, 5)  # bin1: rot 1-9, bin2: rot 10-16, bin3: rot >= 17


def rot_bin(rot):
    if 1 <= rot <= 9:
        return 0
    if 10 <= rot <= 16:
        return 1
    if rot >= 17:
        return 2
    return None  # rot < 1: cannot go into TUNE


def split(df, seed=0):
    rng = random.Random(seed)
    df = df.copy()
    df["_bin"] = df["rot"].apply(rot_bin)

    # One item per target: the vector of how many rows it would add to each bin.
    # A target is eligible for TUNE only if ALL its rows fit in a bin and the
    # vector does not exceed the bin sizes.
    items = []
    for target, g in df.groupby("target"):
        if g["_bin"].isna().any():
            continue
        vec = [int((g["_bin"] == b).sum()) for b in range(3)]
        if all(v <= s for v, s in zip(vec, BIN_SIZES)):
            items.append((target, tuple(vec)))
    rng.shuffle(items)  # different seeds -> different valid splits

    # DP over reachable states (n1, n2, n3); parent pointers to rebuild the choice.
    start = (0, 0, 0)
    parent = {start: None}  # state -> (previous_state, target)
    for target, vec in items:
        for state in list(parent):  # snapshot: each target used at most once
            new = tuple(s + v for s, v in zip(state, vec))
            if new not in parent and all(n <= s for n, s in zip(new, BIN_SIZES)):
                parent[new] = (state, target)

    goal = BIN_SIZES
    if goal not in parent:
        raise ValueError("No split satisfies the constraints (try relaxing bin sizes).")

    chosen, state = set(), goal
    while parent[state] is not None:
        state, target = parent[state]
        chosen.add(target)

    tune = df[df["target"].isin(chosen)].drop(columns="_bin")
    valid = df[~df["target"].isin(chosen)].drop(columns="_bin")

    # sanity checks
    assert len(tune) == sum(BIN_SIZES)
    assert set(tune["target"]).isdisjoint(valid["target"])
    counts = tune["rot"].apply(rot_bin).value_counts().reindex(range(3), fill_value=0)
    assert tuple(counts) == BIN_SIZES
    return tune, valid


if __name__ == "__main__":
    path = sys.argv[1]
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    table = pd.read_csv(path)
    tune, valid = split(table, seed)
    tune.to_csv("TUNE.csv", index=False)
    valid.to_csv("VALID.csv", index=False)
    print(f"TUNE: {len(tune)} rows, VALID: {len(valid)} rows")