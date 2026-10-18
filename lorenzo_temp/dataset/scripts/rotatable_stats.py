
#!/usr/bin/env python3

import argparse
import numpy as np


def main():
    parser = argparse.ArgumentParser(
        description="Compute statistics for rotatable-bond counts."
    )
    parser.add_argument(
        "bonds",
        nargs="+",
        type=int,
        help="List of rotatable-bond counts, e.g. 3 5 8 12 15"
    )
    parser.add_argument(
        "--percentiles",
        nargs="+",
        type=float,
        default=[25, 50, 75, 90, 95],
        help="Percentiles to compute (default: 25 50 75 90 95)"
    )
    args = parser.parse_args()

    bonds = np.array(args.bonds, dtype=float)

    if np.any(bonds < 0):
        parser.error("Rotatable-bond counts cannot be negative.")

    if np.any((np.array(args.percentiles) < 0) |
              (np.array(args.percentiles) > 100)):
        parser.error("Percentiles must be between 0 and 100.")

    print(f"Number of ligands: {len(bonds)}")
    print(f"Mean:               {np.mean(bonds):.2f}")
    print(f"Median:             {np.median(bonds):.2f}")
    print(f"Std dev (pop.):     {np.std(bonds, ddof=0):.2f}")
    print(f"Std dev (sample):   {np.std(bonds, ddof=1):.2f}"
          if len(bonds) > 1 else "Std dev (sample):   N/A")

    print("\nPercentiles:")
    for p in args.percentiles:
        value = np.percentile(bonds, p)
        print(f"  P{p:g}: {value:.2f}")


if __name__ == "__main__":
    main()