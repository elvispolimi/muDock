import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os

def main():
    # 1. Setup command line arguments
    parser = argparse.ArgumentParser(description="Generate heatmaps for LGA vs GA docking scores.")
    parser.add_argument("csv_file", type=str, help="Path to the results CSV file.")
    parser.add_argument("--agg", choices=["mean", "median"], default="mean", 
                        help="Aggregation method to combine seeds (default: mean).")
    parser.add_argument("--outdir", type=str, default=".", 
                        help="Directory to save the generated heatmaps (default: current directory).")
    args = parser.parse_args()

    # 2. Read the CSV file
    # pandas automatically converts "N/A" strings to NaN (Not a Number) by default
    df = pd.read_csv(args.csv_file)

    # Ensure rate and iteration are numeric (NaN for standard GA runs)
    df['rate'] = pd.to_numeric(df['rate'], errors='coerce')
    df['iteration'] = pd.to_numeric(df['iteration'], errors='coerce')

    # 3. Separate GA runs (baseline) and LGA runs
    ga_runs = df[df['rate'].isna() & df['iteration'].isna()]
    lga_runs = df[df['rate'].notna() & df['iteration'].notna()]

    if ga_runs.empty or lga_runs.empty:
        print("Error: Could not find both GA (N/A rate/iteration) and LGA runs in the dataset.")
        return

    # 4. Aggregate the data by ligand (and config) using the chosen metric
    agg_func = args.agg

    # Baseline GA scores per ligand
    ga_baseline = ga_runs.groupby('name')['score'].agg(agg_func).reset_index()
    ga_baseline.rename(columns={'score': 'ga_score'}, inplace=True)

    # LGA scores per ligand, rate, and iteration
    lga_agg = lga_runs.groupby(['name', 'rate', 'iteration'])['score'].agg(agg_func).reset_index()

    # 5. Merge and calculate the difference
    merged_data = pd.merge(lga_agg, ga_baseline, on='name')
    # Difference: LGA - GA (negative means LGA scored lower/better than GA)
    merged_data['score_diff'] = merged_data['score'] - merged_data['ga_score']

    # Create output directory if it doesn't exist
    os.makedirs(args.outdir, exist_ok=True)

    # 6. Generate a heatmap for each ligand
    ligands = merged_data['name'].unique()
    
    for ligand in ligands:
        # Filter data for the specific ligand
        ligand_data = merged_data[merged_data['name'] == ligand]
        
        # Pivot the data to create a 2D matrix for the heatmap
        # Columns = rate (X-axis), Index = iteration (Y-axis)
        pivot_table = ligand_data.pivot(index='iteration', columns='rate', values='score_diff')
        
        # Sort the Y-axis so lower iterations are at the bottom (or reverse if preferred)
        pivot_table = pivot_table.sort_index(ascending=False)

        # Plotting
        plt.figure(figsize=(8, 6))
        sns.heatmap(pivot_table, annot=True, cmap="coolwarm", center=0, fmt=".2f", 
                    cbar_kws={'label': f'Score Difference ({args.agg.capitalize()})'})
        
        plt.title(f"{ligand} - LGA vs GA\n(Score Difference: LGA - GA)")
        plt.xlabel("Local Search Rate")
        plt.ylabel("Local Search Iterations")
        plt.tight_layout()
        
        # Save the figure
        out_path = os.path.join(args.outdir, f"heatmap_{ligand}_{args.agg}.png")
        plt.savefig(out_path, dpi=300)
        plt.close()
        
        print(f"Saved heatmap for {ligand} to {out_path}")

if __name__ == "__main__":
    main()