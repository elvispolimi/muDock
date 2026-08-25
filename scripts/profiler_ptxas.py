import os
import re
import subprocess
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

REPO_DIR = "/work/onedina/muDock_ON"
BUILD_CUDA = os.path.join(REPO_DIR, "build/cuda")
BUILD_ALPAKA = os.path.join(REPO_DIR, "build/alpaka-cuda")

def demangle(name):
    try:
        res = subprocess.run(["c++filt", name], stdout=subprocess.PIPE, text=True)
        return res.stdout.strip()
    except:
        return name

def get_ptxas_output(build_dir):
    print(f"Cleaning {build_dir}...")
    subprocess.run(["make", "clean"], cwd=build_dir, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    
    print(f"Compiling {build_dir} sequentially (-j1) to safely extract ptxas info...")
    # Use -j1 to prevent interleaved ptxas lines from different translation units
    result = subprocess.run(["make", "-j1"], cwd=build_dir, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    return result.stdout

def parse_ptxas(output, backend_name):
    records = []
    current_func = None
    
    for line in output.split('\n'):
        if "ptxas info" not in line:
            continue
            
        if "Compiling entry function" in line:
            m = re.search(r"Compiling entry function '([^']+)'", line)
            if m:
                mangled = m.group(1)
                demangled = demangle(mangled)
                
                # Tag the kernel based on known keywords
                if "adt_score" in demangled or "calc_energy" in demangled:
                    current_func = "adt_score"
                elif "genetic" in demangled or "iterate" in demangled or "finalize" in demangled:
                    current_func = "genetic"
                elif "geom_transform" in demangled or "apply" in demangled:
                    current_func = "geom_transform"
                else:
                    current_func = "other"
                    
        elif "Used" in line and "registers" in line and current_func:
            reg_m = re.search(r"Used (\d+) registers", line)
            cmem_m = re.search(r"(\d+) bytes cmem", line)
            lmem_m = re.search(r"(\d+) bytes lmem", line)
            smem_m = re.search(r"(\d+) bytes smem", line)
            
            if reg_m and current_func != "other":
                records.append({
                    "Backend": backend_name,
                    "Kernel": current_func,
                    "Registers": int(reg_m.group(1)),
                    "Cmem": int(cmem_m.group(1)) if cmem_m else 0,
                    "Lmem": int(lmem_m.group(1)) if lmem_m else 0,
                    "Smem": int(smem_m.group(1)) if smem_m else 0
                })
            current_func = None 
            
    return records

def main():
    print("=== PTXAS Static Profiler ===")
    
    # 1. Compile and parse CUDA
    out_cuda = get_ptxas_output(BUILD_CUDA)
    rec_cuda = parse_ptxas(out_cuda, "CUDA Native")
    
    # 2. Compile and parse Alpaka
    out_alpaka = get_ptxas_output(BUILD_ALPAKA)
    rec_alpaka = parse_ptxas(out_alpaka, "Alpaka CUDA")
    
    # Combine and save to CSV
    df = pd.DataFrame(rec_cuda + rec_alpaka)
    
    # If a kernel has multiple variants (e.g. templates), average or take max.
    # Usually we want the max registers to be safe.
    if not df.empty:
        df = df.groupby(['Backend', 'Kernel']).max().reset_index()
    
        csv_path = "ptxas_metrics.csv"
        df.to_csv(csv_path, index=False)
        print(f"\n[+] Saved metrics to {csv_path}")
        print("\nExtracted Data:")
        print(df.to_string(index=False))
        
        # 3. Plotting
        sns.set_theme(style="whitegrid")
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        
        # Plot Registers
        sns.barplot(data=df, x="Kernel", y="Registers", hue="Backend", ax=axes[0], palette="Set2")
        axes[0].set_title("Register Usage per Kernel (Lower is Better)")
        axes[0].set_ylabel("Registers per Thread")
        
        # Add a red dashed line at 64 registers (common sweet spot for 100% occupancy on A100)
        axes[0].axhline(64, color='red', linestyle='--', alpha=0.5, label='64 Regs (100% Occupancy limit)')
        axes[0].legend()
        for container in axes[0].containers:
            axes[0].bar_label(container, fmt='%.0f', padding=3)
        
        # Plot Constant Memory
        sns.barplot(data=df, x="Kernel", y="Cmem", hue="Backend", ax=axes[1], palette="Set2")
        axes[1].set_title("Constant Memory Usage (Bytes)")
        axes[1].set_ylabel("Bytes")
        for container in axes[1].containers:
            axes[1].bar_label(container, fmt='%.0f', padding=3)
        
        plt.tight_layout()
        pdf_path = "ptxas_analysis.pdf"
        plt.savefig(pdf_path, format="pdf")
        print(f"[+] Saved plot to {pdf_path}")
        
    else:
        print("[-] No PTXAS data extracted. Did CMake enable -Xptxas -v ?")

if __name__ == "__main__":
    main()
