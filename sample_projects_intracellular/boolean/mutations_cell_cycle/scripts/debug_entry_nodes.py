#!/usr/bin/env python3
"""
Debugging script to track entry node states (S_entry, G2M_entry, G0G1_entry) over time.
This script analyzes the MultiCellDS output files to track the actual entry node states
and visualize their dynamics along the simulation.
"""

import os
import re
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from collections import defaultdict

# Try to import pcdl, fall back to alternative if not available
try:
    from pcdl import TimeSeries
    PCDL_AVAILABLE = True
except ImportError:
    print("Warning: pcdl not available. Will try alternative methods.")
    PCDL_AVAILABLE = False

# -----------------------
# Matplotlib global settings for publication quality
# -----------------------
mpl.rcParams.update({
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "font.size": 10
})

output_dir = "output"

# Create results directory structure
script_dir = os.path.dirname(os.path.abspath(__file__))
results_dir = os.path.join(script_dir, "results")
debug_dir = os.path.join(results_dir, "entry_nodes_debug")

# Create directories if they don't exist
os.makedirs(debug_dir, exist_ok=True)
print(f"Results will be saved to: {debug_dir}")

def analyze_entry_nodes_from_csv(output_dir):
    """Analyze entry node states from CSV files"""
    
    entry_nodes = ["S_entry", "G2M_entry", "G0G1_entry"]
    
    # Data structure to store results
    node_state_data = {
        'times': [],
        'S_entry_on': [],
        'S_entry_off': [],
        'G2M_entry_on': [],
        'G2M_entry_off': [],
        'G0G1_entry_on': [],
        'G0G1_entry_off': [],
        'total_cells': [],
        'S_entry_percentage': [],
        'G2M_entry_percentage': [],
        'G0G1_entry_percentage': [],
        'cell_states': defaultdict(list)  # cell_id -> [(time, S_entry, G2M_entry, G0G1_entry), ...]
    }
    
    # Find all _mut.csv files
    csv_files = []
    for file in os.listdir(output_dir):
        if file.endswith("_mut.csv"):
            csv_files.append(file)
    
    csv_files.sort(key=lambda x: float(re.search(r'(\d+\.?\d*)', x).group(1)) if re.search(r'(\d+\.?\d*)', x) else 0)
    
    print(f"Found {len(csv_files)} CSV files to process")
    
    for csv_file in csv_files:
        filepath = os.path.join(output_dir, csv_file)
        
        # Extract time from filename
        time_match = re.search(r'(\d+\.?\d*)', csv_file)
        if not time_match:
            continue
        time_step = float(time_match.group(1))
        
        try:
            df = pd.read_csv(filepath)
        except Exception as e:
            print(f"Error reading {csv_file}: {e}")
            continue
        
        if len(df) == 0:
            continue
        
        total_cells = len(df)
        node_state_data['times'].append(time_step)
        node_state_data['total_cells'].append(total_cells)
        
        # Analyze each entry node
        for node_name in entry_nodes:
            if node_name not in df.columns:
                print(f"Warning: {node_name} not found in {csv_file}")
                continue
            
            node_states = df[node_name]
            if node_states.dtype == bool:
                node_on_count = node_states.sum()
            else:
                node_on_count = (node_states >= 0.5).sum()  # Treat values >= 0.5 as ON
            node_off_count = total_cells - node_on_count
            node_on_pct = (node_on_count / total_cells * 100) if total_cells > 0 else 0
            
            # Store counts
            if node_name == "S_entry":
                node_state_data['S_entry_on'].append(node_on_count)
                node_state_data['S_entry_off'].append(node_off_count)
                node_state_data['S_entry_percentage'].append(node_on_pct)
            elif node_name == "G2M_entry":
                node_state_data['G2M_entry_on'].append(node_on_count)
                node_state_data['G2M_entry_off'].append(node_off_count)
                node_state_data['G2M_entry_percentage'].append(node_on_pct)
            elif node_name == "G0G1_entry":
                node_state_data['G0G1_entry_on'].append(node_on_count)
                node_state_data['G0G1_entry_off'].append(node_off_count)
                node_state_data['G0G1_entry_percentage'].append(node_on_pct)
            
            # Store individual cell states
            if 'ID' in df.columns:
                for _, row in df.iterrows():
                    cell_id = int(row['ID'])
                    if isinstance(row[node_name], bool):
                        node_state = row[node_name]
                    else:
                        node_state = bool(row[node_name] >= 0.5)
                    
                    # Find or create entry for this cell at this time
                    found = False
                    for i, (t, s, g2m, g01) in enumerate(node_state_data['cell_states'][cell_id]):
                        if abs(t - time_step) < 0.1:  # Same time point
                            if node_name == "S_entry":
                                node_state_data['cell_states'][cell_id][i] = (time_step, node_state, g2m, g01)
                            elif node_name == "G2M_entry":
                                node_state_data['cell_states'][cell_id][i] = (time_step, s, node_state, g01)
                            elif node_name == "G0G1_entry":
                                node_state_data['cell_states'][cell_id][i] = (time_step, s, g2m, node_state)
                            found = True
                            break
                    
                    if not found:
                        if node_name == "S_entry":
                            node_state_data['cell_states'][cell_id].append((time_step, node_state, False, False))
                        elif node_name == "G2M_entry":
                            node_state_data['cell_states'][cell_id].append((time_step, False, node_state, False))
                        elif node_name == "G0G1_entry":
                            node_state_data['cell_states'][cell_id].append((time_step, False, False, node_state))
    
    # Sort cell states by time
    for cell_id in node_state_data['cell_states']:
        node_state_data['cell_states'][cell_id].sort(key=lambda x: x[0])
    
    return node_state_data

def plot_entry_node_percentages(node_state_data):
    """Plot percentage of cells with each entry node ON over time"""
    
    fig, axes = plt.subplots(3, 1, figsize=(10, 12), sharex=True)
    
    times = np.array(node_state_data['times'])
    
    # S_entry
    axes[0].plot(times, node_state_data['S_entry_percentage'], 'b-', linewidth=2, label='S_entry ON')
    axes[0].fill_between(times, 0, node_state_data['S_entry_percentage'], alpha=0.3, color='blue')
    axes[0].set_ylabel('Percentage of Cells (%)', fontsize=12)
    axes[0].set_title('S_entry Node State Over Time', fontsize=14, fontweight='bold')
    axes[0].grid(True, alpha=0.3)
    axes[0].legend()
    axes[0].set_ylim([0, 105])
    
    # G2M_entry
    axes[1].plot(times, node_state_data['G2M_entry_percentage'], 'g-', linewidth=2, label='G2M_entry ON')
    axes[1].fill_between(times, 0, node_state_data['G2M_entry_percentage'], alpha=0.3, color='green')
    axes[1].set_ylabel('Percentage of Cells (%)', fontsize=12)
    axes[1].set_title('G2M_entry Node State Over Time', fontsize=14, fontweight='bold')
    axes[1].grid(True, alpha=0.3)
    axes[1].legend()
    axes[1].set_ylim([0, 105])
    
    # G0G1_entry
    axes[2].plot(times, node_state_data['G0G1_entry_percentage'], 'r-', linewidth=2, label='G0G1_entry ON')
    axes[2].fill_between(times, 0, node_state_data['G0G1_entry_percentage'], alpha=0.3, color='red')
    axes[2].set_xlabel('Time (min)', fontsize=12)
    axes[2].set_ylabel('Percentage of Cells (%)', fontsize=12)
    axes[2].set_title('G0G1_entry Node State Over Time', fontsize=14, fontweight='bold')
    axes[2].grid(True, alpha=0.3)
    axes[2].legend()
    axes[2].set_ylim([0, 105])
    
    plt.tight_layout()
    
    # Save figure
    output_path = os.path.join(debug_dir, "entry_node_percentages.png")
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_path}")
    
    output_path_pdf = os.path.join(debug_dir, "entry_node_percentages.pdf")
    plt.savefig(output_path_pdf, bbox_inches='tight')
    print(f"Saved: {output_path_pdf}")
    
    plt.close()

def plot_entry_node_counts(node_state_data):
    """Plot counts of cells with each entry node ON/OFF over time"""
    
    fig, axes = plt.subplots(3, 1, figsize=(10, 12), sharex=True)
    
    times = np.array(node_state_data['times'])
    
    # S_entry
    axes[0].plot(times, node_state_data['S_entry_on'], 'b-', linewidth=2, label='S_entry ON')
    axes[0].plot(times, node_state_data['S_entry_off'], 'b--', linewidth=2, label='S_entry OFF')
    axes[0].set_ylabel('Number of Cells', fontsize=12)
    axes[0].set_title('S_entry Node State Counts Over Time', fontsize=14, fontweight='bold')
    axes[0].grid(True, alpha=0.3)
    axes[0].legend()
    
    # G2M_entry
    axes[1].plot(times, node_state_data['G2M_entry_on'], 'g-', linewidth=2, label='G2M_entry ON')
    axes[1].plot(times, node_state_data['G2M_entry_off'], 'g--', linewidth=2, label='G2M_entry OFF')
    axes[1].set_ylabel('Number of Cells', fontsize=12)
    axes[1].set_title('G2M_entry Node State Counts Over Time', fontsize=14, fontweight='bold')
    axes[1].grid(True, alpha=0.3)
    axes[1].legend()
    
    # G0G1_entry
    axes[2].plot(times, node_state_data['G0G1_entry_on'], 'r-', linewidth=2, label='G0G1_entry ON')
    axes[2].plot(times, node_state_data['G0G1_entry_off'], 'r--', linewidth=2, label='G0G1_entry OFF')
    axes[2].set_xlabel('Time (min)', fontsize=12)
    axes[2].set_ylabel('Number of Cells', fontsize=12)
    axes[2].set_title('G0G1_entry Node State Counts Over Time', fontsize=14, fontweight='bold')
    axes[2].grid(True, alpha=0.3)
    axes[2].legend()
    
    plt.tight_layout()
    
    # Save figure
    output_path = os.path.join(debug_dir, "entry_node_counts.png")
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_path}")
    
    output_path_pdf = os.path.join(debug_dir, "entry_node_counts.pdf")
    plt.savefig(output_path_pdf, bbox_inches='tight')
    print(f"Saved: {output_path_pdf}")
    
    plt.close()

def plot_all_entry_nodes_combined(node_state_data):
    """Plot all entry nodes on the same plot for comparison"""
    
    fig, ax = plt.subplots(1, 1, figsize=(12, 6))
    
    times = np.array(node_state_data['times'])
    
    ax.plot(times, node_state_data['S_entry_percentage'], 'b-', linewidth=2, label='S_entry', marker='o', markersize=3)
    ax.plot(times, node_state_data['G2M_entry_percentage'], 'g-', linewidth=2, label='G2M_entry', marker='s', markersize=3)
    ax.plot(times, node_state_data['G0G1_entry_percentage'], 'r-', linewidth=2, label='G0G1_entry', marker='^', markersize=3)
    
    ax.set_xlabel('Time (min)', fontsize=12)
    ax.set_ylabel('Percentage of Cells with Node ON (%)', fontsize=12)
    ax.set_title('Entry Node States Over Time', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=10)
    ax.set_ylim([0, 105])
    
    plt.tight_layout()
    
    # Save figure
    output_path = os.path.join(debug_dir, "all_entry_nodes_combined.png")
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_path}")
    
    output_path_pdf = os.path.join(debug_dir, "all_entry_nodes_combined.pdf")
    plt.savefig(output_path_pdf, bbox_inches='tight')
    print(f"Saved: {output_path_pdf}")
    
    plt.close()

def main():
    """Main function"""
    
    print("=" * 60)
    print("Entry Node State Analysis")
    print("=" * 60)
    
    # Change to script directory to ensure relative paths work
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(os.path.dirname(script_dir))
    output_path = os.path.join(project_root, output_dir)
    
    if not os.path.exists(output_path):
        print(f"Error: Output directory not found: {output_path}")
        return
    
    print(f"Analyzing entry node states from: {output_path}")
    
    # Analyze entry node states
    print("\nAnalyzing entry node states from CSV files...")
    node_state_data = analyze_entry_nodes_from_csv(output_path)
    
    if not node_state_data['times']:
        print("No data found. Make sure the simulation has been run and output files exist.")
        return
    
    print(f"Processed {len(node_state_data['times'])} time points")
    print(f"Tracked {len(node_state_data['cell_states'])} cells")
    
    # Generate plots
    print("\nGenerating plots...")
    plot_entry_node_percentages(node_state_data)
    plot_entry_node_counts(node_state_data)
    plot_all_entry_nodes_combined(node_state_data)
    
    print("\n" + "=" * 60)
    print("Analysis complete!")
    print(f"All plots saved to: {debug_dir}")
    print("=" * 60)

if __name__ == "__main__":
    main()

