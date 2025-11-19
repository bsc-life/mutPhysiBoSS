#!/usr/bin/env python3
"""
Debugging script to track PLK1 rate mutations over time.
This script analyzes the .mut.csv files to verify that PLK1 rate mutations
are occurring and tracks their distribution throughout the simulation.
"""

import os
import re
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from collections import defaultdict

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
debug_dir = os.path.join(results_dir, "plk1_debug")

# Create directories if they don't exist
os.makedirs(debug_dir, exist_ok=True)
print(f"Results will be saved to: {debug_dir}")

# Extract time index from filename
def extract_time(filename):
    match = re.search(r'output(\d+)_mut', filename)
    return int(match.group(1)) if match else -1

def parse_mutations(mutation_string):
    """Parse mutation string into individual mutations"""
    if pd.isna(mutation_string) or mutation_string == '':
        return []
    # Mutations are separated by periods
    mutations = mutation_string.split('.')
    # Filter out empty strings
    return [m for m in mutations if m.strip()]

def extract_plk1_mutations(mutations_list):
    """Extract PLK1-related mutations from a list of mutations"""
    plk1_mutations = []
    for mut in mutations_list:
        if 'plk1' in mut.lower() and 'rate' in mut.lower():
            plk1_mutations.append(mut)
    return plk1_mutations

def parse_plk1_mutation(mutation_string):
    """Parse a PLK1 mutation string to extract type and effect size
    
    Format: Plk1_rate_{up|down}_{effect_size}
    Example: Plk1_rate_down_1.000000
    """
    pattern = r'Plk1_rate_(up|down)_([\d.]+)'
    match = re.match(pattern, mutation_string, re.IGNORECASE)
    if match:
        mutation_type = match.group(1).lower()
        effect_size = float(match.group(2))
        return {
            'type': mutation_type,
            'effect_size': effect_size,
            'raw': mutation_string
        }
    return None

def analyze_plk1_mutations(output_dir):
    """Analyze PLK1 rate mutations from .mut.csv files"""
    
    # Read all .mut.csv files
    mut_files = [f for f in os.listdir(output_dir) if f.endswith('mut.csv')]
    mut_files = sorted(mut_files, key=extract_time)
    
    print(f"Found {len(mut_files)} mutation files")
    
    # Initialize data structures
    time_data = []
    plk1_stats = {
        'cells_with_plk1_mutations': [],
        'total_plk1_mutations': [],
        'plk1_up_mutations': [],
        'plk1_down_mutations': [],
        'effect_sizes': defaultdict(list),
        'first_mutation_time': None,
        'mutation_details': []
    }
    
    for mut_file in mut_files:
        time_idx = extract_time(mut_file)
        if time_idx == -1:  # Skip non-time-indexed files
            continue
            
        try:
            df = pd.read_csv(os.path.join(output_dir, mut_file))
            
            # Count cells with PLK1 mutations
            cells_with_plk1 = 0
            total_plk1_mut_count = 0
            up_count = 0
            down_count = 0
            effect_sizes_in_timestep = []
            
            for idx, row in df.iterrows():
                mutations_list = parse_mutations(row['mutations'])
                plk1_mutations = extract_plk1_mutations(mutations_list)
                
                if plk1_mutations:
                    cells_with_plk1 += 1
                    total_plk1_mut_count += len(plk1_mutations)
                    
                    # Parse each PLK1 mutation
                    for plk1_mut in plk1_mutations:
                        parsed = parse_plk1_mutation(plk1_mut)
                        if parsed:
                            if parsed['type'] == 'up':
                                up_count += 1
                            else:
                                down_count += 1
                            
                            effect_sizes_in_timestep.append(parsed['effect_size'])
                            
                            # Record first mutation time
                            if plk1_stats['first_mutation_time'] is None:
                                plk1_stats['first_mutation_time'] = time_idx
                            
                            # Store detailed mutation info
                            plk1_stats['mutation_details'].append({
                                'time': time_idx,
                                'cell_id': row['ID'],
                                'generation': row['generation'],
                                'type': parsed['type'],
                                'effect_size': parsed['effect_size'],
                                'mutation_string': plk1_mut
                            })
            
            plk1_stats['cells_with_plk1_mutations'].append(cells_with_plk1)
            plk1_stats['total_plk1_mutations'].append(total_plk1_mut_count)
            plk1_stats['plk1_up_mutations'].append(up_count)
            plk1_stats['plk1_down_mutations'].append(down_count)
            plk1_stats['effect_sizes'][time_idx] = effect_sizes_in_timestep
            
            time_data.append(time_idx)
            
            if time_idx % 10 == 0 or cells_with_plk1 > 0:
                print(f"Time {time_idx}: {cells_with_plk1} cells with PLK1 mutations, "
                      f"{total_plk1_mut_count} total PLK1 mutations "
                      f"(up: {up_count}, down: {down_count})")
            
        except Exception as e:
            print(f"Error processing {mut_file}: {e}")
            continue
    
    return plk1_stats, time_data

def plot_plk1_mutation_timeline(plk1_stats, time_points, output_path):
    """Plot timeline of PLK1 mutations"""
    
    if not time_points:
        print("No time points available for plotting")
        return
    
    fig, axes = plt.subplots(2, 1, figsize=(14, 10))
    
    # Plot 1: Cells with PLK1 mutations over time
    ax1 = axes[0]
    ax1.plot(time_points, plk1_stats['cells_with_plk1_mutations'], 
             'b-', linewidth=2, marker='o', markersize=4, label='Cells with PLK1 mutations')
    ax1.set_xlabel('Time Step', fontsize=12)
    ax1.set_ylabel('Number of Cells', fontsize=12)
    ax1.set_title('Cells with PLK1 Rate Mutations Over Time', fontsize=14)
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    
    # Plot 2: Total PLK1 mutations over time (stacked)
    ax2 = axes[1]
    ax2.plot(time_points, plk1_stats['plk1_up_mutations'], 
             'g-', linewidth=2, marker='s', markersize=4, label='PLK1 rate UP mutations')
    ax2.plot(time_points, plk1_stats['plk1_down_mutations'], 
             'r-', linewidth=2, marker='^', markersize=4, label='PLK1 rate DOWN mutations')
    ax2.set_xlabel('Time Step', fontsize=12)
    ax2.set_ylabel('Number of Mutations', fontsize=12)
    ax2.set_title('PLK1 Rate Mutation Counts Over Time', fontsize=14)
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    
    plt.tight_layout()
    plt.savefig(f"{output_path}_timeline.png", bbox_inches='tight', dpi=300)
    plt.savefig(f"{output_path}_timeline.pdf", bbox_inches='tight')
    plt.close()
    
    print(f"Timeline plot saved to: {output_path}_timeline.png/pdf")

def get_base_plk1_rate(config_dir="config"):
    """Try to read the base $d_Plk1 rate from MaBoSS config files"""
    import glob
    
    # Look for .cfg files
    cfg_files = []
    for pattern in ["*.cfg", "**/*.cfg"]:
        cfg_files.extend(glob.glob(os.path.join(config_dir, pattern), recursive=True))
    
    base_rate = None
    for cfg_file in cfg_files:
        try:
            with open(cfg_file, 'r') as f:
                for line in f:
                    # Look for $d_Plk1 parameter definition (format: $d_Plk1 = 1;)
                    line_stripped = line.strip()
                    if '$d_Plk1' in line_stripped and '=' in line_stripped:
                        # Try to extract value (format: $d_Plk1 = value; or $d_Plk1 = value)
                        parts = line_stripped.split('=')
                        if len(parts) > 1:
                            value_str = parts[1].strip().rstrip(';').strip()
                            try:
                                base_rate = float(value_str)
                                print(f"Found base $d_Plk1 rate: {base_rate} in {cfg_file}")
                                return base_rate
                            except ValueError:
                                continue
        except Exception as e:
            continue
    
    # Default fallback - based on config files found, default is 1.0
    print("Warning: Could not find base $d_Plk1 rate in config files. Using default: 1.0")
    return 1.0  # Default from config files

def calculate_actual_rate_values(plk1_stats, base_rate):
    """Calculate actual rate values from mutations"""
    # Track rate for each cell over time
    cell_rates = {}  # {cell_id: [(time, rate), ...]}
    
    for mutation_detail in plk1_stats['mutation_details']:
        cell_id = mutation_detail['cell_id']
        time = mutation_detail['time']
        effect_size = mutation_detail['effect_size']
        
        if cell_id not in cell_rates:
            cell_rates[cell_id] = []
            # Initialize with base rate at time 0
            cell_rates[cell_id].append((0, base_rate))
        
        # Get current rate (last rate for this cell)
        current_rate = cell_rates[cell_id][-1][1]
        
        # Calculate new rate: new_rate = current_rate * effect_size
        new_rate = current_rate * effect_size
        
        # Add the new rate at this time
        cell_rates[cell_id].append((time, new_rate))
    
    return cell_rates

def plot_effect_size_distribution(plk1_stats, output_path):
    """Plot distribution of effect sizes"""
    
    # Collect all effect sizes
    all_effect_sizes = []
    for time, sizes in plk1_stats['effect_sizes'].items():
        all_effect_sizes.extend(sizes)
    
    if not all_effect_sizes:
        print("No effect size data available for plotting")
        return
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Plot 1: Histogram
    ax1 = axes[0]
    ax1.hist(all_effect_sizes, bins=50, edgecolor='black', alpha=0.7)
    ax1.set_xlabel('Effect Size (Multiplier)', fontsize=12)
    ax1.set_ylabel('Frequency', fontsize=12)
    ax1.set_title('Distribution of PLK1 Rate Effect Sizes\n(Note: This shows the multiplier, not actual rates)', fontsize=14)
    ax1.grid(True, alpha=0.3)
    ax1.axvline(x=1.0, color='r', linestyle='--', linewidth=2, 
                label='Effect size = 1.0 (no change!)')
    ax1.legend()
    
    # Plot 2: Box plot
    ax2 = axes[1]
    ax2.boxplot(all_effect_sizes, vert=True)
    ax2.set_ylabel('Effect Size (Multiplier)', fontsize=12)
    ax2.set_title('PLK1 Rate Effect Size Statistics', fontsize=14)
    ax2.grid(True, alpha=0.3)
    ax2.axhline(y=1.0, color='r', linestyle='--', linewidth=2, 
                label='Effect size = 1.0 (no change!)')
    ax2.legend()
    
    plt.tight_layout()
    plt.savefig(f"{output_path}_effect_sizes.png", bbox_inches='tight', dpi=300)
    plt.savefig(f"{output_path}_effect_sizes.pdf", bbox_inches='tight')
    plt.close()
    
    print(f"Effect size distribution plot saved to: {output_path}_effect_sizes.png/pdf")
    
    # Print statistics
    print(f"\nEffect Size Statistics:")
    print(f"  Mean: {np.mean(all_effect_sizes):.6f}")
    print(f"  Median: {np.median(all_effect_sizes):.6f}")
    print(f"  Min: {np.min(all_effect_sizes):.6f}")
    print(f"  Max: {np.max(all_effect_sizes):.6f}")
    print(f"  Std: {np.std(all_effect_sizes):.6f}")
    
    # Warning if effect size is 1.0
    if np.allclose(all_effect_sizes, 1.0, atol=1e-6):
        print("\n⚠️  WARNING: All effect sizes are 1.0!")
        print("   This means: new_rate = current_rate * 1.0 = current_rate")
        print("   The PLK1 rates are NOT actually changing!")
        print("   Check plk1_rate_effect_size in PhysiCell_settings.xml")

def plot_actual_rate_values(plk1_stats, time_points, output_path, base_rate):
    """Plot actual rate values over time (calculated from mutations)"""
    
    # Calculate actual rate values
    cell_rates = calculate_actual_rate_values(plk1_stats, base_rate)
    
    if not cell_rates:
        print("No rate data available for plotting")
        return
    
    # Collect all rate values for statistics
    all_rates = []
    for cell_id, rate_history in cell_rates.items():
        for time, rate in rate_history:
            all_rates.append(rate)
    
    if not all_rates:
        return
    
    # Create plots
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    # Plot 1: Distribution of actual rate values
    ax1 = axes[0, 0]
    ax1.hist(all_rates, bins=50, edgecolor='black', alpha=0.7, color='green')
    ax1.set_xlabel('Actual $d_Plk1 Rate Value', fontsize=12)
    ax1.set_ylabel('Frequency', fontsize=12)
    ax1.set_title('Distribution of Actual PLK1 Down Rate Values\n(Shows rate changes from mutations)', fontsize=14)
    ax1.grid(True, alpha=0.3)
    ax1.axvline(x=base_rate, color='r', linestyle='--', linewidth=2, 
                label=f'Base rate = {base_rate}')
    ax1.legend()
    
    # Plot 2: Rate values over time (sample cells)
    ax2 = axes[0, 1]
    # Sample some cells to plot
    sample_cells = list(cell_rates.keys())[:20]  # Plot first 20 cells
    colors = plt.cm.viridis(np.linspace(0, 1, len(sample_cells)))
    
    for i, cell_id in enumerate(sample_cells):
        rate_history = cell_rates[cell_id]
        if len(rate_history) > 1:  # Only plot if cell has mutations
            times = [t for t, r in rate_history]
            rates = [r for t, r in rate_history]
            ax2.plot(times, rates, 'o-', color=colors[i], alpha=0.6, 
                    linewidth=1.5, markersize=3, label=f'Cell {cell_id}' if i < 5 else '')
    
    ax2.axhline(y=base_rate, color='r', linestyle='--', linewidth=2, 
                label=f'Base rate = {base_rate}')
    ax2.set_xlabel('Time Step', fontsize=12)
    ax2.set_ylabel('Actual $d_Plk1 Rate Value', fontsize=12)
    ax2.set_title('PLK1 Rate Values Over Time (Sample Cells)', fontsize=14)
    ax2.grid(True, alpha=0.3)
    if len(sample_cells) <= 5:
        ax2.legend()
    
    # Plot 3: Average rate over time
    ax3 = axes[1, 0]
    # Calculate average rate at each time point
    time_rate_avg = {}
    time_rate_counts = {}
    
    for cell_id, rate_history in cell_rates.items():
        for time, rate in rate_history:
            if time not in time_rate_avg:
                time_rate_avg[time] = 0.0
                time_rate_counts[time] = 0
            time_rate_avg[time] += rate
            time_rate_counts[time] += 1
    
    # Calculate averages
    times_avg = sorted(time_rate_avg.keys())
    avg_rates = [time_rate_avg[t] / time_rate_counts[t] for t in times_avg]
    
    ax3.plot(times_avg, avg_rates, 'b-', linewidth=2, marker='o', markersize=4)
    ax3.axhline(y=base_rate, color='r', linestyle='--', linewidth=2, 
                label=f'Base rate = {base_rate}')
    ax3.set_xlabel('Time Step', fontsize=12)
    ax3.set_ylabel('Average $d_Plk1 Rate Value', fontsize=12)
    ax3.set_title('Average PLK1 Rate Over Time', fontsize=14)
    ax3.grid(True, alpha=0.3)
    ax3.legend()
    
    # Plot 4: Rate statistics
    ax4 = axes[1, 1]
    rate_stats = {
        'Min': np.min(all_rates),
        '25th %ile': np.percentile(all_rates, 25),
        'Median': np.median(all_rates),
        '75th %ile': np.percentile(all_rates, 75),
        'Max': np.max(all_rates),
        'Base': base_rate
    }
    
    bars = ax4.bar(range(len(rate_stats)), list(rate_stats.values()), 
                   color=['blue', 'cyan', 'green', 'orange', 'red', 'purple'])
    ax4.set_xticks(range(len(rate_stats)))
    ax4.set_xticklabels(list(rate_stats.keys()), rotation=45, ha='right')
    ax4.set_ylabel('$d_Plk1 Rate Value', fontsize=12)
    ax4.set_title('PLK1 Rate Value Statistics', fontsize=14)
    ax4.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.savefig(f"{output_path}_actual_rates.png", bbox_inches='tight', dpi=300)
    plt.savefig(f"{output_path}_actual_rates.pdf", bbox_inches='tight')
    plt.close()
    
    print(f"Actual rate values plot saved to: {output_path}_actual_rates.png/pdf")
    
    # Print statistics
    print(f"\nActual Rate Value Statistics:")
    print(f"  Base rate: {base_rate:.6f}")
    print(f"  Min: {np.min(all_rates):.6f}")
    print(f"  Max: {np.max(all_rates):.6f}")
    print(f"  Mean: {np.mean(all_rates):.6f}")
    print(f"  Median: {np.median(all_rates):.6f}")
    print(f"  Std: {np.std(all_rates):.6f}")
    print(f"  Fold change (max/base): {np.max(all_rates)/base_rate:.2f}x")

def create_mutation_summary_table(plk1_stats, time_points, output_path):
    """Create a summary table of PLK1 mutations"""
    
    if not plk1_stats['mutation_details']:
        print("No mutation details available for summary table")
        return
    
    # Create DataFrame from mutation details
    df = pd.DataFrame(plk1_stats['mutation_details'])
    
    # Save to CSV
    csv_path = f"{output_path}_mutation_details.csv"
    df.to_csv(csv_path, index=False)
    print(f"Mutation details table saved to: {csv_path}")
    
    # Create summary statistics
    summary_data = []
    for time in time_points:
        time_mutations = [m for m in plk1_stats['mutation_details'] if m['time'] == time]
        if time_mutations:
            summary_data.append({
                'Time': time,
                'Cells_with_PLK1_mutations': len(set(m['cell_id'] for m in time_mutations)),
                'Total_PLK1_mutations': len(time_mutations),
                'Up_mutations': sum(1 for m in time_mutations if m['type'] == 'up'),
                'Down_mutations': sum(1 for m in time_mutations if m['type'] == 'down'),
                'Mean_effect_size': np.mean([m['effect_size'] for m in time_mutations]),
                'Min_effect_size': np.min([m['effect_size'] for m in time_mutations]),
                'Max_effect_size': np.max([m['effect_size'] for m in time_mutations])
            })
    
    if summary_data:
        summary_df = pd.DataFrame(summary_data)
        summary_csv_path = f"{output_path}_summary.csv"
        summary_df.to_csv(summary_csv_path, index=False)
        print(f"Summary table saved to: {summary_csv_path}")
    
    return df

def analyze_parameter_changes(plk1_stats, time_points):
    """Analyze how parameter values would change based on mutations"""
    
    # Try to read the base parameter value from config
    # For now, we'll simulate the changes assuming a base value
    # In reality, we'd need to read this from the MaBoSS model files
    
    base_rate = 1.0  # Default assumption - would need to read from .bnd/.cfg files
    
    parameter_trajectories = {}
    
    for mutation_detail in plk1_stats['mutation_details']:
        cell_id = mutation_detail['cell_id']
        time = mutation_detail['time']
        effect_size = mutation_detail['effect_size']
        mut_type = mutation_detail['type']
        
        if cell_id not in parameter_trajectories:
            parameter_trajectories[cell_id] = {
                'times': [],
                'values': [],
                'mutations': []
            }
        
        # Calculate new value: new_rate = old_rate * effect_size
        if len(parameter_trajectories[cell_id]['values']) == 0:
            current_value = base_rate
        else:
            current_value = parameter_trajectories[cell_id]['values'][-1]
        
        new_value = current_value * effect_size
        parameter_trajectories[cell_id]['times'].append(time)
        parameter_trajectories[cell_id]['values'].append(new_value)
        parameter_trajectories[cell_id]['mutations'].append(mutation_detail)
    
    return parameter_trajectories

def plot_parameter_trajectories(parameter_trajectories, output_path, max_cells=50):
    """Plot parameter value trajectories for individual cells"""
    
    if not parameter_trajectories:
        print("No parameter trajectories available for plotting")
        return
    
    # Select a subset of cells to plot (to avoid overcrowding)
    cell_ids = list(parameter_trajectories.keys())
    if len(cell_ids) > max_cells:
        # Sample cells
        import random
        random.seed(42)  # For reproducibility
        cell_ids = random.sample(cell_ids, max_cells)
        print(f"Plotting {max_cells} randomly selected cells out of {len(parameter_trajectories)} total")
    
    fig, ax = plt.subplots(figsize=(14, 8))
    
    colors = plt.cm.viridis(np.linspace(0, 1, len(cell_ids)))
    
    for i, cell_id in enumerate(cell_ids):
        traj = parameter_trajectories[cell_id]
        if traj['times']:
            ax.plot(traj['times'], traj['values'], 
                   color=colors[i], alpha=0.6, linewidth=1.5,
                   label=f'Cell {cell_id}' if i < 10 else '')
    
    ax.axhline(y=1.0, color='r', linestyle='--', linewidth=2, 
               label='Base rate (1.0) - No change if effect_size=1.0')
    ax.set_xlabel('Time Step', fontsize=12)
    ax.set_ylabel('Parameter Value (relative to base)', fontsize=12)
    ax.set_title(f'PLK1 Rate Parameter Trajectories (showing {len(cell_ids)} cells)', fontsize=14)
    ax.grid(True, alpha=0.3)
    if len(cell_ids) <= 10:
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
    
    plt.tight_layout()
    plt.savefig(f"{output_path}_parameter_trajectories.png", bbox_inches='tight', dpi=300)
    plt.savefig(f"{output_path}_parameter_trajectories.pdf", bbox_inches='tight')
    plt.close()
    
    print(f"Parameter trajectories plot saved to: {output_path}_parameter_trajectories.png/pdf")

def plot_mutation_frequency(plk1_stats, time_points, output_path):
    """Plot frequency of mutations per cell over time"""
    
    if not plk1_stats['mutation_details']:
        return
    
    # Calculate mutations per cell at each time point
    mutations_per_cell = []
    for time in time_points:
        time_mutations = [m for m in plk1_stats['mutation_details'] if m['time'] == time]
        cells_with_mutations = len(set(m['cell_id'] for m in time_mutations))
        total_mutations = len(time_mutations)
        
        if cells_with_mutations > 0:
            avg_mutations = total_mutations / cells_with_mutations
        else:
            avg_mutations = 0
        
        mutations_per_cell.append(avg_mutations)
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 10))
    
    # Plot 1: Average mutations per cell
    ax1.plot(time_points, mutations_per_cell, 'b-', linewidth=2, marker='o', markersize=4)
    ax1.set_xlabel('Time Step', fontsize=12)
    ax1.set_ylabel('Average PLK1 Mutations per Cell', fontsize=12)
    ax1.set_title('Average Number of PLK1 Rate Mutations per Cell Over Time', fontsize=14)
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Cumulative mutations
    cumulative = np.cumsum(plk1_stats['total_plk1_mutations'])
    ax2.plot(time_points, cumulative, 'g-', linewidth=2, marker='s', markersize=4)
    ax2.set_xlabel('Time Step', fontsize=12)
    ax2.set_ylabel('Cumulative PLK1 Mutations', fontsize=12)
    ax2.set_title('Cumulative PLK1 Rate Mutations Over Time', fontsize=14)
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f"{output_path}_mutation_frequency.png", bbox_inches='tight', dpi=300)
    plt.savefig(f"{output_path}_mutation_frequency.pdf", bbox_inches='tight')
    plt.close()
    
    print(f"Mutation frequency plot saved to: {output_path}_mutation_frequency.png/pdf")

def print_diagnostic_info(plk1_stats, time_points):
    """Print diagnostic information"""
    
    print("\n" + "="*60)
    print("PLK1 RATE MUTATION DIAGNOSTICS")
    print("="*60)
    
    total_cells_with_plk1 = sum(plk1_stats['cells_with_plk1_mutations'])
    total_plk1_mutations = sum(plk1_stats['total_plk1_mutations'])
    
    print(f"\nOverall Statistics:")
    print(f"  Time range: {min(time_points)} - {max(time_points)} steps")
    print(f"  Total time points analyzed: {len(time_points)}")
    print(f"  Total cells with PLK1 mutations (cumulative): {total_cells_with_plk1}")
    print(f"  Total PLK1 mutations (cumulative): {total_plk1_mutations}")
    
    if plk1_stats['first_mutation_time'] is not None:
        print(f"  First PLK1 mutation occurred at time: {plk1_stats['first_mutation_time']}")
    else:
        print(f"  ⚠️  WARNING: No PLK1 mutations found in the simulation!")
    
    # Check effect sizes
    all_effect_sizes = []
    for sizes in plk1_stats['effect_sizes'].values():
        all_effect_sizes.extend(sizes)
    
    if all_effect_sizes:
        unique_effect_sizes = set(all_effect_sizes)
        print(f"\nEffect Size Analysis:")
        print(f"  Unique effect sizes found: {sorted(unique_effect_sizes)}")
        if len(unique_effect_sizes) == 1 and 1.0 in unique_effect_sizes:
            print(f"  ⚠️  CRITICAL: Only effect size = 1.0 found!")
            print(f"     This means PLK1 rates are NOT changing (multiplied by 1.0)")
            print(f"     Check plk1_rate_effect_size in PhysiCell_settings.xml")
            print(f"     Expected: effect_size != 1.0 (e.g., 0.5 for 50% reduction, 2.0 for doubling)")
    
    # Mutation type distribution
    total_up = sum(plk1_stats['plk1_up_mutations'])
    total_down = sum(plk1_stats['plk1_down_mutations'])
    print(f"\nMutation Type Distribution:")
    print(f"  UP mutations: {total_up}")
    print(f"  DOWN mutations: {total_down}")
    print(f"  Target parameter: $d_Plk1 (down/inhibition rate)")
    
    # Calculate mutation rate
    if len(time_points) > 1:
        time_span = max(time_points) - min(time_points)
        mutations_per_step = total_plk1_mutations / len(time_points) if time_points else 0
        print(f"\nMutation Rate:")
        print(f"  Average mutations per time step: {mutations_per_step:.2f}")
        print(f"  Mutation probability setting: 0.5 (from config)")
    
    # Check configuration
    print(f"\nConfiguration Check:")
    print(f"  ⚠️  Please verify PhysiCell_settings.xml:")
    print(f"     - plk1_rate_mutation_probability: should be > 0 (currently 0.5)")
    print(f"     - plk1_rate_effect_size: should NOT be 1.0 (currently 1.0 - THIS IS THE PROBLEM!)")
    print(f"     - plk1_rate_mutation_target: 'down' (affecting $d_Plk1)")
    print(f"     - Recommended: Set effect_size to 0.5 (50% reduction) or 2.0 (doubling)")
    
    print("="*60)

# Main execution
if __name__ == "__main__":
    print("=== PLK1 Rate Mutation Debugging Analysis ===")
    
    # Analyze PLK1 mutations
    plk1_stats, time_points = analyze_plk1_mutations(output_dir)
    
    if time_points:
        # Set output path
        base_path = os.path.join(debug_dir, "plk1")
        
        # Get base rate from config
        config_dir = os.path.join(os.path.dirname(os.path.dirname(script_dir)), "config")
        base_rate = get_base_plk1_rate(config_dir)
        
        # Create plots
        plot_plk1_mutation_timeline(plk1_stats, time_points, base_path)
        plot_effect_size_distribution(plk1_stats, base_path)
        plot_actual_rate_values(plk1_stats, time_points, base_path, base_rate)
        plot_mutation_frequency(plk1_stats, time_points, base_path)
        
        # Analyze parameter trajectories
        parameter_trajectories = analyze_parameter_changes(plk1_stats, time_points)
        if parameter_trajectories:
            plot_parameter_trajectories(parameter_trajectories, base_path)
        
        # Create summary tables
        mutation_df = create_mutation_summary_table(plk1_stats, time_points, base_path)
        
        # Print diagnostics
        print_diagnostic_info(plk1_stats, time_points)
        
        print(f"\nAll PLK1 debugging plots saved to: {debug_dir}")
    else:
        print("Failed to load mutation data. Please check:")
        print("1. Output directory contains valid .mut.csv files")
        print("2. Files are properly formatted")

