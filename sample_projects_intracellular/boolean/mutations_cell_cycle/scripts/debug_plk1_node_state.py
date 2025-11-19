#!/usr/bin/env python3
"""
Debugging script to track PLK1 node state (Boolean ON/OFF) over time.
This script analyzes the MultiCellDS output files to track the actual
PLK1 node state and its relationship to rate mutations.
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
debug_dir = os.path.join(results_dir, "plk1_debug")

# Create directories if they don't exist
os.makedirs(debug_dir, exist_ok=True)
print(f"Results will be saved to: {debug_dir}")

def analyze_plk1_node_states_pcdl(output_dir):
    """Analyze PLK1 node states using pcdl"""
    
    if not PCDL_AVAILABLE:
        print("pcdl not available. Cannot analyze node states.")
        return None, None
    
    try:
        print("Loading simulation data with pcdl...")
        ts = TimeSeries(output_dir)
        
        try:
            mcds_list = ts.get_mcds_list()
        except Exception as e:
            print(f"Error loading MultiCellDS time steps: {e}")
            return None, None

        if not mcds_list:
            print("No MultiCellDS snapshots found in the output directory.")
            return None, None

        time_points = [float(step.get_time()) for step in mcds_list]
        print(f"Found {len(time_points)} time steps")
        
        # Initialize data structures
        node_state_data = {
            'times': [],
            'cells_with_plk1_on': [],
            'cells_with_plk1_off': [],
            'total_cells': [],
            'plk1_on_percentage': [],
            'cell_states': {}  # {cell_id: [(time, plk1_state), ...]}
        }
        
        # Process each time point
        sorted_indices = sorted(range(len(time_points)), key=lambda i: time_points[i])
        for idx in sorted_indices:
            timestep = time_points[idx]
            timestep_data = mcds_list[idx]
            
            try:
                print(f"Processing time step {timestep}...")
                
                # Get cell data
                cell_df = timestep_data.get_cell_df()
                
                if cell_df.empty:
                    print(f"  No cell data found at time step {timestep}")
                    continue
                
                # Check for PLK1 node state
                plk1_on_count = 0
                plk1_off_count = 0
                
                if 'Plk1' in cell_df.columns:
                    # PLK1 node state is directly in the dataframe
                    plk1_states = cell_df['Plk1']
                    plk1_on_count = (plk1_states == True).sum() if plk1_states.dtype == bool else (plk1_states == 1).sum()
                    plk1_off_count = len(cell_df) - plk1_on_count
                    
                    # Store individual cell states
                    for _, row in cell_df.iterrows():
                        cell_id = int(row.get('ID', -1))
                        if cell_id >= 0:
                            plk1_state = bool(row['Plk1']) if row['Plk1'].dtype == bool else bool(row['Plk1'] == 1)
                            if cell_id not in node_state_data['cell_states']:
                                node_state_data['cell_states'][cell_id] = []
                            node_state_data['cell_states'][cell_id].append((timestep, plk1_state))
                
                elif 'boolean_intracellular' in cell_df.columns:
                    # Try to parse boolean states from intracellular data
                    print(f"  Found boolean_intracellular column, attempting to parse...")
                    # This would require parsing the intracellular state string
                    # For now, skip this approach
                    continue
                else:
                    print(f"  PLK1 node state not found in columns: {list(cell_df.columns)[:10]}...")
                    # Try to get it from intracellular model
                    print(f"  Attempting to access intracellular model directly...")
                    continue
                
                total_cells = len(cell_df)
                plk1_on_pct = (plk1_on_count / total_cells * 100) if total_cells > 0 else 0
                
                node_state_data['times'].append(timestep)
                node_state_data['cells_with_plk1_on'].append(plk1_on_count)
                node_state_data['cells_with_plk1_off'].append(plk1_off_count)
                node_state_data['total_cells'].append(total_cells)
                node_state_data['plk1_on_percentage'].append(plk1_on_pct)
                
                if idx % 10 == 0 or plk1_on_count > 0:
                    print(f"  Time {timestep}: {plk1_on_count} cells with PLK1 ON, "
                          f"{plk1_off_count} cells with PLK1 OFF ({plk1_on_pct:.1f}% ON)")
                
            except Exception as e:
                print(f"Error processing time step {timestep}: {e}")
                import traceback
                traceback.print_exc()
                continue
        
        if not node_state_data['times']:
            print("No valid PLK1 node state data found!")
            return None, None
        
        print(f"Successfully processed {len(node_state_data['times'])} time points")
        return node_state_data, node_state_data['times']
        
    except Exception as e:
        print(f"Error loading simulation data: {e}")
        import traceback
        traceback.print_exc()
        return None, None

def analyze_plk1_node_states_from_csv(output_dir):
    """Alternative: Try to read PLK1 node states from mut.csv files if available"""
    
    # Look for mut.csv files (these contain custom data including Plk1)
    csv_files = [f for f in os.listdir(output_dir) if f.endswith('_mut.csv')]
    csv_files = sorted([f for f in csv_files if 'output' in f])  # Only output files, sorted
    
    if not csv_files:
        print("No mut.csv files found")
        return None, None
    
    print(f"Found {len(csv_files)} boolean_intracellular CSV files")
    
    node_state_data = {
        'times': [],
        'cells_with_plk1_on': [],
        'cells_with_plk1_off': [],
        'total_cells': [],
        'plk1_on_percentage': [],
        'cell_states': {}
    }
    
    for csv_file in csv_files:
        try:
            # Extract time from filename (format: output00000080_mut.csv)
            time_match = re.search(r'output(\d+)_mut\.csv', csv_file)
            if not time_match:
                continue
            time_step = int(time_match.group(1))
            
            df = pd.read_csv(os.path.join(output_dir, csv_file))
            
            # Check if PLK1 column exists
            if 'Plk1' not in df.columns:
                if len(csv_files) > 0 and csv_file == csv_files[0]:
                    print(f"Warning: 'Plk1' column not found in {csv_file}. Available columns: {list(df.columns)}")
                continue
            
            plk1_states = df['Plk1']
            # Handle both boolean and numeric (0.0/1.0) formats
            if plk1_states.dtype == bool:
                plk1_on_count = plk1_states.sum()
            else:
                # Treat values >= 0.5 as ON, < 0.5 as OFF
                plk1_on_count = (plk1_states >= 0.5).sum()
            plk1_off_count = len(df) - plk1_on_count
            total_cells = len(df)
            plk1_on_pct = (plk1_on_count / total_cells * 100) if total_cells > 0 else 0
            
            node_state_data['times'].append(time_step)
            node_state_data['cells_with_plk1_on'].append(plk1_on_count)
            node_state_data['cells_with_plk1_off'].append(plk1_off_count)
            node_state_data['total_cells'].append(total_cells)
            node_state_data['plk1_on_percentage'].append(plk1_on_pct)
            
            # Store individual cell states
            if 'ID' in df.columns:
                for _, row in df.iterrows():
                    cell_id = int(row['ID'])
                    # Handle both boolean and numeric formats
                    if isinstance(row['Plk1'], bool):
                        plk1_state = row['Plk1']
                    else:
                        plk1_state = bool(row['Plk1'] >= 0.5)
                    if cell_id not in node_state_data['cell_states']:
                        node_state_data['cell_states'][cell_id] = []
                    node_state_data['cell_states'][cell_id].append((time_step, plk1_state))
            
            if len(node_state_data['times']) % 10 == 0:
                print(f"  Processed time step {time_step}: {plk1_on_count} cells with PLK1 ON ({plk1_on_pct:.1f}%)")
                
        except Exception as e:
            print(f"Error processing {csv_file}: {e}")
            continue
    
    if not node_state_data['times']:
        return None, None
    
    return node_state_data, node_state_data['times']

def plot_plk1_node_state_timeline(node_state_data, time_points, output_path):
    """Plot PLK1 node state timeline"""
    
    if not node_state_data or not time_points:
        print("No node state data available for plotting")
        return
    
    fig, axes = plt.subplots(2, 1, figsize=(14, 10))
    
    # Plot 1: Cell counts
    ax1 = axes[0]
    ax1.plot(time_points, node_state_data['cells_with_plk1_on'], 
             'g-', linewidth=2, marker='o', markersize=4, label='PLK1 ON')
    ax1.plot(time_points, node_state_data['cells_with_plk1_off'], 
             'r-', linewidth=2, marker='s', markersize=4, label='PLK1 OFF')
    ax1.set_xlabel('Time Step', fontsize=12)
    ax1.set_ylabel('Number of Cells', fontsize=12)
    ax1.set_title('PLK1 Node State Distribution Over Time', fontsize=14)
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Percentage
    ax2 = axes[1]
    ax2.plot(time_points, node_state_data['plk1_on_percentage'], 
             'b-', linewidth=2, marker='o', markersize=4, label='% Cells with PLK1 ON')
    ax2.set_xlabel('Time Step', fontsize=12)
    ax2.set_ylabel('Percentage of Cells (%)', fontsize=12)
    ax2.set_title('Percentage of Cells with PLK1 Node ON Over Time', fontsize=14)
    ax2.set_ylim(0, 100)
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f"{output_path}_node_state_timeline.png", bbox_inches='tight', dpi=300)
    plt.savefig(f"{output_path}_node_state_timeline.pdf", bbox_inches='tight')
    plt.close()
    
    print(f"Node state timeline plot saved to: {output_path}_node_state_timeline.png/pdf")

def plot_plk1_node_state_trajectories(node_state_data, output_path, max_cells=50):
    """Plot PLK1 node state trajectories for individual cells"""
    
    if not node_state_data or not node_state_data['cell_states']:
        print("No individual cell state data available for plotting")
        return
    
    cell_ids = list(node_state_data['cell_states'].keys())
    if len(cell_ids) > max_cells:
        import random
        random.seed(42)
        cell_ids = random.sample(cell_ids, max_cells)
        print(f"Plotting {max_cells} randomly selected cells out of {len(node_state_data['cell_states'])} total")
    
    fig, ax = plt.subplots(figsize=(14, 8))
    
    colors = plt.cm.viridis(np.linspace(0, 1, len(cell_ids)))
    
    for i, cell_id in enumerate(cell_ids):
        state_history = node_state_data['cell_states'][cell_id]
        if state_history:
            times = [t for t, s in state_history]
            states = [1.0 if s else 0.0 for t, s in state_history]  # Convert bool to 0/1
            ax.plot(times, states, 'o-', color=colors[i], alpha=0.6, 
                   linewidth=1.5, markersize=2, label=f'Cell {cell_id}' if i < 10 else '')
    
    ax.axhline(y=0.5, color='k', linestyle='--', linewidth=1, alpha=0.3)
    ax.set_xlabel('Time Step', fontsize=12)
    ax.set_ylabel('PLK1 Node State (0=OFF, 1=ON)', fontsize=12)
    ax.set_title(f'PLK1 Node State Trajectories (showing {len(cell_ids)} cells)', fontsize=14)
    ax.set_ylim(-0.1, 1.1)
    ax.set_yticks([0, 1])
    ax.set_yticklabels(['OFF', 'ON'])
    ax.grid(True, alpha=0.3)
    if len(cell_ids) <= 10:
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
    
    plt.tight_layout()
    plt.savefig(f"{output_path}_node_state_trajectories.png", bbox_inches='tight', dpi=300)
    plt.savefig(f"{output_path}_node_state_trajectories.pdf", bbox_inches='tight')
    plt.close()
    
    print(f"Node state trajectories plot saved to: {output_path}_node_state_trajectories.png/pdf")

def correlate_plk1_rate_and_state(plk1_stats, node_state_data, time_points, output_path):
    """Correlate PLK1 rate mutations with PLK1 node state"""
    
    if not plk1_stats or not node_state_data:
        print("Insufficient data for correlation analysis")
        return
    
    # Match time points between mutation data and node state data
    # Mutation data uses time indices (0, 1, 2, ...)
    # Node state data uses actual time values
    
    # Create correlation data
    correlation_data = []
    
    # Get mutation times
    mutation_times = set()
    for mutation_detail in plk1_stats['mutation_details']:
        mutation_times.add(mutation_detail['time'])
    
    # Match with node state data
    for time_idx in sorted(mutation_times):
        # Find corresponding node state time
        node_state_idx = None
        for i, t in enumerate(node_state_data['times']):
            if abs(t - time_idx * 100) < 50:  # Approximate matching (time_idx * 100 = actual time)
                node_state_idx = i
                break
        
        if node_state_idx is not None:
            plk1_on_pct = node_state_data['plk1_on_percentage'][node_state_idx]
            cells_with_mutations = sum(1 for m in plk1_stats['mutation_details'] if m['time'] == time_idx)
            
            correlation_data.append({
                'time': time_idx,
                'plk1_on_percentage': plk1_on_pct,
                'cells_with_mutations': cells_with_mutations
            })
    
    if not correlation_data:
        print("No matching time points for correlation")
        return
    
    df_corr = pd.DataFrame(correlation_data)
    
    # Create plot
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 10))
    
    # Plot 1: Overlay
    ax1_twin = ax1.twinx()
    
    line1 = ax1.plot(df_corr['time'], df_corr['plk1_on_percentage'], 
                     'b-', linewidth=2, marker='o', markersize=4, label='% Cells with PLK1 ON')
    line2 = ax1_twin.plot(df_corr['time'], df_corr['cells_with_mutations'], 
                          'r-', linewidth=2, marker='s', markersize=4, label='Cells with PLK1 Rate Mutations')
    
    ax1.set_xlabel('Time Step', fontsize=12)
    ax1.set_ylabel('PLK1 ON Percentage (%)', fontsize=12, color='b')
    ax1_twin.set_ylabel('Cells with PLK1 Rate Mutations', fontsize=12, color='r')
    ax1.set_title('PLK1 Node State vs Rate Mutations Over Time', fontsize=14)
    ax1.tick_params(axis='y', labelcolor='b')
    ax1_twin.tick_params(axis='y', labelcolor='r')
    ax1.grid(True, alpha=0.3)
    
    # Combined legend
    lines = line1 + line2
    labels = [l.get_label() for l in lines]
    ax1.legend(lines, labels, loc='upper left')
    
    # Plot 2: Scatter correlation
    if len(df_corr) > 1:
        ax2.scatter(df_corr['cells_with_mutations'], df_corr['plk1_on_percentage'], 
                   alpha=0.6, s=50, edgecolors='black')
        ax2.set_xlabel('Cells with PLK1 Rate Mutations', fontsize=12)
        ax2.set_ylabel('PLK1 ON Percentage (%)', fontsize=12)
        ax2.set_title('Correlation: PLK1 Rate Mutations vs Node State', fontsize=14)
        ax2.grid(True, alpha=0.3)
        
        # Calculate correlation coefficient
        if len(df_corr) > 2:
            corr_coef = np.corrcoef(df_corr['cells_with_mutations'], df_corr['plk1_on_percentage'])[0, 1]
            ax2.text(0.05, 0.95, f'Correlation: {corr_coef:.3f}', 
                    transform=ax2.transAxes, fontsize=12,
                    verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout()
    plt.savefig(f"{output_path}_rate_state_correlation.png", bbox_inches='tight', dpi=300)
    plt.savefig(f"{output_path}_rate_state_correlation.pdf", bbox_inches='tight')
    plt.close()
    
    print(f"Rate-state correlation plot saved to: {output_path}_rate_state_correlation.png/pdf")

def print_node_state_summary(node_state_data, time_points):
    """Print summary of PLK1 node state data"""
    
    if not node_state_data or not time_points:
        return
    
    print("\n" + "="*60)
    print("PLK1 NODE STATE SUMMARY")
    print("="*60)
    
    print(f"\nTime Range:")
    print(f"  {min(time_points):.1f} - {max(time_points):.1f} minutes")
    print(f"  {len(time_points)} time points")
    
    print(f"\nPLK1 Node State Statistics:")
    avg_on_pct = np.mean(node_state_data['plk1_on_percentage'])
    max_on_pct = np.max(node_state_data['plk1_on_percentage'])
    min_on_pct = np.min(node_state_data['plk1_on_percentage'])
    
    print(f"  Average % cells with PLK1 ON: {avg_on_pct:.2f}%")
    print(f"  Maximum % cells with PLK1 ON: {max_on_pct:.2f}%")
    print(f"  Minimum % cells with PLK1 ON: {min_on_pct:.2f}%")
    
    if node_state_data['cell_states']:
        print(f"\nIndividual Cell Tracking:")
        print(f"  Cells tracked: {len(node_state_data['cell_states'])}")
        
        # Count cells that switch states
        cells_switching = 0
        for cell_id, state_history in node_state_data['cell_states'].items():
            if len(state_history) > 1:
                states = [s for t, s in state_history]
                if len(set(states)) > 1:  # Has both ON and OFF
                    cells_switching += 1
        
        print(f"  Cells that switch PLK1 state: {cells_switching}")
    
    print("="*60)

# Main execution
if __name__ == "__main__":
    print("=== PLK1 Node State Debugging Analysis ===")
    
    # Try CSV files first (faster and more direct)
    print("Trying to read PLK1 node state from mut.csv files...")
    node_state_data, time_points = analyze_plk1_node_states_from_csv(output_dir)
    
    # Fallback to pcdl if CSV doesn't work
    if not node_state_data:
        print("\nCSV method failed. Trying pcdl method...")
        node_state_data, time_points = analyze_plk1_node_states_pcdl(output_dir)
    
    if node_state_data and time_points:
        # Set output path
        base_path = os.path.join(debug_dir, "plk1_node_state")
        
        # Create plots
        plot_plk1_node_state_timeline(node_state_data, time_points, base_path)
        plot_plk1_node_state_trajectories(node_state_data, base_path)
        
        # Try to correlate with rate mutations
        try:
            # Import the mutation analysis
            import sys
            sys.path.append(os.path.dirname(__file__))
            from debug_plk1_rates import analyze_plk1_mutations
            
            print("\nAnalyzing PLK1 rate mutations for correlation...")
            plk1_stats, mutation_times = analyze_plk1_mutations(output_dir)
            if plk1_stats:
                correlate_plk1_rate_and_state(plk1_stats, node_state_data, time_points, base_path)
        except Exception as e:
            print(f"Could not correlate with rate mutations: {e}")
        
        # Print summary
        print_node_state_summary(node_state_data, time_points)
        
        print(f"\nAll PLK1 node state debugging plots saved to: {debug_dir}")
    else:
        print("Failed to load PLK1 node state data. Please check:")
        print("1. pcdl package is installed: pip install pcdl")
        print("2. Output directory contains valid MultiCellDS data")
        print("3. PLK1 node state is available in the data")

