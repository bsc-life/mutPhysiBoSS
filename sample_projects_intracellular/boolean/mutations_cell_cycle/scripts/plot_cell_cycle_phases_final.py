# this script will work with pcdl v4.0.5 to analyze cell cycle phases
# It uses the correct API methods for this version

import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pandas as pd
from pcdl import TimeSeries

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

# Define the directory containing the simulation output files
output_dir = "/home/oth/BSC/physicell_mutations/mutPhysiBoSS/output"

# Create results directory structure
script_dir = os.path.dirname(os.path.abspath(__file__))
results_dir = os.path.join(script_dir, "results")
cell_cycle_dir = os.path.join(results_dir, "cell_cycle_phases")

# Create directories if they don't exist
os.makedirs(cell_cycle_dir, exist_ok=True)
print(f"Results will be saved to: {cell_cycle_dir}")

def analyze_cell_cycle_phases_pcdl_v4(output_dir):
    """Analyze cell cycle phase distributions using pcdl v4.0.5"""
    
    try:
        # Load the simulation data
        print("Loading simulation data with pcdl v4.0.5...")
        ts = TimeSeries(output_dir)
        
        # Get available time points using the correct v4.0.5 API
        time_points = []
        try:
            # Method 1: Check if timesteps attribute exists
            if hasattr(ts, 'timesteps'):
                time_points = ts.timesteps
                print(f"Found {len(time_points)} time steps using timesteps attribute")
            # Method 2: Check if get_times method exists
            elif hasattr(ts, 'get_times'):
                time_points = ts.get_times()
                print(f"Found {len(time_points)} time steps using get_times method")
            # Method 3: Try to get from l_mcds (list of mcds objects)
            elif hasattr(ts, 'l_mcds') and ts.l_mcds:
                time_points = [mcds.time for mcds in ts.l_mcds]
                print(f"Found {len(time_points)} time steps using l_mcds")
            else:
                print("Could not determine time points from TimeSeries object")
                return None, None
        except Exception as e:
            print(f"Error getting time points: {e}")
            return None, None
        
        # Initialize data structures
        phase_counts = {}
        valid_time_points = []
        
        # Process each time point
        for timestep in time_points:
            try:
                print(f"Processing time step {timestep}...")
                
                # Get cell data using the correct v4.0.5 API
                cell_df = ts.get_cell_df(timestep)
                
                if cell_df.empty:
                    print(f"  No cell data found at time step {timestep}")
                    continue
                
                print(f"  Found {len(cell_df)} cells")
                print(f"  Available columns: {list(cell_df.columns)}")
                
                # Check for cell cycle phase information
                phase_column = None
                for col in ['current_phase', 'phase', 'cell_phase', 'cycle_phase', 'cycle_model']:
                    if col in cell_df.columns:
                        phase_column = col
                        break
                
                if phase_column:
                    phase_distribution = cell_df[phase_column].value_counts().to_dict()
                    print(f"  Found phase data in column '{phase_column}': {phase_distribution}")
                else:
                    # If no phase column, use total cell count
                    total_cells = len(cell_df)
                    phase_distribution = {'total_cells': total_cells}
                    print(f"  No phase data found, using total cell count: {total_cells}")
                    
                    # Check for other relevant columns
                    relevant_cols = [col for col in cell_df.columns if any(keyword in col.lower() 
                                    for keyword in ['cycle', 'phase', 'state', 'model', 'type'])]
                    if relevant_cols:
                        print(f"  Found potentially relevant columns: {relevant_cols}")
                        for col in relevant_cols:
                            if cell_df[col].dtype == 'object' or cell_df[col].dtype.name == 'category':
                                col_dist = cell_df[col].value_counts().to_dict()
                                phase_distribution[col] = col_dist
                                print(f"  Column '{col}' distribution: {col_dist}")
                
                # Store the data
                phase_counts[timestep] = phase_distribution
                valid_time_points.append(timestep)
                
            except Exception as e:
                print(f"Error processing time step {timestep}: {e}")
                continue
        
        if not phase_counts:
            print("No valid cell cycle data found!")
            return None, None
        
        print(f"Successfully processed {len(valid_time_points)} time points")
        return phase_counts, valid_time_points
        
    except Exception as e:
        print(f"Error loading simulation data: {e}")
        return None, None

def plot_phase_distribution(phase_counts, time_points, output_path):
    """Plot cell cycle phase distribution over time"""
    
    if not phase_counts or not time_points:
        print("No data available for plotting")
        return
    
    # Get all unique phases across all time points
    all_phases = set()
    for phase_dict in phase_counts.values():
        all_phases.update(phase_dict.keys())
    all_phases = sorted(list(all_phases))
    
    print(f"Found data categories: {all_phases}")
    
    # Convert the phase counts to a structured format for plotting
    phase_data = {phase: [] for phase in all_phases}
    
    for time in time_points:
        for phase in all_phases:
            count = phase_counts[time].get(phase, 0)
            phase_data[phase].append(count)
    
    # Create the plot
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Define colors for different phases
    colors = plt.cm.tab10(np.linspace(0, 1, len(all_phases)))
    
    # Plot each phase/category
    for i, phase in enumerate(all_phases):
        if phase == 'total_cells':
            ax.plot(time_points, phase_data[phase], 
                   label='Total Cells', 
                   color=colors[i], 
                   linewidth=3, 
                   marker='o', 
                   markersize=6)
        else:
            ax.plot(time_points, phase_data[phase], 
                   label=f'{phase}', 
                   color=colors[i], 
                   linewidth=2, 
                   marker='o', 
                   markersize=4)
    
    # Formatting
    ax.set_xlabel('Time (hours)', fontsize=12)
    ax.set_ylabel('Number of Cells', fontsize=12)
    if 'total_cells' in all_phases and len(all_phases) == 1:
        ax.set_title('Total Cell Population Over Time', fontsize=14)
    else:
        ax.set_title('Cell Cycle Phase Distribution Over Time', fontsize=14)
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f"{output_path}_phase_distribution.png", bbox_inches='tight', dpi=300)
    plt.savefig(f"{output_path}_phase_distribution.pdf", bbox_inches='tight')
    plt.close()
    
    print(f"Phase distribution plot saved to: {output_path}_phase_distribution.png/pdf")

def plot_phase_percentages(phase_counts, time_points, output_path):
    """Plot cell cycle phase percentages over time (stacked area plot)"""
    
    if not phase_counts or not time_points:
        print("No data available for plotting")
        return
    
    # Get all unique phases across all time points
    all_phases = set()
    for phase_dict in phase_counts.values():
        all_phases.update(phase_dict.keys())
    all_phases = sorted(list(all_phases))
    
    # Convert to percentages
    phase_percentages = {phase: [] for phase in all_phases}
    
    for time in time_points:
        total_cells = sum(phase_counts[time].values())
        if total_cells > 0:
            for phase in all_phases:
                count = phase_counts[time].get(phase, 0)
                percentage = (count / total_cells) * 100
                phase_percentages[phase].append(percentage)
        else:
            for phase in all_phases:
                phase_percentages[phase].append(0)
    
    # Create the plot
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Define colors for different phases
    colors = plt.cm.tab10(np.linspace(0, 1, len(all_phases)))
    
    # Create stacked area plot
    ax.stackplot(time_points, 
                [phase_percentages[phase] for phase in all_phases],
                labels=[f'{phase}' for phase in all_phases],
                colors=colors,
                alpha=0.7)
    
    # Formatting
    ax.set_xlabel('Time (hours)', fontsize=12)
    ax.set_ylabel('Percentage of Cells (%)', fontsize=12)
    ax.set_title('Cell Cycle Phase Distribution Over Time (Percentages)', fontsize=14)
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    ax.set_ylim(0, 100)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f"{output_path}_phase_percentages.png", bbox_inches='tight', dpi=300)
    plt.savefig(f"{output_path}_phase_percentages.pdf", bbox_inches='tight')
    plt.close()
    
    print(f"Phase percentages plot saved to: {output_path}_phase_percentages.png/pdf")

def create_summary_table(phase_counts, time_points, output_path):
    """Create a summary table of cell cycle phase data"""
    
    if not phase_counts or not time_points:
        print("No data available for summary table")
        return
    
    # Get all unique phases
    all_phases = set()
    for phase_dict in phase_counts.values():
        all_phases.update(phase_dict.keys())
    all_phases = sorted(list(all_phases))
    
    # Create DataFrame
    data = []
    for time in time_points:
        row = {'Time': time}
        total_cells = sum(phase_counts[time].values())
        row['Total_Cells'] = total_cells
        
        for phase in all_phases:
            count = phase_counts[time].get(phase, 0)
            percentage = (count / total_cells * 100) if total_cells > 0 else 0
            row[f'{phase}_Count'] = count
            row[f'{phase}_Percentage'] = round(percentage, 2)
        
        data.append(row)
    
    df = pd.DataFrame(data)
    
    # Save to CSV
    csv_path = f"{output_path}_phase_summary.csv"
    df.to_csv(csv_path, index=False)
    print(f"Summary table saved to: {csv_path}")
    
    # Print summary statistics
    print("\n=== Cell Cycle Phase Summary ===")
    print(f"Time range: {min(time_points):.1f} - {max(time_points):.1f} hours")
    print(f"Total time points: {len(time_points)}")
    print(f"Phases observed: {all_phases}")
    print(f"Peak cell count: {max(df['Total_Cells'])} cells")
    
    return df

# Main execution
if __name__ == "__main__":
    print("=== Cell Cycle Phase Analysis (pcdl v4.0.5) ===")
    
    # Analyze cell cycle phases
    phase_counts, time_points = analyze_cell_cycle_phases_pcdl_v4(output_dir)
    
    if phase_counts and time_points:
        # Set output path
        base_path = os.path.join(cell_cycle_dir, "cell_cycle")
        
        # Create different types of plots
        plot_phase_distribution(phase_counts, time_points, base_path)
        plot_phase_percentages(phase_counts, time_points, base_path)
        
        # Create summary table
        summary_df = create_summary_table(phase_counts, time_points, base_path)
        
        print(f"\nAll cell cycle analysis plots saved to: {cell_cycle_dir}")
    else:
        print("Failed to load cell cycle data. Please check:")
        print("1. pcdl package is installed: pip install pcdl")
        print("2. Output directory contains valid PhysiCell data")
        print("3. Cell data contains phase information")

