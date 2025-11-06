# this script will use pctk (PhysiCell ToolKit) to analyze cell cycle phases
# pctk is specifically designed for PhysiCell/PhysiBoSS simulations

import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pandas as pd
from pctk import multicellds

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

def analyze_cell_cycle_phases_pctk(output_dir):
    """Analyze cell cycle phase distributions using pctk"""
    
    try:
        # Create a MultiCellDS reader
        print("Loading simulation data with pctk...")
        reader = multicellds.MultiCellDS(output_folder=output_dir)
        
        # Create an iterator to load cell DataFrames for each stored simulation time step
        print("Creating cell data iterator...")
        try:
            df_iterator = reader.cells_as_frames_iterator()
        except Exception as e:
            print(f"Error creating iterator: {e}")
            print("\nNOTE: This error typically occurs when .cells.mat files are corrupted.")
            print("This is a known issue with PhysiCell MultiCellDS version 0.5.")
            print("The cell cycle phase data cannot be extracted from corrupted files.")
            return None, None
        
        # Initialize data structures
        phase_counts = {}
        valid_time_points = []
        error_count = 0
        
        # Iterate over all simulation outputs
        for t, df_cells in df_iterator:
            try:
                print(f"Processing time step {t}...")
                
                if df_cells.empty:
                    print(f"  No cell data found at time step {t}")
                    continue
                
                print(f"  Found {len(df_cells)} cells")
                print(f"  Available columns: {list(df_cells.columns)}")
                
                # Check for cell cycle phase information
                phase_column = None
                for col in ['current_phase', 'phase', 'cell_phase', 'cycle_phase', 'cycle_model']:
                    if col in df_cells.columns:
                        phase_column = col
                        break
                
                if phase_column:
                    phase_distribution = df_cells[phase_column].value_counts().to_dict()
                    print(f"  Found phase data in column '{phase_column}': {phase_distribution}")
                else:
                    # If no phase column, use total cell count
                    total_cells = len(df_cells)
                    phase_distribution = {'total_cells': total_cells}
                    print(f"  No phase data found, using total cell count: {total_cells}")
                    
                    # Check for other relevant columns
                    relevant_cols = [col for col in df_cells.columns if any(keyword in col.lower() 
                                    for keyword in ['cycle', 'phase', 'state', 'model', 'type'])]
                    if relevant_cols:
                        print(f"  Found potentially relevant columns: {relevant_cols}")
                        for col in relevant_cols:
                            if df_cells[col].dtype == 'object' or df_cells[col].dtype.name == 'category':
                                col_dist = df_cells[col].value_counts().to_dict()
                                phase_distribution[col] = col_dist
                                print(f"  Column '{col}' distribution: {col_dist}")
                
                # Store the data
                phase_counts[t] = phase_distribution
                valid_time_points.append(t)
                
            except Exception as e:
                print(f"Error processing time step {t}: {e}")
                continue
        
        if not phase_counts:
            print("No valid cell cycle data found!")
            return None, None
        
        print(f"Successfully processed {len(valid_time_points)} time points")
        return phase_counts, valid_time_points
        
    except Exception as e:
        print(f"Error loading simulation data: {e}")
        import traceback
        traceback.print_exc()
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

def plot_alive_dead_cells(phase_counts, time_points, output_path):
    """Plot alive vs dead cells over time (if current_phase is available)"""
    
    if not phase_counts or not time_points:
        print("No data available for plotting")
        return
    
    # Check if we have current_phase data
    has_phase_data = any('current_phase' in str(phase) for phase in phase_counts.values())
    
    if not has_phase_data:
        # Try to extract from the phase_counts structure
        # In pctk, alive cells typically have current_phase <= 14
        alive_counts = []
        dead_counts = []
        
        for time in time_points:
            # Try to get phase information
            phase_dict = phase_counts[time]
            alive = 0
            dead = 0
            
            # Check if we have numeric phase values
            for phase_key, count in phase_dict.items():
                if isinstance(phase_key, (int, float)):
                    if phase_key <= 14:
                        alive += count
                    else:
                        dead += count
                elif 'total_cells' in str(phase_key):
                    # If we only have total cells, we can't distinguish alive/dead
                    pass
            
            alive_counts.append(alive)
            dead_counts.append(dead)
        
        if sum(alive_counts) == 0 and sum(dead_counts) == 0:
            print("Cannot determine alive/dead cells from available data")
            return
    
    # Create the plot
    fig, ax = plt.subplots(figsize=(12, 8))
    
    ax.plot(time_points, alive_counts, 'g-', label='Alive Cells', linewidth=2, marker='o', markersize=4)
    ax.plot(time_points, dead_counts, 'r-', label='Dead Cells', linewidth=2, marker='s', markersize=4)
    
    # Formatting
    ax.set_xlabel('Time (hours)', fontsize=12)
    ax.set_ylabel('Number of Cells', fontsize=12)
    ax.set_title('Alive vs Dead Cells Over Time', fontsize=14)
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f"{output_path}_alive_dead.png", bbox_inches='tight', dpi=300)
    plt.savefig(f"{output_path}_alive_dead.pdf", bbox_inches='tight')
    plt.close()
    
    print(f"Alive/dead cells plot saved to: {output_path}_alive_dead.png/pdf")

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
    if 'Total_Cells' in df.columns:
        print(f"Peak cell count: {max(df['Total_Cells'])} cells")
    
    return df

# Main execution
if __name__ == "__main__":
    print("=== Cell Cycle Phase Analysis (pctk) ===")
    
    # Analyze cell cycle phases using pctk
    phase_counts, time_points = analyze_cell_cycle_phases_pctk(output_dir)
    
    if phase_counts and time_points:
        # Set output path
        base_path = os.path.join(cell_cycle_dir, "cell_cycle")
        
        # Create different types of plots
        plot_phase_distribution(phase_counts, time_points, base_path)
        plot_phase_percentages(phase_counts, time_points, base_path)
        plot_alive_dead_cells(phase_counts, time_points, base_path)
        
        # Create summary table
        summary_df = create_summary_table(phase_counts, time_points, base_path)
        
        print(f"\nAll cell cycle analysis plots saved to: {cell_cycle_dir}")
    else:
        print("\n" + "="*60)
        print("FAILED TO LOAD CELL CYCLE DATA")
        print("="*60)
        print("\nPossible reasons:")
        print("1. Corrupted .cells.mat files (known PhysiCell MultiCellDS v0.5 bug)")
        print("2. Missing cell cycle phase information in output")
        print("3. Incompatible PhysiCell version")
        print("\nRecommendations:")
        print("- Check PhysiCell_settings.xml for cell cycle model configuration")
        print("- Verify that cell cycle data is being written to output files")
        print("- Consider using plot_cell_population.py for population dynamics analysis")
        print("- The .mut.csv files work correctly and can be analyzed separately")
        print("="*60)
