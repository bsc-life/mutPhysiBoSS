# this script will read the .mut files and plot cell population dynamics over time
# This is a simpler approach that works with the available data

import os
import re
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

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
population_dir = os.path.join(results_dir, "cell_population")

# Create directories if they don't exist
os.makedirs(population_dir, exist_ok=True)
print(f"Results will be saved to: {population_dir}")

# Extract time index from filename
def extract_time(filename):
    match = re.search(r'output(\d+)_mut', filename)
    return int(match.group(1)) if match else -1

def analyze_cell_population(output_dir):
    """Analyze cell population dynamics from .mut.csv files"""
    
    # Read all .mut.csv files
    mut_files = [f for f in os.listdir(output_dir) if f.endswith('mut.csv')]
    mut_files = sorted(mut_files, key=extract_time)
    
    print(f"Found {len(mut_files)} mutation files")
    
    # Initialize data structures
    population_data = []
    time_points = []
    
    for mut_file in mut_files:
        time_idx = extract_time(mut_file)
        if time_idx == -1:  # Skip non-time-indexed files
            continue
            
        try:
            df = pd.read_csv(os.path.join(output_dir, mut_file))
            
            # Count total cells
            total_cells = len(df)
            
            # Count cells by generation
            generation_counts = df['generation'].value_counts().to_dict()
            
            # Count cells by type
            type_counts = df['type'].value_counts().to_dict()
            
            population_data.append({
                'time': time_idx,
                'total_cells': total_cells,
                'generation_counts': generation_counts,
                'type_counts': type_counts
            })
            time_points.append(time_idx)
            
            print(f"Time {time_idx}: {total_cells} cells, generations: {generation_counts}")
            
        except Exception as e:
            print(f"Error processing {mut_file}: {e}")
            continue
    
    return population_data, time_points

def plot_total_population(population_data, time_points, output_path):
    """Plot total cell population over time"""
    
    if not population_data:
        print("No population data available for plotting")
        return
    
    # Extract total cell counts
    total_cells = [data['total_cells'] for data in population_data]
    
    # Create the plot
    fig, ax = plt.subplots(figsize=(12, 8))
    
    ax.plot(time_points, total_cells, 'b-', linewidth=3, marker='o', markersize=6, label='Total Cells')
    
    # Formatting
    ax.set_xlabel('Time Step', fontsize=12)
    ax.set_ylabel('Number of Cells', fontsize=12)
    ax.set_title('Total Cell Population Over Time', fontsize=14)
    ax.grid(True, alpha=0.3)
    ax.legend()
    
    plt.tight_layout()
    plt.savefig(f"{output_path}_total_population.png", bbox_inches='tight', dpi=300)
    plt.savefig(f"{output_path}_total_population.pdf", bbox_inches='tight')
    plt.close()
    
    print(f"Total population plot saved to: {output_path}_total_population.png/pdf")

def plot_generation_distribution(population_data, time_points, output_path):
    """Plot cell generation distribution over time"""
    
    if not population_data:
        print("No population data available for plotting")
        return
    
    # Get all unique generations
    all_generations = set()
    for data in population_data:
        all_generations.update(data['generation_counts'].keys())
    all_generations = sorted(list(all_generations))
    
    print(f"Found generations: {all_generations}")
    
    # Create the plot
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Define colors for different generations
    colors = plt.cm.tab10(np.linspace(0, 1, len(all_generations)))
    
    # Plot each generation
    for i, generation in enumerate(all_generations):
        gen_counts = []
        for data in population_data:
            gen_counts.append(data['generation_counts'].get(generation, 0))
        
        ax.plot(time_points, gen_counts, 
               label=f'Generation {generation}', 
               color=colors[i], 
               linewidth=2, 
               marker='o', 
               markersize=4)
    
    # Formatting
    ax.set_xlabel('Time Step', fontsize=12)
    ax.set_ylabel('Number of Cells', fontsize=12)
    ax.set_title('Cell Generation Distribution Over Time', fontsize=14)
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f"{output_path}_generation_distribution.png", bbox_inches='tight', dpi=300)
    plt.savefig(f"{output_path}_generation_distribution.pdf", bbox_inches='tight')
    plt.close()
    
    print(f"Generation distribution plot saved to: {output_path}_generation_distribution.png/pdf")

def plot_generation_percentages(population_data, time_points, output_path):
    """Plot cell generation percentages over time (stacked area plot)"""
    
    if not population_data:
        print("No population data available for plotting")
        return
    
    # Get all unique generations
    all_generations = set()
    for data in population_data:
        all_generations.update(data['generation_counts'].keys())
    all_generations = sorted(list(all_generations))
    
    # Convert to percentages
    generation_percentages = {generation: [] for generation in all_generations}
    
    for data in population_data:
        total_cells = data['total_cells']
        if total_cells > 0:
            for generation in all_generations:
                count = data['generation_counts'].get(generation, 0)
                percentage = (count / total_cells) * 100
                generation_percentages[generation].append(percentage)
        else:
            for generation in all_generations:
                generation_percentages[generation].append(0)
    
    # Create the plot
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Define colors for different generations
    colors = plt.cm.tab10(np.linspace(0, 1, len(all_generations)))
    
    # Create stacked area plot
    ax.stackplot(time_points, 
                [generation_percentages[generation] for generation in all_generations],
                labels=[f'Generation {generation}' for generation in all_generations],
                colors=colors,
                alpha=0.7)
    
    # Formatting
    ax.set_xlabel('Time Step', fontsize=12)
    ax.set_ylabel('Percentage of Cells (%)', fontsize=12)
    ax.set_title('Cell Generation Distribution Over Time (Percentages)', fontsize=14)
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    ax.set_ylim(0, 100)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f"{output_path}_generation_percentages.png", bbox_inches='tight', dpi=300)
    plt.savefig(f"{output_path}_generation_percentages.pdf", bbox_inches='tight')
    plt.close()
    
    print(f"Generation percentages plot saved to: {output_path}_generation_percentages.png/pdf")

def plot_division_rates(population_data, time_points, output_path):
    """Plot cell division rates over time"""
    
    if not population_data or len(population_data) < 2:
        print("Insufficient data for division rate analysis")
        return
    
    # Calculate total cell counts
    total_cells = [data['total_cells'] for data in population_data]
    
    # Calculate division rates (derivative of total cell count)
    division_rates = np.diff(total_cells)
    
    # Create the plot
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))
    
    # Plot 1: Total cell count over time
    ax1.plot(time_points, total_cells, 'b-', linewidth=3, marker='o', markersize=6)
    ax1.set_xlabel('Time Step', fontsize=12)
    ax1.set_ylabel('Total Number of Cells', fontsize=12)
    ax1.set_title('Total Cell Population Over Time', fontsize=14)
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Cell division rate
    ax2.plot(time_points[1:], division_rates, 'r-', linewidth=2, marker='s', markersize=4)
    ax2.set_xlabel('Time Step', fontsize=12)
    ax2.set_ylabel('Cell Division Rate (cells/step)', fontsize=12)
    ax2.set_title('Cell Division Rate Over Time', fontsize=14)
    ax2.grid(True, alpha=0.3)
    ax2.axhline(y=0, color='k', linestyle='--', alpha=0.5)
    
    plt.tight_layout()
    plt.savefig(f"{output_path}_division_rates.png", bbox_inches='tight', dpi=300)
    plt.savefig(f"{output_path}_division_rates.pdf", bbox_inches='tight')
    plt.close()
    
    print(f"Division rates plot saved to: {output_path}_division_rates.png/pdf")

def create_summary_table(population_data, time_points, output_path):
    """Create a summary table of cell population data"""
    
    if not population_data:
        print("No population data available for summary table")
        return
    
    # Get all unique generations
    all_generations = set()
    for data in population_data:
        all_generations.update(data['generation_counts'].keys())
    all_generations = sorted(list(all_generations))
    
    # Create DataFrame
    data_rows = []
    for i, data in enumerate(population_data):
        row = {'Time': time_points[i], 'Total_Cells': data['total_cells']}
        
        for generation in all_generations:
            count = data['generation_counts'].get(generation, 0)
            percentage = (count / data['total_cells'] * 100) if data['total_cells'] > 0 else 0
            row[f'Gen_{generation}_Count'] = count
            row[f'Gen_{generation}_Percentage'] = round(percentage, 2)
        
        data_rows.append(row)
    
    df = pd.DataFrame(data_rows)
    
    # Save to CSV
    csv_path = f"{output_path}_population_summary.csv"
    df.to_csv(csv_path, index=False)
    print(f"Summary table saved to: {csv_path}")
    
    # Print summary statistics
    print("\n=== Cell Population Summary ===")
    print(f"Time range: {min(time_points)} - {max(time_points)} steps")
    print(f"Total time points: {len(time_points)}")
    print(f"Generations observed: {all_generations}")
    print(f"Peak cell count: {max(df['Total_Cells'])} cells")
    print(f"Final cell count: {df['Total_Cells'].iloc[-1]} cells")
    
    return df

# Main execution
if __name__ == "__main__":
    print("=== Cell Population Analysis ===")
    
    # Analyze cell population
    population_data, time_points = analyze_cell_population(output_dir)
    
    if population_data and time_points:
        # Set output path
        base_path = os.path.join(population_dir, "population")
        
        # Create different types of plots
        plot_total_population(population_data, time_points, base_path)
        plot_generation_distribution(population_data, time_points, base_path)
        plot_generation_percentages(population_data, time_points, base_path)
        plot_division_rates(population_data, time_points, base_path)
        
        # Create summary table
        summary_df = create_summary_table(population_data, time_points, base_path)
        
        print(f"\nAll population analysis plots saved to: {population_dir}")
    else:
        print("Failed to load population data. Please check:")
        print("1. Output directory contains valid .mut.csv files")
        print("2. Files are properly formatted")

