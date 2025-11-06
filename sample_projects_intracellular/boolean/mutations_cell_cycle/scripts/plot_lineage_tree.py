# this script will read the .mut files from the output directory and plot lineage trees
# showing parent-daughter relationships

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
lineage_dir = os.path.join(results_dir, "lineage_trees")

# Create directories if they don't exist
os.makedirs(lineage_dir, exist_ok=True)
print(f"Results will be saved to: {lineage_dir}")

# Extract time index from filename
def extract_time(filename):
    match = re.search(r'output(\d+)_mut', filename)
    return int(match.group(1)) if match else -1

# Transform full mutation string into final genotype
def final_genotype(mutation_str):
    if not mutation_str or pd.isna(mutation_str):
        return 'No_mutations'
    mutations = mutation_str.split('.')
    gene_state = {}
    for mut in mutations:
        if '_' in mut:
            gene, state = mut.split('_')
            gene_state[gene] = state  # keep latest mutation per gene
    # Sort genes alphabetically for consistency
    final = [f"{gene}_{gene_state[gene]}" for gene in sorted(gene_state.keys())]
    return '.'.join(final)

# Build lineage tree from mutation data
def build_lineage_tree(mut_files, output_dir):
    """Build a lineage tree from mutation files showing parent-daughter relationships"""
    lineage_data = []
    
    for mut_file in mut_files:
        time_idx = extract_time(mut_file)
        if time_idx == -1:  # Skip non-time-indexed files
            continue
            
        df = pd.read_csv(os.path.join(output_dir, mut_file))
        df['mutations'] = df['mutations'].fillna('')
        df['final_lineage'] = df['mutations'].apply(final_genotype)
        
        for _, row in df.iterrows():
            lineage_data.append({
                'time': time_idx,
                'cell_id': row['ID'],
                'parent_id': row['parent_ID'],
                'generation': row['generation'],
                'mutations': row['mutations'],
                'final_lineage': row['final_lineage']
            })
    
    return pd.DataFrame(lineage_data)

# Plot lineage tree
def plot_lineage_tree(lineage_df, output_path):
    """Create a lineage tree plot showing parent-daughter relationships"""
    if lineage_df.empty:
        print("No lineage data available for tree plot")
        return
    
    # Create figure
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Get unique lineages and assign colors
    unique_lineages = lineage_df['final_lineage'].unique()
    colors = plt.cm.tab20(np.linspace(0, 1, len(unique_lineages)))
    lineage_colors = dict(zip(unique_lineages, colors))
    
    # Plot each cell as a point
    for lineage in unique_lineages:
        lineage_data = lineage_df[lineage_df['final_lineage'] == lineage]
        ax.scatter(lineage_data['time'], lineage_data['cell_id'], 
                  c=[lineage_colors[lineage]], label=lineage, s=50, alpha=0.7)
    
    # Draw parent-daughter connections
    for _, row in lineage_df.iterrows():
        if row['parent_id'] != -1:  # Not a root cell
            parent_data = lineage_df[
                (lineage_df['cell_id'] == row['parent_id']) & 
                (lineage_df['time'] < row['time'])
            ]
            if not parent_data.empty:
                parent_time = parent_data.iloc[-1]['time']  # Get latest parent time
                ax.plot([parent_time, row['time']], 
                       [row['parent_id'], row['cell_id']], 
                       'k-', alpha=0.3, linewidth=0.5)
    
    # Formatting
    ax.set_xlabel('Time Step', fontsize=12)
    ax.set_ylabel('Cell ID', fontsize=12)
    ax.set_title('Cell Lineage Tree', fontsize=14)
    ax.grid(True, alpha=0.3)
    
    # Legend
    handles = [plt.Line2D([0], [0], marker='o', color='w', 
                         markerfacecolor=lineage_colors[lineage], 
                         markersize=8, label=lineage) 
               for lineage in unique_lineages]
    ax.legend(handles=handles, title='Lineage', bbox_to_anchor=(1.05, 1), 
              loc='upper left', fontsize=8)
    
    plt.tight_layout()
    plt.savefig(f"{output_path}_lineage_tree.png", bbox_inches='tight', dpi=300)
    plt.savefig(f"{output_path}_lineage_tree.pdf", bbox_inches='tight')
    plt.close()
    
    print(f"Lineage tree plot saved to: {output_path}_lineage_tree.png/pdf")

# Plot lineage tree with generation levels
def plot_lineage_generations(lineage_df, output_path):
    """Create a lineage plot organized by generations"""
    if lineage_df.empty:
        print("No lineage data available for generation plot")
        return
    
    # Create figure
    fig, ax = plt.subplots(figsize=(14, 8))
    
    # Get unique lineages and assign colors
    unique_lineages = lineage_df['final_lineage'].unique()
    colors = plt.cm.tab20(np.linspace(0, 1, len(unique_lineages)))
    lineage_colors = dict(zip(unique_lineages, colors))
    
    # Plot cells organized by generation
    max_generation = lineage_df['generation'].max()
    
    for generation in range(max_generation + 1):
        gen_data = lineage_df[lineage_df['generation'] == generation]
        if gen_data.empty:
            continue
            
        # Spread cells horizontally within generation
        y_pos = generation
        x_positions = np.linspace(0, 1, len(gen_data))
        
        for i, (_, row) in enumerate(gen_data.iterrows()):
            ax.scatter(row['time'], y_pos + x_positions[i] * 0.1, 
                      c=[lineage_colors[row['final_lineage']]], 
                      s=50, alpha=0.7, edgecolors='black', linewidth=0.5)
            
            # Add cell ID as text
            ax.annotate(f"ID:{row['cell_id']}", 
                       (row['time'], y_pos + x_positions[i] * 0.1),
                       xytext=(5, 5), textcoords='offset points',
                       fontsize=6, alpha=0.8)
    
    # Draw parent-daughter connections
    for _, row in lineage_df.iterrows():
        if row['parent_id'] != -1:  # Not a root cell
            parent_data = lineage_df[
                (lineage_df['cell_id'] == row['parent_id']) & 
                (lineage_df['time'] < row['time'])
            ]
            if not parent_data.empty:
                parent_time = parent_data.iloc[-1]['time']
                parent_gen = parent_data.iloc[-1]['generation']
                
                # Calculate positions
                parent_y = parent_gen
                daughter_y = row['generation']
                
                ax.plot([parent_time, row['time']], 
                       [parent_y, daughter_y], 
                       'k-', alpha=0.4, linewidth=1)
    
    # Formatting
    ax.set_xlabel('Time Step', fontsize=12)
    ax.set_ylabel('Generation', fontsize=12)
    ax.set_title('Cell Lineage by Generation', fontsize=14)
    ax.set_yticks(range(max_generation + 1))
    ax.grid(True, alpha=0.3)
    
    # Legend
    handles = [plt.Line2D([0], [0], marker='o', color='w', 
                         markerfacecolor=lineage_colors[lineage], 
                         markersize=8, label=lineage) 
               for lineage in unique_lineages]
    ax.legend(handles=handles, title='Lineage', bbox_to_anchor=(1.05, 1), 
              loc='upper left', fontsize=8)
    
    plt.tight_layout()
    plt.savefig(f"{output_path}_lineage_generations.png", bbox_inches='tight', dpi=300)
    plt.savefig(f"{output_path}_lineage_generations.pdf", bbox_inches='tight')
    plt.close()
    
    print(f"Lineage generations plot saved to: {output_path}_lineage_generations.png/pdf")

# Plot lineage tree as a proper tree structure
def plot_lineage_tree_structure(lineage_df, output_path):
    """Create a proper tree structure plot showing lineage hierarchy"""
    if lineage_df.empty:
        print("No lineage data available for tree structure plot")
        return
    
    # Build tree structure
    tree = defaultdict(list)
    root_cells = []
    
    for _, row in lineage_df.iterrows():
        if row['parent_id'] == -1:
            root_cells.append(row['cell_id'])
        else:
            tree[row['parent_id']].append(row['cell_id'])
    
    # Create figure
    fig, ax = plt.subplots(figsize=(16, 10))
    
    # Get unique lineages and assign colors
    unique_lineages = lineage_df['final_lineage'].unique()
    colors = plt.cm.tab20(np.linspace(0, 1, len(unique_lineages)))
    lineage_colors = dict(zip(unique_lineages, colors))
    
    # Position nodes
    node_positions = {}
    y_levels = {}
    
    def position_nodes(cell_id, level=0, x_offset=0):
        if cell_id in node_positions:
            return x_offset
        
        # Get cell data
        cell_data = lineage_df[lineage_df['cell_id'] == cell_id]
        if cell_data.empty:
            return x_offset
        
        cell_data = cell_data.iloc[0]
        
        # Position this node
        node_positions[cell_id] = (x_offset, level)
        y_levels[level] = y_levels.get(level, [])
        y_levels[level].append(x_offset)
        
        # Position children
        children = tree.get(cell_id, [])
        if children:
            child_x = x_offset
            for child in children:
                child_x = position_nodes(child, level + 1, child_x) + 1
        else:
            return x_offset + 1
        
        return child_x
    
    # Position all root nodes
    x_offset = 0
    for root in root_cells:
        x_offset = position_nodes(root, 0, x_offset) + 2
    
    # Plot nodes
    for cell_id, (x, y) in node_positions.items():
        cell_data = lineage_df[lineage_df['cell_id'] == cell_id].iloc[0]
        lineage = cell_data['final_lineage']
        
        ax.scatter(x, y, c=[lineage_colors[lineage]], s=100, 
                  alpha=0.8, edgecolors='black', linewidth=1)
        
        # Add cell ID and generation
        ax.annotate(f"ID:{cell_id}\nGen:{cell_data['generation']}", 
                   (x, y), xytext=(0, 10), textcoords='offset points',
                   ha='center', fontsize=8, alpha=0.9)
    
    # Draw edges
    for parent_id, children in tree.items():
        if parent_id in node_positions:
            parent_x, parent_y = node_positions[parent_id]
            for child_id in children:
                if child_id in node_positions:
                    child_x, child_y = node_positions[child_id]
                    ax.plot([parent_x, child_x], [parent_y, child_y], 
                           'k-', alpha=0.6, linewidth=1)
    
    # Formatting
    ax.set_xlabel('Lineage Branch', fontsize=12)
    ax.set_ylabel('Generation Level', fontsize=12)
    ax.set_title('Cell Lineage Tree Structure', fontsize=14)
    ax.set_yticks(range(max(y_levels.keys()) + 1))
    ax.grid(True, alpha=0.3)
    
    # Legend
    handles = [plt.Line2D([0], [0], marker='o', color='w', 
                         markerfacecolor=lineage_colors[lineage], 
                         markersize=10, label=lineage) 
               for lineage in unique_lineages]
    ax.legend(handles=handles, title='Lineage', bbox_to_anchor=(1.05, 1), 
              loc='upper left', fontsize=8)
    
    plt.tight_layout()
    plt.savefig(f"{output_path}_lineage_structure.png", bbox_inches='tight', dpi=300)
    plt.savefig(f"{output_path}_lineage_structure.pdf", bbox_inches='tight')
    plt.close()
    
    print(f"Lineage structure plot saved to: {output_path}_lineage_structure.png/pdf")

# Main execution
if __name__ == "__main__":
    # Read all .mut.csv files
    mut_files = [f for f in os.listdir(output_dir) if f.endswith('mut.csv')]
    mut_files = sorted(mut_files, key=extract_time)
    
    print(f"Found {len(mut_files)} mutation files")
    
    # Build lineage tree
    print("Building lineage tree...")
    lineage_df = build_lineage_tree(mut_files, output_dir)
    
    if not lineage_df.empty:
        print(f"Lineage data: {len(lineage_df)} cell records")
        print(f"Time range: {lineage_df['time'].min()} to {lineage_df['time'].max()}")
        print(f"Generations: {lineage_df['generation'].min()} to {lineage_df['generation'].max()}")
        print(f"Unique lineages: {lineage_df['final_lineage'].nunique()}")
        
        # Set output path
        base_path = os.path.join(lineage_dir, "lineage")
        
        # Create different types of lineage plots
        plot_lineage_tree(lineage_df, base_path)
        plot_lineage_generations(lineage_df, base_path)
        plot_lineage_tree_structure(lineage_df, base_path)
        
        print(f"\nAll lineage plots saved to: {lineage_dir}")
    else:
        print("No lineage data found for plotting")

