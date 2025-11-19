import os
import re
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

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

script_dir = os.path.dirname(os.path.abspath(__file__))
repo_root = os.path.abspath(os.path.join(script_dir, os.pardir, os.pardir, os.pardir, os.pardir))
output_dir = os.path.join(repo_root, "output")

# Create results directory structure
results_dir = os.path.join(script_dir, "results")
aggregated_states_dir = os.path.join(results_dir, "aggregated_boolean_states")

# Create directories if they don't exist
os.makedirs(aggregated_states_dir, exist_ok=True)
print(f"Results will be saved to: {aggregated_states_dir}")


# Extract time index from filename
def extract_time(filename):
    match = re.search(r'output(\d+)_boolean_intracellular\.csv', filename)
    return int(match.group(1)) if match else -1


# Sanitize state labels
def clean_state(state_value):
    if pd.isna(state_value):
        return "No_state"
    state_value = str(state_value).strip()
    if not state_value or state_value == "<nil>":
        return "No_state"
    return state_value


# Read all boolean intracellular state files
state_files = [
    f for f in os.listdir(output_dir)
    if f.endswith("boolean_intracellular.csv")
]
state_files = sorted(state_files, key=extract_time)

records = []
for state_file in state_files:
    time_idx = extract_time(state_file)
    file_path = os.path.join(output_dir, state_file)

    try:
        df = pd.read_csv(file_path)
    except pd.errors.EmptyDataError:
        print(f"Warning: {state_file} is empty.")
        continue

    if "state" not in df.columns:
        print(f"Warning: {state_file} missing 'state' column. Skipping.")
        continue

    df["state"] = df["state"].apply(clean_state)

    print(f"File: {state_file}")
    print(df["state"].value_counts())

    for state, count in df["state"].value_counts().items():
        records.append({"time": time_idx, "state": state, "count": count})


# Combine into DataFrame
state_df = pd.DataFrame(records)
print(state_df.head())

if state_df.empty:
    raise SystemExit("No boolean state data found. Check your boolean intracellular CSV files.")

# Pivot: time x state
pivot_df = state_df.pivot_table(index="time", columns="state", values="count", fill_value=0)
pivot_df = pivot_df.sort_index()

# Sort by total abundance
state_totals = pivot_df.sum().sort_values(ascending=False)
top_n = 8
majority_states = list(state_totals.head(top_n).index)

# Plot absolute counts
fig, ax = plt.subplots(figsize=(10, 5))
pivot_df.plot.area(ax=ax, colormap="tab20")
ax.set_xlabel("Time Step", fontsize=12)
ax.set_ylabel("Number of Cells", fontsize=12)
ax.set_title("Boolean State Abundance Over Time", fontsize=14)

# Legend
handles, labels = ax.get_legend_handles_labels()
label_to_handle = dict(zip(labels, handles))
handles_major = [label_to_handle[l] for l in majority_states if l in label_to_handle]
labels_major = [l for l in majority_states if l in label_to_handle]

ax.legend(handles_major, labels_major, title="State (top)", bbox_to_anchor=(1.02, 1),
          loc="upper left", borderaxespad=0, fontsize=8, title_fontsize=9)

# Save paths
base_path = os.path.join(aggregated_states_dir, "boolean_states")

# Full legend
patch_handles = []
for h, lab in zip(handles, labels):
    if hasattr(h, "get_facecolor"):
        fc = h.get_facecolor()
        color = fc[0] if len(fc) else fc
    elif hasattr(h, "get_color"):
        color = h.get_color()
    else:
        color = "gray"
    patch_handles.append(mpl.patches.Patch(facecolor=color, label=lab))

fig_leg = plt.figure(figsize=(3, max(1, 0.3 * len(patch_handles))))
fig_leg.legend(handles=patch_handles, labels=[ph.get_label() for ph in patch_handles],
               ncol=1, frameon=False, loc="center", title="State")
fig_leg.tight_layout()
fig_leg.savefig(f"{base_path}_legend.png", bbox_inches="tight")
fig_leg.savefig(f"{base_path}_legend.pdf", bbox_inches="tight")
plt.close(fig_leg)

# Save main figure
fig.tight_layout()
fig.savefig(f"{base_path}_aggregated.png", bbox_inches="tight")
fig.savefig(f"{base_path}_aggregated.pdf", bbox_inches="tight")

print(f"Figures saved to: {base_path}_aggregated.png and {base_path}_aggregated.pdf")

# --------------------------------------------------
# Normalized (100%) stacked AREA plot of state proportions
# --------------------------------------------------
# Compute percentage per time step
pivot_pct = pivot_df.div(pivot_df.sum(axis=1), axis=0).fillna(0) * 100

# Consistent colors: map each state to its area color
state_colors = {}
for lbl, hndl in label_to_handle.items():
    if hasattr(hndl, "get_facecolor"):
        fc = hndl.get_facecolor()
        state_colors[lbl] = fc[0] if len(fc) else fc
    elif hasattr(hndl, "get_color"):
        state_colors[lbl] = hndl.get_color()

# Define color list in column order for area plot
color_list = [state_colors.get(col, None) for col in pivot_pct.columns]

fig_prop, ax_prop = plt.subplots(figsize=(10, 5))
pivot_pct.plot.area(ax=ax_prop, color=color_list)

# Formatting – journal style
ax_prop.set_xlabel("Time Step", fontsize=12)
ax_prop.set_ylabel("State Proportion (%)", fontsize=12)
ax_prop.set_title("Normalized Boolean State Proportions Over Time", fontsize=14, pad=12)
ax_prop.set_ylim(0, 100)

# Clean up spines for a minimalist look
for spine in ["top", "right", "left", "bottom"]:
    ax_prop.spines[spine].set_visible(False)
ax_prop.grid(axis="y", color="0.9", linewidth=0.7)

# Legend outside – top states only (fresh Patch handles)
prop_patch_handles = [
    mpl.patches.Patch(facecolor=state_colors.get(s, "gray"), label=s)
    for s in majority_states if s in state_colors
]
ax_prop.legend(handles=prop_patch_handles,
               labels=[p.get_label() for p in prop_patch_handles],
               title="State (top)", bbox_to_anchor=(1.02, 1), loc="upper left",
               borderaxespad=0, fontsize=8, title_fontsize=9)

fig_prop.tight_layout()
fig_prop.savefig(f"{base_path}_aggregated_proportion.png", bbox_inches="tight")
fig_prop.savefig(f"{base_path}_aggregated_proportion.pdf", bbox_inches="tight")
plt.close(fig_prop)

# Build fishplot-compatible DataFrames
state_ids = {state: i for i, state in enumerate(pivot_df.columns)}
populations = []
for time, row in pivot_df.iterrows():
    for state, count in row.items():
        populations.append({"Id": state_ids[state], "Step": time, "Pop": count})
pop_df = pd.DataFrame(populations)

# Set all parent IDs to 0 for now
parent_df = pd.DataFrame({"ParentId": [0] * len(state_ids), "ChildId": list(state_ids.values())})

print("Fishplot dataframes (pop_df, parent_df) prepared for downstream use.")


