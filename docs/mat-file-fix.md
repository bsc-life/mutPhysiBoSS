# MAT Export Fix — November 2025

## Summary
- **Issue**: `pcdl` flagged every `outputXXXX_cells.mat` snapshot as corrupt (`Not enough bytes to read matrix 'cells'`), so the cell-cycle analysis script could not load any timesteps.
- **Root Cause**: We added a `mutations` entry to the cells legend but never wrote a corresponding scalar to the MAT stream. The column count in the header no longer matched the number of doubles appended, truncating each file.
- **Fix**: Reintroduced a single `std::fwrite` with the per-cell mutation count inside `modules/PhysiCell_MultiCellDS.cpp` (`add_PhysiCell_cells_to_open_xml_pugi_v2`). This restores parity between the header declaration and the actual payload.
- **Outcome**: Regenerated snapshots load with `scipy.io.loadmat`, and `plot_cell_cycle_phases_final.py` now progresses past the MAT read (remaining warnings are only about missing `states_*.csv` PhysioBOSS exports).

## Diagnosis Steps
- Ran `plot_cell_cycle_phases_final.py` → repeated `corrupt ... cells.mat` warnings.
- Used `scipy.io.whosmat` to confirm the expected shape (104×1), but `loadmat` failed for insufficient bytes.
- Compared MAT header/size expectations to the stream written in `PhysiCell_MultiCellDS.cpp` and noticed the missing `fwrite` for the new label.

## Remediation
```diff
// custom vector variables
// ...
// mutations
- // (previously empty)
+ double mutation_count = static_cast<double>(pCell->custom_data.mutations.size());
+ std::fwrite(&(mutation_count), sizeof(double), 1, fp);
```
Full context committed in `modules/PhysiCell_MultiCellDS.cpp`.

## Verification
- `make` to rebuild `mutations_cell_cycle`.
- Reran the simulation to regenerate MAT snapshots.
- `python3 - <<'PY' ... loadmat('output00000000_cells.mat') ...` → succeeded, showing a 104×1 array.
- `plot_cell_cycle_phases_final.py` now ingests the MAT outputs without corruption errors.

## Follow-Up
- Consider exporting the PhysioBOSS `states_*.csv` files to silence the remaining warnings and expose per-time-step metadata such as `TimeStep.time`.
- Add a regression unit test that opens a newly written `cells.mat` with SciPy (or MATLAB) to catch header/body mismatches.


