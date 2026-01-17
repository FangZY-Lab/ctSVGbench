# ctSVGbench
Benchmarking existing computational approaches for detecting cell-type-specific spatially variable genes (ctSVGs) in spatial transcriptomics data.

## Project Structure

### Real Data Analysis
Scripts for analyzing and visualizing real datasets.

#### Core analysis
- **`real/fig4_rotate_sc.R`**: Rotate spatial coordinates and run ctSVG methods (single-cell data).
- **`real/fig4_rotate_spot.R`**: Rotate spatial coordinates and run ctSVG methods (spot-level data).
- **`real/real datasets1.R`**: Run real-data benchmarking (spot level data).
- **`real/real datasets2.R`**: Run real-data benchmarking (spot level data with spot subsetting).
- **`real/real datasets3sc.R`**: Run real-data benchmarking (single-cell level data).
- **`real/scalability.R`**: Scalability evaluation on real datasets (runtime/memory under varying spot sizes).

#### Preprocessing
- **`real/preprocess_1.2_public_data.R`**: Preprocess public datasets used in benchmarks.

- **`real/preprocess_2_bounary.R`**: Performs boundary preprocessing for real data. Required by spVC.
- **`real/preprocess_3_subset_spots.R`**: Preprocesses data for subset spots.

#### Plotting (figures)
- **`real/plot_fig2_realdata.R`**: Figure 2 — consistency analysis on real data.
- **`real/plot_fig4_rotate_sc.R`**: Figure 4 — rotation results (sc data).
- **`real/plot_fig4_rotate_spot.R`**: Figure 4 — rotation results (spot data).
- **`real/plot_fig5_time_mem.R`**: Figure 5 — runtime and memory benchmarking.
- **`real/plot_fig6_lung cancer.R`**: Figure 6 — lung cancer results.
- **`real/plot_fig7_mbm.R`**: Figure 7 — MBM dataset results.
- **`real/plot_fig8_summary.R`**: Figure 8 — summary of overall performance.
- **`real/plot_figs1_datasets.R`**: Supplementary Fig S1 — dataset summary used in Fig S1.
- **`real/plot_lung_stance_conserv.R`**: Additional plot for lung analysis (conserved gene for STANCE).

#### Utilities (real)
- **`real/utils/real-bench.R`**: Real-data benchmark helpers.
- **`real/utils/real_expr_bench.R`**: Expression / p-value correlation benchmark helpers.
- **`real/utils/rotate_bench.R`**: Rotation benchmark helpers.
- **`real/utils/rotate_subset.R`**:Subsetting helper functions for rotation.

---

### Simulation
Scripts for generating and analyzing simulated datasets.

#### Running simulations
- **`sim/1-runsim_sc.R`**: Run simulation pipeline for sc data and call ctSVG methods.
- **`sim/1-runsim_spot.R`**: Run simulation pipeline for spot data and call ctSVG methods.
- **`sim/1-runsim_noRCTD.R`**: Run simulation pipeline without RCTD.
- **`sim/fig3-permutation.R`**: Build permutation datasets and run ctSVG methods (for Fig 3).

#### Plotting (simulation & supplementary)
- **`sim/plot_fig3.R`**: Figure 3 — main simulation results.
- **`sim/plot_fig3f_permutation.R`**: Figure 3f — permutation results visualization.
- **`sim/plot_nodeconv.R`**: Plot deconvolution influence.
- **`sim/plot_Supplementary fig prop_fpr.R`**: Supplementary figure — plots the relationship between false positive rate (FPR) and cell-type proportion.
- **`sim/plot_Supplementary figs auc.R`**: Supplementary figures — AUC with different drop-outs.

#### Utilities (simulation)
- **`sim/utils/generate_sc.R`**: Generate simulated single-cell reference data.
- **`sim/utils/generate_st_P1.R`**: Generates simulated spatial transcriptomics data with 10% dropout.
- **`sim/utils/generate_st_P2.R`**: Generates simulated spatial transcriptomics data with 20% dropout.
- **`sim/utils/generate_st_P3.R`**: Generates simulated spatial transcriptomics data with 30% dropout.
- **`sim/utils/run_analysis_for_pattern.R`**: Runs ctSVG analysis for simulated data.
- **`sim/utils/run_analysis_for_pattern_sc.R`**: Runs ctSVG analysis for simulated single cell level spatial data.
- **`sim/utils/run_analysis_for_pattern_sp_noRCTD.R`**: Runs ctSVG analysis for simulated data (without deconvolution).
- **`sim/utils/calc_false_positive_rate.R`**: Compute false positive rate metrics.
- **`sim/utils/sim-bench.R`**: Simulation benchmark helpers.
- **`sim/utils/sim-bench-sc.R`**: Simulation benchmark helpers (sc spatial data).
- **`sim/utils/sim-bench-nodeconv.R`**: Simulation benchmark helpers (without deconvolution).

---

### Theme and Styles
- **`my_theme.R`**: Custom ggplot theme and style settings used across most figures.
