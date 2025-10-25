# ctSVGbench: A systematic benchmarking of existing computational approaches for detecting cell-type-specific spatially variable genes in spatial transcriptomics data

## Project Structure

### Root Files

my_theme.R - Custom ggplot theme for consistent plotting


### **Real Data Analysis**  
Scripts for analyzing and visualizing real-world data.

- **`real/fig4_rotate.R`**: Rotates spots and then calls ctSVG.
- **`real/plot_fig2_realdata.R`**: Plots figure 2 for consistency analysis.
- **`real/plot_fig4_rotate.R`**: Plots figure 4 related to data rotation.
- **`real/plot_fig5_time_mem.R`**: Plots time and memory performance for figure 5.
- **`real/plot_fig6_lung_cancer.R`**: Plots figure 6 for lung cancer data.
- **`real/plot_fig7_mbm.R`**: Plots figure 7 for MBM data.
- **`real/plot_fig8_summary.R`**: Summary plot for method's performance.
- **`real/plot_figs1_dataset.R`**: Summarizes the dataset used in figure 2 and figure 5.
- **`real/preprocess_3_subset_spots.R`**: Preprocesses data for subset spots.
- **`real/preprocess_2_bounary.R`**: Performs boundary preprocessing for real data. Required by spVC.
- **`real/preprocess_1.3_lung_cancer.R`**: Lung cancer data preprocessing and subset creation.
- **`real/preprocess_1.2_public_data.R`**: Preprocesses public datasets.
- **`real/real_datasets1.R`**: Loads and analyzes real datasets without subsetting spots.
- **`real/real_datasets2.R`**: Loads and analyzes real datasets with spot subsetting.
- **`real/utils/real-bench.R`**: Benchmarking utilities for real data.
- **`real/utils/real_expr_bench.R`**: Expression and p-value correlation benchmarking utilities.
- **`real/utils/rotate_bench.R`**: Rotation benchmarking utilities.

### **Simulation**  
Scripts for generating and analyzing simulated data.

- **`sim/1-runsim.R`**: Runs simulation and calls ctSVG.
- **`sim/fig3-permutation.R`**: Generates permutation-based dataset and calls ctSVG.
- **`sim/plot_fig3f_permutation.R`**: Plots permutation results for figure 3f.
- **`sim/plot_sim.R`**: General plotting of simulation results.
- **`sim/plot_sim_den_Sfig.R`**: Density plot for supplementary figure about specificity in minority and majority cell types.
- **`sim/utils/generate_sc.R`**: Generates simulated single-cell data.
- **`sim/utils/generate_st_P1.R`**: Generates simulated spatial transcriptomics data with 10% dropout.
- **`sim/utils/generate_st_P2.R`**: Generates simulated spatial transcriptomics data with 20% dropout.
- **`sim/utils/generate_st_P3.R`**: Generates simulated spatial transcriptomics data with 30% dropout.
- **`sim/utils/run_analysis_for_pattern.R`**: Runs analysis for simulated data patterns.
- **`sim/utils/sim-bench.R`**: Benchmarking utilities for simulation.

### **Theme and Styles**
- **`my_theme.R`**: Custom theme settings for consistent plotting styles across all figures.
