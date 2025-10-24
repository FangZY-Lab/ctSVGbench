# ctSVGbench: A systematic benchmarking of existing computational approaches for detecting cell-type-specific spatially variable genes in spatial transcriptomics data

## Project Structure

### Root Files

my_theme.R - Custom ggplot theme for consistent plotting


### Real Data Analysis (real/)

#### Main Analysis Scripts:

real datasets1.R, real datasets2.R - Main real data analysis pipelines

plot_fig2_realdata.R - Figure 2: Real data SVG detection results

plot_fig3f_null.R - Figure 3: Null model comparison

plot_fig4_rotate.R - Figure 4: Rotation robustness test

plot_fig5_time_mem.R - Figure 5: Computational performance (time/memory)

plot_fig6_lusc.R - Figure 6: LUSC cancer analysis

plot_fig7_mbm.R - Figure 7: Melanoma brain metastasis analysis

plot_fig8_summary.R - Figure 8: Summary figure

plot_figs1_dataset.R - Supplementary Figure 1: Dataset overview

#### Utility Functions (real/utils/):

real-bench.R - Main benchmarking function for real data

real_expr_bench.R - Benchmarking p-value and expression relation

rotate_bench.R - Rotation-based robustness testing

get_boundry_spVC.R - Boundary prepare for spVC method

get_pos.R - Save spatial positions 

### Simulation Analysis (sim/)

#### Main Scripts:

1-runsim.R - Main simulation pipeline

plot_sim.R - Main simulation results plotting

plot_sim_den_Sfig.R - Supplementary density figures for simulations

#### Utility Functions (sim/utils/):

sim-bench.R - Compare simulation results

generate_sc.R - Single-cell simulation data generation

generate_st_P1.R, generate_st_P2.R, generate_st_P3.R - Spatial pattern data generation (3 drop out paramter sets)

run_analysis_for_pattern.R - Run per simulation data

subset.R - Spots subsetting utilities

