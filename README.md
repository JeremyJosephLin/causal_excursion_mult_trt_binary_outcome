# paper_Causal_excursion_mult_trt_binary_outcome

Reproducible code for the paper "Causal excursion effect for continuous outcome with multilevel treatment with binary Outcome" by Jeremy Lin and Tianchen Qian. For an example of using the wcls estimator for continuous treatment, run the R codes in the "misc/toy-example" folder.

## File structure

- application code: code to reproduce results in Section 5 "Data Example".
- functions: code containing relevant user-defined R functions that includes sample size calculator and CEE estimator
- simulation code: code to reproduce simulation results in Section 4 "Sample Size Formula with Categorical Treatment" (Sec. 4.3) and Appendix C and D.

## How to reproduce the results

To reproduce the results in Section 5, "Data Example," run all R scripts in the folder "application code." The R scripts do not depend on one another, so there is no particular order for running.
- To reproduce simulation results in Section  "Detailed Simulation Results" (Appendix D), do the following for each subfolder inside "simulation code":
    - First, run the R script(s) inside all subfolder(s) "simulationX(.X)". This conducts Monte Carlo simulations and saves result file. **Caution: each R script may take a long time (days) to finish. Also, X.X in simulationX.X does not correspond to Sections in the paper; these indices were created during the development of the paper.**
    - Second, run the Rmd file(s) named "sim_X_X_report.RMD". This makes plots and results using the simulation result files. **Caution: X in the "make figure X.R" does not correspond to figure index in the paper. See table below for the figure index correspondence.**
    - For example, to reproduce everything in Section D.4 in the appendix, go inside the folder "2. WA-a violated", then run the R scripts in subfolders "sim 2.1a", "sim 2.2a", "sim 2.3a", "sim 2.4a", "sim 2.5a", "sim 2.6a", and "sim 3.1a". Then, run the Rmd scripts "sim_X_X_report.RMD".


| Figure in paper | R script to make the figure                                                            |
|-----------------|----------------------------------------------------------------------------------------|
| 1.a             | application code/graph_constant.Rmd                                                    |
| 1.b             | application code/theta_sensitivity_graph.Rmd                                           |
| 1.c,d           | application code/AA_sensitivity_graph.Rmd                                              |
| S1, S2          | simulation/power simulation/1. all working assumptions hold/sim_1a_report.Rmd          |
| S3              | simulation/power simulation/2.WA-a violated/sim 3.1a/sim_3_1a_report.Rmd               |
| S4a             | simulation/power simulation/2.WA-a violated/sim 2.1a/sim_2_1_report.Rmd                |
| S4b             | simulation/power simulation/2.WA-a violated/sim 2.2a/sim_2_2_report.Rmd                |
| S4c             | simulation/power simulation/2.WA-a violated/sim 2.5a/sim_2_5_report.Rmd                |
| S4d             | simulation/power simulation/2.WA-a violated/sim 2.3a/sim_2_3_report.Rmd                |
| S4e             | simulation/power simulation/2.WA-a violated/sim 2.4a/sim_2_4_report.Rmd                |
| S4f             | simulation/power simulation/2.WA-a violated/sim 2.6a/sim_2_6_report.Rmd                |
| S5a             | simulation/power simulation/3.WA-b violated/result/sim_4_1_report.Rmd                  |
| S5b             | simulation/power simulation/3.WA-b violated/sim_const_true_eo.Rmd                      |
| S5c             | simulation/power simulation/3.WA-b violated/sim_const_working_eo.Rmd                   |
| S6a, b          | simulation/power simulation/3.WA-b violated/sim_4_7a_report.Rmd                        |
| S6c, d          | simulation/power simulation/3.WA-b violated/sim_4_7b_report.Rmd                        |
| S7              | simulation/power simulation/4. WA-c violated/sim 5.3/sim_5_3a_report.Rmd               |
| S8a             | simulation/power simulation/4. WA-c violated/sim 5.3/sim_5_2a_report.Rmd               |
| S8b             | simulation/power simulation/4. WA-c violated/sim 5.3/sim_5_1a_report.Rmd               |
| S9              | simulation/power simulation/5. WA-d violated/sim_wa_D_report.Rmd                       |  
| S10a, b, c      | simulation/power simulation/5. WA-e violated/sim_wa_E_report.Rmd                       |

