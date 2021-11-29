#!/bin/bash

# simulate the data sets
cd semi_parametric_sim
Rscript -e "rmarkdown::render('datasets_simulation.Rmd')"
# assess the quality of the simulated data sets
Rscript -e "rmarkdown::render('assess_simdata_qual.Rmd')"
cd ..

## preprocess the simulated data sets
cd summarized_benchmark
Rscript -e "rmarkdown::render('prepare_datasets.Rmd')"
## compute the benchmark performance
Rscript -e "rmarkdown::render('benchmark.Rmd')"
## collect and plot the performance evaluations
Rscript -e "rmarkdown::render('evaluations.Rmd')"
## compute quality statistics on the simulated data sets after the pre-processing
Rscript -e "rmarkdown::render('preprocessed_dataset_sim_qual.Rmd')"

################ MS simulations ########
# simulate the data sets
cd semi_parametric_sim
Rscript -e "rmarkdown::render('datasets_simulation_MS.Rmd')"
# assess the quality of the simulated data sets
# Rscript -e "rmarkdown::render('assess_simdata_qual.Rmd')"
cd ..

## preprocess the simulated data sets
cd summarized_benchmark
Rscript -e "rmarkdown::render('prepare_datasets_tmpl.Rmd', params = list(outdir = \"MS_prep_ds\", simulated_datasets_qs = \"../semi_parametric_sim/MS_sims/trimmed_simulated_datasets.qs\"), output_file = \"prepare_datasets_MS.html\")"
## compute the benchmark performance
# Rscript -e "rmarkdown::render('benchmark.Rmd')"
## collect and plot the performance evaluations
# Rscript -e "rmarkdown::render('evaluations.Rmd')"
## compute quality statistics on the simulated data sets after the pre-processing
# Rscript -e "rmarkdown::render('preprocessed_dataset_sim_qual.Rmd')"

################ DM1 simulations ########
# simulate the data sets
cd semi_parametric_sim
Rscript -e "rmarkdown::render('datasets_simulation_DM1.Rmd')"
# assess the quality of the simulated data sets
# Rscript -e "rmarkdown::render('assess_simdata_qual.Rmd')"
cd ..

## preprocess the simulated data sets
cd summarized_benchmark
Rscript -e "rmarkdown::render('prepare_datasets_tmpl.Rmd', params = list(outdir = \"DM1_prep_ds\", simulated_datasets_qs = \"../semi_parametric_sim/DM1_sims/trimmed_simulated_datasets.qs\"), output_file = \"prepare_datasets_DM1.html\")"
## compute the benchmark performance
# Rscript -e "rmarkdown::render('benchmark.Rmd')"
## collect and plot the performance evaluations
# Rscript -e "rmarkdown::render('evaluations.Rmd')"
## compute quality statistics on the simulated data sets after the pre-processing
# Rscript -e "rmarkdown::render('preprocessed_dataset_sim_qual.Rmd')"

################ IDC simulations ########
# simulate the data sets
cd semi_parametric_sim
Rscript -e "rmarkdown::render('datasets_simulation_IDC.Rmd')"
# assess the quality of the simulated data sets
# Rscript -e "rmarkdown::render('assess_simdata_qual.Rmd')"
cd ..

## preprocess the simulated data sets
cd summarized_benchmark
Rscript -e "rmarkdown::render('prepare_datasets_tmpl.Rmd', params = list(outdir = \"IDC_prep_ds\", simulated_datasets_qs = \"../semi_parametric_sim/IDC_sims/trimmed_simulated_datasets.qs\"), output_file = \"prepare_datasets_IDC.html\")"
## compute the benchmark performance
# Rscript -e "rmarkdown::render('benchmark.Rmd')"
## collect and plot the performance evaluations
# Rscript -e "rmarkdown::render('evaluations.Rmd')"
## compute quality statistics on the simulated data sets after the pre-processing
# Rscript -e "rmarkdown::render('preprocessed_dataset_sim_qual.Rmd')"

################ IPF simulations ########
# simulate the data sets
cd semi_parametric_sim
Rscript -e "rmarkdown::render('datasets_simulation_IPF.Rmd')"
# assess the quality of the simulated data sets
# Rscript -e "rmarkdown::render('assess_simdata_qual.Rmd')"
cd ..

## preprocess the simulated data sets
cd summarized_benchmark
Rscript -e "rmarkdown::render('prepare_datasets_tmpl.Rmd', params = list(outdir = \"IPF_prep_ds\", simulated_datasets_qs = \"../semi_parametric_sim/IPF_sims/trimmed_simulated_datasets.qs\"), output_file = \"prepare_datasets_IPF.html\")"
## compute the benchmark performance
# Rscript -e "rmarkdown::render('benchmark.Rmd')"
## collect and plot the performance evaluations
# Rscript -e "rmarkdown::render('evaluations.Rmd')"
## compute quality statistics on the simulated data sets after the pre-processing
# Rscript -e "rmarkdown::render('preprocessed_dataset_sim_qual.Rmd')"
