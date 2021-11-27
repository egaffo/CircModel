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
