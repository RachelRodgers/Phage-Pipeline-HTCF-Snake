#!/bin/bash

# submit_phage_pipeline_wrapper.sbatch

mkdir -p ./slurm_out

sbatch submit_phage_pipeline_snake.sbatch
