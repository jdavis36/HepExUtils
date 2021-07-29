#!/bin/bash
sbatch -o slurm/slurm-%A_%a.out --array=1-164 decay_same_sign_run.sh
