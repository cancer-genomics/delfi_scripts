#!/bin/bash

#$ -cwd
#$ -j y
#$ -l mem_free=100G
#$ -l h_vmem=100G
#$ -l cancergen

#--------
# Input
#--------------------------------------------
rlib=/users/dbruhm/Library/R/3.4-bioc-devel
#--------------------------------------------

R_LIBS_USER=$rlib Rscript ./getScores.R
