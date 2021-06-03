#!/bin/bash

#$ -cwd
#$ -j y
#$ -l mem_free=100G
#$ -l h_vmem=100G

#--------
# Input
#--------------------------------------------
rlib=<PATH TO R LIBRARY>
#--------------------------------------------

R_LIBS_USER=$rlib Rscript ./getScores.R
