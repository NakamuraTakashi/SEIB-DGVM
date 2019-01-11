#!/bin/sh

############################################################################
#   Script for running a MPI job @ GCOE saver
#
# ___ Notice ___
#   (1) Execusion command -> sh go.sh &
#   (2) Attribute of this file must be '744'
#   (3) Once executed, environmental variables will be retained until logout
#   (4) Use Makefile for making a execution file for MPI jobs
############################################################################

#Submit MPI program, and wait for complete
   qsub go.bat
   echo "Waiting for MPI comuptation"
   date
   while [[ `qstat | grep hsato | wc -l` -ge 1  ]]; do
     sleep 5
   done
   date
   
#Visualizing result files with R script
   echo "Visualizing result files"
   R --vanilla --slave --quiet < code_visualize/view1d.R
   
#message
   echo "Normal End"
