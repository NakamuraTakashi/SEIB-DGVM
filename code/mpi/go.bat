#!/bin/sh 
#$ -S /bin/sh
#$ -cwd            # Execute at current directory
#$ -V              # Inherit environmental variables for all machines in a cluster server
#$ -q all.q        # 
#$ -N seib_wide    # Define job name
#$ -j y            # Standard output and error message will be written in a same file
#$ -pe openmpi8 40 # Define processor number employing this job (per each node , total)

cd /home/hsato/work
mpirun -np $NSLOTS ./go.out  >& log.txt

