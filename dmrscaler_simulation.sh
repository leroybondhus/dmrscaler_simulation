#### dmrscaler_simulation.sh START ####
#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o joblog.$JOB_ID.$TASK_ID
#$ -j y
## Edit the line below as needed:
#$ -l highp,h_rt=2:00:00,h_data=4G
## Modify the parallel environment
## and the number of cores as needed:
#$ -pe shared 100
# Email address to notify
#$ -M $USER@mail
# Notify when
#$ -m bea
# Job array indexes

NUM_SIMUL_SETS=`wc -l < simul_table.csv`
NUM_METHOD_SETS=`wc -l < method_table.csv`
UPPER_LIM=$(expr $NUM_SIMUL_SETS \* $NUM_METHOD_SETS )
#$ -t 1-${UPPER_LIM}:1

# echo job info on joblog:
echo "Job $JOB_ID started on:   " `hostname -s`
echo "Job $JOB_ID started on:   " `date `
echo " "

# load the job environment:
. /u/local/Modules/default/init/modules.sh
## Edit the line below as needed:
module load R/4.1.0

## substitute the command to run your code
## in the two lines below:

## run setup to generate 1. Rdata that feeds into everything 2. simulation table 3. method table



Rscript simul_setup.R

SIMUL_SET_ID=$(expr $SGE_TASK_ID % $NUM_SIMUL_SETS )
METHOD_SET_ID=$(expr 1 + $(expr $SGE_TASK_ID / $NUM_SIMUL_SETS ))

Rscript simul_individual_run.R $SIMUL_SET_ID $METHOD_SET_ID

Rscript


Rscript test_dopar.R

# echo job info on joblog:
echo "Job $JOB_ID ended on:   " `hostname -s`
echo "Job $JOB_ID ended on:   " `date `
echo " "
#### submit_job.sh STOP ####
