#### dmrscaler_simulation_job_submitter.sh START ####
#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o joblog.$JOB_ID.$TASK_ID
#$ -j y
## Edit the line below as needed:
#$ -l highp,h_rt=2:00:00,h_data=4G
## Modify the parallel environment
## and the number of cores as needed:
#$ -pe shared 1
# Email address to notify
#$ -M $USER@mail
# Notify when
#$ -m bea



# echo job info on joblog:
echo "Job $JOB_ID started on:   " `hostname -s`
echo "Job $JOB_ID started on:   " `date `
echo " "

# load the job environment:
. /u/local/Modules/default/init/modules.sh
## Edit the line below as needed:
module load R/4.1.0

## run setup to generate 1. Rdata that feeds into everything 2. simulation table 3. method table
#Rscript simul_setup.R

NUM_SIMUL_SETS=`wc -l < simul_table.csv`
NUM_METHOD_SETS=`wc -l < method_table.csv`


for SIMUL_SET_ID in {1..$NUM_SIMUL_SETS}
  do
  for METHOD_SET_ID in {1..$NUM_METHOD_SETS}
  do
    echo "SIMUL SET: $SIMUL_SET_ID, METHOD_SET: $METHOD_SET_ID"
    qsub dmrscaler_simulation_job_individual.sh $SIMUL_SET_ID $METHOD_SET_ID
  done
done




# echo job info on joblog:
echo "Job $JOB_ID ended on:   " `hostname -s`
echo "Job $JOB_ID ended on:   " `date `
echo " "
#### submit_job.sh STOP ####
