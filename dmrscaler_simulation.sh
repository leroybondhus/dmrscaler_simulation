#### submit_job.sh START ####
#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o joblog.$JOB_ID
#$ -j y
## Edit the line below as needed:
#$ -l highp,h_rt=1:00:00,h_data=2G
## Modify the parallel environment
## and the number of cores as needed:
#$ -pe shared 100
# Email address to notify
#$ -M bondh185@ucla.edu
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

## substitute the command to run your code
## in the two lines below:

## run setup to generate 1. Rdata that feeds into everything 2. simulation table 3. method table



Rscript simul_setup.R

NUM_SIMUL_SETS=`wc -l < simul_table.csv`
NUM_METHOD_SETS=`wc -l < method_table.csv`

for i in {1..$NUM_SIMUL_SETS}
do
		for j in {1..$NUM_METHOD_SETS}
		do
			Rscript simul_individual_run.R $NUM_SIMUL_SETS $NUM_METHOD_SETS  & ## include ampersand to run in parallel	
		done
done

Rscript 


Rscript test_dopar.R

# echo job info on joblog:
echo "Job $JOB_ID ended on:   " `hostname -s`
echo "Job $JOB_ID ended on:   " `date `
echo " "
#### submit_job.sh STOP ####
