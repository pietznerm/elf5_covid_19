#!/bin/bash

## script to obtain SNP summary statistics for SomaScan proteins and Covid-19

#!#############################################################
#!#### Modify the options in this section as appropriate ######
#!#############################################################

#! Name of the job:
#SBATCH -J snp_stats_SL_covid

#! Which project should be charged:
#SBATCH -A MRC-EPID-SL0-CPU

#! How many whole nodes should be allocated?
# In your case you should leave this at 1
#SBATCH --nodes=1

#SBATCH --ntasks=1

#SBATCH --cpus-per-task 4

##SBATCH --exclusive

##SBATCH --exclude=cpu-d-[1-5]

#SBATCH --array=1-1126%50

#! Specify required run time
#SBATCH --time=72:00:00

#! What types of email messages do you wish to receive?
#SBATCH --mail-type=FAIL

#SBATCH --output=slurm-%x-%j.out

#! Uncomment this to prevent the job from being requeued (e.g. if
#! interrupted by node failure or system downtime):
##SBATCH --no-requeue

#! Do not change:
#SBATCH -p epid

#! ############################################################
#! Modify the settings below to specify the application's environment, location
#! and launch method:

#! Optionally modify the environment seen by the application
#! (note that SLURM reproduces the environment at submission irrespective of ~/.bashrc):
. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module purge                               # Removes all modules still loaded
module load rhel7/default-peta4            # REQUIRED - loads the basic environment

#! Insert additional module load commands after this line if needed:

module load gcc/5
module load r-3.6.0-gcc-5.4.0-bzuuksv

## assign directories used for the analysis
cd ~/rfs/rfs-epid-rfs-mpB3sSsgAn4//Studies//People/Maik/COVID19_release6/02_naive_cis_coloc/

## get file name
export FL="input/${1}"

echo ${FL}

## run as array job
echo "Job ID: $SLURM_ARRAY_TASK_ID"
pheno="$(awk -v var="$SLURM_ARRAY_TASK_ID" -F '\t' 'NR == var {print $1}' ${FL})"
chr="$(awk -v var="$SLURM_ARRAY_TASK_ID" -F '\t' 'NR == var {print $2}' ${FL})"
LOWPOS="$(awk -v var="$SLURM_ARRAY_TASK_ID" -F '\t' 'NR == var {print $3}' ${FL})"
UPPOS="$(awk -v var="$SLURM_ARRAY_TASK_ID" -F '\t' 'NR == var {print $4}' ${FL})"
rsid="$(awk -v var="$SLURM_ARRAY_TASK_ID" -F '\t' 'NR == var {print $5}' ${FL})"


echo "Phenotype ${pheno} : Chromosome ${chr} : Locus start ${LOWPOS} : Locus end ${UPPOS} : rsID(s) ${rsid}"

#------------------------------#
## -->      run coloc     <-- ##
#------------------------------#

## run simple coloc with original association statistics
scripts/08_SNP_stats_SomaScan_COVID19.R ${pheno} ${chr} ${LOWPOS} ${UPPOS} "${rsid}"
