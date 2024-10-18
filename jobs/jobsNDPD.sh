#!/bin/bash
#! This line is a comment
#! Make sure you only have comments and #SBATCH directives between here and the end of the #SBATCH directives, or things will break
#! Name of the job:
#SBATCH -J test_Porter5
#! Account name for group, use SL2 for paying queue:
#SBATCH -A AHNERT-SL3-CPU
#! Output filename:
#! %A means slurm job ID and %a means array index
#SBATCH --output=test_dmtcp_%A_%a.out
#! Errors filename:
#SBATCH --error=test_dmtcp_%A_%a.err
#! memory
#SBATCH --mem=120G
#! Number of nodes to be allocated for the job (for single core jobs always leave this at 1)
#SBATCH --nodes=1
#! Number of tasks. By default SLURM assumes 1 task per node and 1 CPU per task. (for single core jobs always leave this at 1)
#SBATCH --ntasks=1
#! How many many cores will be allocated per task? (for single core jobs always leave this at 1)
#SBATCH --cpus-per-task=4
#! Estimated runtime: hh:mm:ss (job is force-stopped after if exceeded):
#SBATCH --time=10:00:00
#! Estimated maximum memory needed (job is force-stopped if exceeded):
#! RAM is allocated in ~5980mb blocks, you are charged per block used,
#! and unused fractions of blocks will not be usable by others.
#! Submit a job array with index values between 0 and 31
#! NOTE: This must be a range, not a single number (i.e. specifying '32' here would only run one job, with index 32)
#!SBATCH --array=1500-5000
#!SBATCH --array=48-270
#SBATCH --array=0-166

#! This is the partition name.
#SBATCH -p icelake-himem

#! mail alert at start, end and abortion of execution
#! emails will default to going to your email address
#! you can specify a different email address manually if needed.
##SBATCH --mail-type=ALL

#! Don't put any #SBATCH directives below this line

#! Modify the environment seen by the application. For this example we need the default modules.
. /etc/profile.d/modules.sh                # This line enables the module command
module purge                               # Removes all modules still loaded
module load rhel7/default-peta4            # REQUIRED - loads the basic environment

source .bashrc
conda init bash
conda activate biopython

#! The variable $SLURM_ARRAY_TASK_ID contains the array index for each job.
#! In this example, each job will be passed its index, so each output file will contain a different value
echo "This is job" $SLURM_ARRAY_TASK_ID

#! Command line that we want to run:
jobDir=Job_$SLURM_ARRAY_TASK_ID
echo $SLURM_ARRAY_TASK_ID
python3 quantitiesND4PD.py $SLURM_ARRAY_TASK_ID
