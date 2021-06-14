#$ -l h_vmem=5G
module load nextflow

SCRIPT=scripts/vHRR_single.nf
CONFIG=nextflow_config/vHRR_single.config

nextflow run $SCRIPT -c $CONFIG -resume
