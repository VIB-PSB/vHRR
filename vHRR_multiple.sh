#$ -l h_vmem=5G
module load nextflow

SCRIPT=scripts/vHRR_multiple.nf
CONFIG=nextflow_config/vHRR_multiple.config

nextflow run $SCRIPT -c $CONFIG -resume
