# Bash script to merge callsets given in an input panel file (first input argument)
# and path to SURVIVOR executable (second input argument)

PANEL=$1
SURVIVOR_PATH=$2
# path on cornell cluster: /athena/ihlab/scratch/vpopic/software/SURVIVOR/Debug/SURVIVOR
DIST=5000
MIN_MERGE_SIZE=1000
#MIN_FILT_SIZE=5000
DIR="$(dirname ${PANEL})"

VCFS=${DIR}/vcf_panel.txt
SAMPLES=${DIR}/samples_panel.txt

cut -f 1 ${PANEL} > ${SAMPLES}
cut -f 2 ${PANEL} > ${VCFS}

${SURVIVOR_PATH} merge ${VCFS} ${DIST} 1 1 -1 -1 ${MIN_MERGE_SIZE} ${DIR}/merged.vcf
#python /athena/ihlab/scratch/vpopic/SVNet/scripts/postproc_vcf.py --in_vcf ${DIR}/merged.vcf --out_vcf ${DIR}/merged_filtered.vcf  --sv_size ${MIN_FILT_SIZE} --pass_only --sv_types DEL,INV,DUP --as-is