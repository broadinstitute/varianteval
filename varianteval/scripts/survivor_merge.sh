# Bash script to merge callsets given in an input panel file (first input argument)
# and path to SURVIVOR executable (second input argument)

PANEL=$1
SURVIVOR_PATH=$2
DIST=5000
MIN_MERGE_SIZE=1000
MIN_FILT_SIZE=5000
DIR="$(dirname ${PANEL})"

VCFS=${DIR}/vcf_panel.txt
SAMPLES=${DIR}/samples_panel.txt

cut -f 1 ${PANEL} > ${SAMPLES}
cut -f 2 ${PANEL} > ${VCFS}

SAMPLE_LIST=$(tr "\n" " " < ${SAMPLES})

${SURVIVOR_PATH} merge ${VCFS} ${DIST} 1 1 -1 -1 ${MIN_MERGE_SIZE} ${DIR}/merged.vcf

python upset_from_annotations.py --annot_vcf ${DIR}/merged.vcf --output_path ${DIR}/upset.png --tools ${SAMPLE_LIST} --plot_type upset