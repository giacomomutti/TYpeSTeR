#!/bin/bash

# input required by user are:
# Y fasta reference
Y_ref=""
# list of bam files
samples=""
HipSTR="HipSTR"

usage() {
  echo "Usage: $(basename $0) [ -r Ref_Y ]
                  [ -s samples ]
                  [ -t HipSTR ]
                  [ -e exclude ]
                  [ -o output ]"
  exit 2
}

exit_abnormal() {
  usage
  exit 1
}

# check if dependencies exists
# https://raymii.org/s/snippets/Bash_Bits_Check_if_command_is_available.html
command_exists() {
    # check if command exists and fail otherwise
    command -v "$1" >/dev/null 2>&1
    if [[ $? -ne 0 ]]; then
        echo "$1 is required but it's not installed. Abort."
        exit 1
    fi
}


while getopts ":r:s:t:e:o:" options; do
  case "${options}" in
    r)
      Y_ref=${OPTARG}
      ;;
    s)
      samples=${OPTARG}
      ;;
    t)
      HipSTR=${OPTARG}
      ;;
    e)
      exclude=${OPTARG}
      ;;
    o)
      outname=${OPTARG}
      ;;
    :)
      echo "Error: -${OPTARG} requires an argument."
      exit_abnormal
      ;;
    *)
      exit_abnormal
      ;;
  esac
done


for COMMAND in "trf"; do
    command_exists "${COMMAND}"
done

if [ ! $exclude ]; then
  command_exists "bedtools"
fi

if [ -z $Y_ref ]; then echo "Y reference not given"; usage; exit 1; fi
if [ -z $samples ]; then echo "samples file not given"; usage; exit 1; fi
if [ -z $HipSTR ] && command_exists $HipSTR; then echo "HipSTR not found"; exit 1; fi

if [ ! -f $Y_ref ]; then echo "Y reference file $Y_ref not found"; exit 1; fi
if [ ! -f $samples ]; then echo "samples file $samples not found"; exit 1; fi

if [ ! $outname ]; then
  outname="Y_STRs"
fi


echo "Looking for STR across $Y_ref on $samples with trf" 1>&2

Y_id=$(head -1 $Y_ref | cut -f1 -d' ' | sed 's/>//')

# run trf, these parameters may be changed but then you'll need to change
# the next line
trf $Y_ref 2 5 7 80 10 80 6 -h -d

if [ ! $exclude ]; then
  sed '1,15d' ${Y_ref}*.dat | \
  awk -v Y_id="$Y_id" 'BEGIN{OFS="\t"} {if ($3 > 2 && $3 < 7 && $2-$1 < 101) { print Y_id, $1, $2, $3, $4}}' | \
  awk '!seen[$3]++' > ${Y_id}_filtered_STR.bed
else
  sed '1,15d' ${Y_ref}*.dat | \
  awk -v Y_id="$Y_id" 'BEGIN{OFS="\t"} {if ($3 > 2 && $3 < 7 && $2-$1 < 101) { print Y_id, $1, $2, $3, $4}}' | \
  awk '!seen[$3]++' |\
  bedtools intersect -v -a stdin -b $exclude > ${Y_id}_filtered_STR.bed
fi
# here you would remove regions if user gives bed with regions to exclude
STRs_number=$(wc -l < ${Y_id}_filtered_STR.bed)

echo -e "\n$STRs_number STRs were found after filtering" 1>&2

echo -e "\nRunning HipSTR to genotype the samples"

# eventually you could let the user decide if less or more min reads or prop of
# flanking indels, ESPECIALLY DEF STUTTER MODEL!
nsamp=$(wc -l < $samples)

$HipSTR --bam-files $samples --fasta $Y_ref --regions ${Y_id}_filtered_STR.bed --haploid-chrs $Y_id \
--str-vcf ${outname}.vcf.gz --log ${outname}.log --lib-from-samp --output-filters \
--max-flank-indel 1 --def-stutter-model --min-reads $nsamp

exit 0
