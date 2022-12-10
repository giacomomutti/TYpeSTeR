#!/bin/bash

# input required by user are:
# Y fasta reference
Y_ref=""
# list of bam files
samples=""
HipSTR="HipSTR"
max_motif=6
max_len=100
skip=false

# usage print if -h
usage() {
  echo "Usage: $(basename $0) -r ref_y -s samples.txt" 2>&1
  echo 'Identify and genotype Y-STRS from WGS samples.'
  echo "[ -r Ref_Y ]            Y chromosome reference fasta"
  echo "[ -s samples ]          file with bam file paths"
  echo "[ -t HipSTR ]           path to compiled HipSTR (optional if in \$PATH)"
  echo "[ -e exclude ]          BED of regions to exclude (optional)"
  echo "[ -o outprefix ]        prefix of output files (optional, default: Y_STRs)"
  echo "[ -m max_motif ]        maximum motif length (default: 6)"
  echo "[ -M max_len ]          maximum STR length (default: 100)"
  echo "[ -n no_geno ]          skip genotyping with HipSTR"
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


while getopts ":r:s:t:e:o:m:M:n" options; do
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
      outprefix=${OPTARG}
      ;;
    m)
      max_motif=${OPTARG}
      ;;
    M)
      max_len=${OPTARG}
      ;;
    n)
      skip=true
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

# check if trf installed
for COMMAND in "trf"; do
    command_exists "${COMMAND}"
done

# if exclude is active check if bedtools is installed
if [ ! $exclude ]; then
  command_exists "bedtools"
fi

# fasta is required and must exists
if [ -z $Y_ref ]; then echo "Y reference not given"; usage; exit 1; fi
if [ ! -f $Y_ref ]; then echo "Y reference file $Y_ref not found"; exit 1; fi

# can hipstr be found?
if [[ "$skip" == false ]]; then
  # samples is required and must exists
  if [ -z $samples ]; then echo "samples file not given"; usage; exit 1; fi
  if [ ! -f $samples ]; then echo "samples file $samples not found"; exit 1; fi
  command_exists $HipSTR
fi

# if prexif not given default to Y_STRss
if [ ! $outprefix ]; then
  outprefix="Y_STRs"
fi

# exit if max motif or max len are not integers
if ! [[ "$max_motif" =~ ^[0-9]+$ ]];
   then echo "error: -m $max_motif is not an integer" >&2; exit 1
fi

if ! [[ "$max_len" =~ ^[0-9]+$ ]];
   then echo "error: -M $max_len is not an integer" >&2; exit 1
fi

echo "Looking for STR across $Y_ref on $samples with trf" 1>&2
echo "filtering for motifs shorter than $max_motif and regions shorter then $max_len" 1>&2

# take fasta id to produce regions file
Y_id=$(head -1 $Y_ref | cut -f1 -d' ' | sed 's/>//')
Y_name=$(basename $Y_ref)
Y_name="${Y_name%.*}"

# run trf with reccomended parameters and max motif
trf $Y_ref 2 5 7 80 10 80 $max_motif -h -d
# name of regions file
filtered_bed=${Y_name}_${max_motif}_${max_len}_filtered_STR.bed

# Filter regions, if exclude is given remove STR overlapping excluded regions
if [ ! $exclude ]; then
  sed '1,15d' $(basename ${Y_ref}).2.5.7.80.10.80.${max_motif}.dat | \
  awk -v Y_id="$Y_id" -v max_motif="$max_motif" -v max_len="$max_len" \
   'BEGIN{OFS="\t"} {if ($3 > 2 && $3 <= max_motif && $2-$1 <= max_len)  \
   {print Y_id, $1, $2, $3, $4}}' | \
  awk '!seen[$3]++' > $filtered_bed
else
  sed '1,15d' $(basename ${Y_ref}).2.5.7.80.10.80.${max_motif}.dat | \
  awk -v Y_id="$Y_id" -v max_motif="$max_motif" -v max_len="$max_len" \
   'BEGIN{OFS="\t"} {if ($3 > 2 && $3 <= max_motif && $2-$1 <= max_len)  \
   {print Y_id, $1, $2, $3, $4}}' | \
  awk '!seen[$3]++' |\
  bedtools intersect -v -a stdin -b $exclude > $filtered_bed
fi

# how many STRs?
STRs_number=$(wc -l < $filtered_bed)

echo -e "\n$STRs_number STRs were found after filtering" 1>&2

if [[ "$skip" == true ]]; then
    echo "Skipping genotyping"
    exit 0
fi

echo -e "\nRunning HipSTR to genotype the samples"

# relaxed filter of min reads per str
nsamp=$(wc -l < $samples)

# if the str is exactly as max_len hipstr will ignore it therefore add 1
max_len_h=$((max_len + 1))

# run hipstr
$HipSTR --bam-files $samples --fasta $Y_ref --regions $filtered_bed --haploid-chrs $Y_id \
--str-vcf ${outprefix}.vcf.gz --log ${outprefix}.log --lib-from-samp --output-filters \
--max-flank-indel 1 --def-stutter-model --min-reads $nsamp --max-str-len $max_len_h

# print some useful stats on hipstr log
no_reads=$(grep -c "Skipping locus with too few reads" ${outprefix}.log)
too_rep=$(grep -c "too repetitive"  ${outprefix}.log)

grep "succeeded" ${outprefix}.log

echo "$no_reads had too few reads (less than the number of samples: $nsamp), in $too_rep the sequence upstream of the repeat is too repetitive for accurate genotyping"

exit 0
