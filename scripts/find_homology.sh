#!/bin/bash
# find homology between how many references you may like

# requires:
# samtools faidx
# blastn

# if interested in human homology assume the user is using the STR panel
# from LineageTracker this is the only flag the other will be positional arguments

# length of flanking regions to detect homology might be a param
flank=200
max_motif=6
max_len=100

# usage print if -h
usage() {
  echo "Usage: $(basename $0) species_a_Y.fa [...] species_n_Y.fa" 2>&1
  echo 'Identify homology of Y-STRS from Y chromosome references.'
  echo "[ -f flanking ]         how many bp to add (default to 200)"
  echo "[ -e exclude ]          BED of regions to exclude (optional)"
  echo "[ -m max_motif ]        maximum motif length (default: 6)"
  echo "[ -M max_len ]          maximum STR length (default: 100)"
	echo "[ -p panel ]            also include human named panel from YLineageTracker"
	echo "[ -r reference ]        if -p is given the user must also give the human grch38 Y reference!"
	echo "[ References ]          as many fastas as you are interested in"
  exit 2
}

exit_abnormal() {
  usage
  exit 1
}

command_exists() {
    # check if command exists and fail otherwise
    command -v "$1" >/dev/null 2>&1
    if [[ $? -ne 0 ]]; then
        echo "$1 is required but it's not installed. Abort."
        exit 1
    fi
}

fastas=()

while [ $# -gt 0 ]
do
	unset OPTIND
	unset OPTARG
	while getopts :f:m:M:e:p:r: options
		do
			case $options in
				f) flank=${OPTARG}
				;;
				p) panel=${OPTARG}
				;;
				m) max_motif=${OPTARG}
				;;
				M) max_len=${OPTARG}
				;;
				e)
					exclude=${OPTARG}
				;;
				r)
					reference=${OPTARG}
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
		fastas+=("${!OPTIND}")
		shift $((OPTIND-1))
	shift
done

# check if panel file exists
if [ ! -f $panel ]; then echo "panel file $panel not found"; exit 1; fi

# check if blastn is installed
command_exists "blastn"

# if exclude is active check if bedtools is installed
if [ ! $exclude ]; then
  command_exists "bedtools"
fi

# exit if user does not give any fasta
if [ -z "$fastas" ]; then
    echo "No references given"
		exit 1
fi

# in this file we'll save various length, name and position
# of each STR
> regions_metadata.txt
# # and we will concatenate all the fasta in this file
> all_Y_strs.fa


get_fasta_str(){
	Y_id=$(head -1 $1 | cut -f1 -d' ' | sed 's/>//')
	# take name
	Y_name=$(basename $1)
	Y_name="${Y_name%.*}"
	num=1
	# for every region
	cat $2 | while read line; do
	start=$(echo $line | cut -f2 -d ' ');
	start_new="$(($start - $flank))";
	end=$(echo $line | cut -f3 -d' ');
	end_new="$(($end + $flank))";
	name=$(echo "${Y_name}_$num")
	length=$(expr $end - $start)
	echo -e "$Y_id\t$start\t$end\t$length\t$name" >> regions_metadata.txt
	samtools faidx $1 $Y_id:${start_new}-$end_new | sed "s/$Y_id:$start_new-$end_new/$name/g" >> all_Y_strs.fa
	let num++
	done
}

if [ ! -z $panel ]; then
	if [ -z $reference ]; then
		echo "panel option activated but no human reference given!" 1>&2
		exit 1
	fi
	echo "Including human named STRs" 1>&2
	# get human id
	Y_id=$(head -1 $reference | cut -f1 -d' ' | sed 's/>//')
	Y_name=$(basename $reference)
	Y_name="${Y_name%.*}"

	# the human panel is peculiar. It alread has names and the script is a little bit different
	awk 'BEGIN{FS=","} { if ($2!="\.") { print } }' $panel | cut -f2,5,6 -d',' | grep -v Start | tr ',' '\t'  | while read line; do
		start=$(echo $line | cut -f2 -d' ');
		start_new="$(($start - $flank))";
		end=$(echo $line | cut -f3 -d' ');
		end_new="$(($end + $flank))";
		name=$(echo $line | cut -f1 -d' ');
		length=$(expr $end - $start)
		echo -e "$Y_id\t$start\t$end\t$length\t$name" >> regions_metadata.txt
		samtools faidx $reference $Y_id:${start_new}-$end_new | sed "s/$Y_id:$start_new-$end_new/$name/g" >> all_Y_strs.fa;
		done

fi



for Y_ref in "${fastas[@]}"; do
	echo $Y_ref
	# take fasta id to produce regions file
	Y_id=$(head -1 $Y_ref | cut -f1 -d' ' | sed 's/>//')
	Y_name=$(basename $Y_ref)
	Y_name="${Y_name%.*}"

	echo "Looking for STR across $Y_name with trf" 1>&2
	echo "filtering for motifs shorter than $max_motif and regions shorter then $max_len" 1>&2

	# run trf with reccomended parameters and max motif
	trf $Y_ref 2 5 7 80 10 80 $max_motif -h -d
	# name of regions file
	filtered_bed=temp_regions.bed

	# Filter regions, if exclude is given remove STR overlapping excluded regions
	if [ ! $exclude ]; then
		sed '1,15d' $(basename ${Y_ref}).2.5.7.80.10.80.${max_motif}.dat | \
		awk -v Y_id="$Y_id" -v max_motif="$max_motif" -v max_len="$max_len" \
		 'BEGIN{OFS="\t"} {if ($3 > 2 && $3 <= max_motif && $2-$1 <= max_len)  \
		 {print Y_id, $1, $2, $3, $4}}' | \
		awk '!seen[$3]++' > $filtered_bed
		rm $(basename ${Y_ref}).2.5.7.80.10.80.${max_motif}.dat
	else
		sed '1,15d' $(basename ${Y_ref}).2.5.7.80.10.80.${max_motif}.dat | \
		awk -v Y_id="$Y_id" -v max_motif="$max_motif" -v max_len="$max_len" \
		 'BEGIN{OFS="\t"} {if ($3 > 2 && $3 <= max_motif && $2-$1 <= max_len)  \
		 {print Y_id, $1, $2, $3, $4}}' | \
		awk '!seen[$3]++' |\
		bedtools intersect -v -a stdin -b $exclude > $filtered_bed
		rm $(basename ${Y_ref}).2.5.7.80.10.80.${max_motif}.dat
	fi
	get_fasta_str $Y_ref $filtered_bed
	rm $filtered_bed
done



blastn -query all_Y_strs.fa -subject all_Y_strs.fa -outfmt 6 > blast_all_v_all_Y_str.csv

#
