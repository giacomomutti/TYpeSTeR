# TYpeSTeR

### 1. Map samples (optional)

If you don't already have the BAMs you can reproduce our analysis by running the snakemake pipeline in the pipeline folder.

You'll need the paired fastq files, a config file where to specify different directories and the reference genome and a json file where each sample is associated to a list of runs.

Check the example files in order to have an idea of what you may need.

The dependencies are:
* snakemake
* picard
* samtools
* bwa-mem2

The reference fasta must be indexed with `bwa-mem2 index $reference` and `samtools faidx $reference`.

### 2. Usage

TYpeSTeR can be used to find, filter and genotype Y chromosome Short Tandem Repeats on a set of samples.

explain dependencies
In order to get the putative Y STR region you'll need the y chromosome fasta and TRF software (available through conda, [Benson, 1999](https://github.com/Benson-Genomics-Lab/TRF)).

For this step you have to install [HipSTR](https://hipstr-tool.github.io/HipSTR/) ([Willems et al, 2017](https://www.nature.com/articles/nmeth.4267)).


In order to run TYpeSTeR you need:
1. the Y chromosome reference in fasta format (`-r`). The reference must have been indexed with `samtools faidx $reference` for HipSTR to work.
2. a file with the paths to your samples bam files, one per line (`-s`). The script will work as long as the Y ref id can be found in the BAM headers. You'll also need to add all the paths to your bam files in a single file for example like this: `find *Y.bam > Y_bams.txt`
3. the path to HipSTR executable file (`-t`). The programme assumes it's `HipSTR` but you may specify the full path if it's not in $PATH.
4. (**optional**) a bed file of regions to exclude (`-e`), such as Pseudo Autosomal Regions (this step requires bedtools)
5. (**optional**) a name for the output files of your analysis, the default is `Y_STRs`
6. Maximum motif length (`-m`), default to 6
7. Maximum STR length (`-M`), default to 100

Then you can run the script with:

`./typester.sh -r $ref_Y.fa -s $samples.txt -t /usr/bin/HipSTR -e $exclude.bed`

The script will first use trf with these parameters:
Match=2, Mismatch=5, Delta=7, PM=80, PI=10, MinScore=80 (higher than reccomended 50), MaxPeriod=6 (can be changed with `-m`).

Then it will filter these regions removing dinucleotides and repeats with motif larger than ` -m` (default: 6). Further, as HipSTR by default does not genotype motif longer than 100bp, repeats longer then that will be filtered. You can modify this parameter by setting the maximum length flag `-M`.

The HipSTR maximum flanking indel parameter (`--max-flank-indel`) is set to 1 as it can be eventually filtered downstream and you might be genotyping not so closely related species. Also the `--min-reads` parameter is the same as your number of samples. This is a very relaxed filtering but again it's better to eventually filter results downstream.

Run `./typester.sh -h` to get additional details on usage.

*Depending on your project you may decide which kind of stutter model to use, a nice guide can be found [here](https://hipstr-tool.github.io/HipSTR/#in-depth-usage). In this case we will use the default stutter model (`--def-stutter-model`) as our samples are from different sequencing projects and have different coverages.*

### 3. Analyze results

explain Rscripts

### 4. Detect homology between STRs (optional)

In case you have multiple y chromosomes which are more or less closely related you can use the script in `scripts/find_homology.sh` to look for putative homology between them.

This script is based on the same `TRF` command as before and you can use the same flags `-m -M` as in `typester.sh` to filter results. You can use this command with how many references you like as in:

`find_homology.sh species_a_Y.fa [...] species_n_Y.fa`

The script will produce three files:
1. `all_Y_strs.fa`: the fasta files of the filtered STRs found in your references. The sequences will be named as filename+ an increasing integer based on their postion (i.e. given an input called Hsap_Y.fa the regions will ba named Hsap_Y_1,2...n). Further, each sequence will have 200 bp additional base paires in the flanking regions (you can change this parameter with `-f`).
2. `regions_metadata.txt`: file with 5 columns: reference_id, start, end, length, name
3. `blast_all_v_all_Y_str.csv`: blast results of the all vs all search (outfmt 6).

You can check additional parameters with `./scripts/find_homology.sh -h`

If you are interested in homology with Human Y-STRs, you might use the Y_STR_panel.csv from [YLineageTracker](https://github.com/Shuhua-Group/Y-LineageTracker) with the `-p` option. This is good as human STRs are well studied and have official names. If you use this you have to specify the fasta Human Y chromosome GRCh38 reference (CM000686.2 or NC_000024.10) with `-r`!

This script also allows to search for multicopies STR within a single reference!

Finally, you might visualize the blast results with:

`Rscript scripts/comparative_STR.R`

This will produce a pdf named `blast_Y_STRs.pdf` with a pairwise comparison of blast results of all the references analyzed.

# TODOs

* fix README
* add scripts to viz vcf and regions
* nicer automatic color palette
* smarter height width plots
