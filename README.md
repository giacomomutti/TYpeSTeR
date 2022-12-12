# TYpeSTeR

### Map samples (optional)

If you don't already have the mapped samples (in BAM format) you can reproduce our analysis by running the snakemake pipeline in the `pipeline` folder.

You'll need the paired fastq files, a config file where to specify different output directories, the reference genome and a json file where each sample is associated to a list of runs.

Check the example files in order to have an idea of what you may need.

The dependencies are:
* snakemake
* picard
* samtools
* bwa-mem2

The reference fasta must be indexed with `bwa-mem2 index $reference` and `samtools faidx $reference`.

### 1. Installation

First clone the github repository.

```
git clone https://github.com/giacomomutti/TYpeSTeR
cd TYpeSTeR

```

It depends on `TRF` in order to identify the putative Y STR regions in the Y chromosome reference ([Benson, 1999](https://github.com/Benson-Genomics-Lab/TRF)). `TRF` is available through conda:

```
conda install -c bioconda trf
```

[HipSTR](https://hipstr-tool.github.io/HipSTR/) is used to genotype the samples in the newly found regions ([Willems et al, 2017](https://www.nature.com/articles/nmeth.4267)). You can install this tool following the [installation instructions](https://github.com/HipSTR-Tool/HipSTR#installation).

### 2. Usage

TYpeSTeR can be used to find, filter and genotype Y chromosome Short Tandem Repeats (Y-STRs) on a set of WGS Illumina samples.

You can simply run the script with:

```
./typester.sh -r ref_Y.fa -s samples.txt -t /usr/bin/HipSTR -e regions_exclude.bed
```


In order to run TYpeSTeR you need:
1. the Y chromosome reference in fasta format (`-r`). The reference must have been indexed with `samtools faidx $reference` for HipSTR to work.
2. a file with the paths to your samples bam files, one per line (`-s`). The script will work as long as the Y reference sequence ID can be found in the BAM headers. You can add the paths in a single file for example like this: `find *Y.bam > Y_bams.txt`
3. the path to HipSTR executable file (`-t`). The program assumes it's `HipSTR` but you may specify the full path if it's not in `$PATH`.
4. (**optional**) a bed file of regions to exclude (`-e`), such as Pseudo Autosomal Regions (this step requires bedtools)
5. (**optional**) a name for the output files of your analysis, the default is `Y_STRs`
6. Maximum motif length (`-m`), default to 6
7. Maximum STR length (`-M`), default to 100


The script will first use `TRF` with these parameters:
`Match=2`, `Mismatch=5`, `Delta=7`, `PM=80`, `PI=10`, `MinScore=80` (higher than the reccomended 50), `MaxPeriod=6` (can be changed with `-m`).

Then it will filter these regions removing dinucleotides and repeats with motif larger than ` -m` (default: 6). Further, as HipSTR does not genotype motif longer than 100bp by default, repeats longer then that will be filtered. You can modify this parameter by setting the maximum length flag `-M`. If you are only interested in identifying the putative STRs you could use the flag `-n` to skip genotyping without providing any sample file `-s`.

The HipSTR maximum flanking indel parameter (`--max-flank-indel`) is set to 1 as it can be eventually filtered downstream and you might be genotyping not so closely related species. Also, the `--min-reads` parameter is the same as your number of samples. This is a very relaxed filtering but again, it's better to eventually filter results downstream.


`typester.sh` will output:

1. the raw `TRF` data file which will be named as `${Y_ref}.2.5.7.80.10.80.${max_motif}.dat`
2. the filtered STR regions named `${Y_name}_${max_motif}_${max_str}_filtered_STR.bed`
3. The HipSTR genotypes and log in `${outprefix}.vcf.gz` and `.log` . `${outprefix}` can be set with `-o` (Default: Y_STRs)

Run `./typester.sh -h` to get additional details on usage.

### 3. Analyze results

After you obtain the VCF and the TRF regions you can explore the results with the R script:

```
Rscript scripts/plot_STR_vcf.R -r $regions -t $trf -i $vcf -d $taxonomy
```

You will need to input the three output of TYpeSTeR: `-r` filtered regions, `-t` trf output and `-i` the VCF. Further, you have to input a file with taxonomic division. This file has to have two column: sample and taxon and the sample code must be the same as found in the VCF file. Check `Rscript scripts/plot_STR_vcf.R -h` for usage.

The script will output three pdf files in the `plots` directory: one plot will be the distribution of the motif size and the location across the reference, an exploratory plot of the VCF and the heatmap of each alleles length of the filtered STR regions.


### 4. Detect homology between STRs (optional)

In case you have multiple y chromosomes which are more or less closely related you can use the script in `scripts/find_homology.sh` to look for putative homology between their Y-STRs.

This script is based on the same `TRF` command as before and you can use the same flags `-m -M -e` exactly as in `typester.sh` to filter results. You can use this command with how many references you like as in:

```
find_homology.sh species_a_Y.fa [...] species_n_Y.fa
```

The script will produce three files:
1. `all_Y_strs.fa`: the fasta files of the filtered STRs found in your references. The sequences will be named as filename+an increasing integer based on their postion (i.e. given an input called Hsap_Y.fa the regions will be named Hsap_Y_1,2...n). Further, each sequence will have 200 bp additional base paires in the flanking regions (you can change this parameter with `-f`).
2. `regions_metadata.txt`: file with 5 columns: reference_id, start, end, length, name
3. `blast_all_v_all_Y_str.csv`: blast results of the all vs all search (outfmt 6).

You can check additional parameters with `./scripts/find_homology.sh -h`

If you are interested in homology with Human Y-STRs, you might use the Y_STR_panel.csv from [YLineageTracker](https://github.com/Shuhua-Group/Y-LineageTracker) with the `-p` option. You may obtain this file by running: `wget https://raw.githubusercontent.com/Shuhua-Group/Y-LineageTracker/main/LineageTracker/Data/Y_STR_panel.csv`. This is good as human STRs are well studied and have official names. If you use this you have to specify the fasta Human Y chromosome GRCh38 reference (CM000686.2 or NC_000024.10) with `-r`!

This script also allows to search for multicopies STR within a single reference!

Finally, you might visualize the blast results with:

`Rscript scripts/comparative_STR.R`

This will produce a pdf named `plots/blast_Y_STRs.pdf` with a pairwise comparison of blast results of all the references analyzed. Check `Rscript scripts/comparative_STR.R -h` for usage.

# TODOs

* add test data and output
* hipstr installation (docker?)
* discuss if adding snakemake pipeline or useless, I think it's kind of useless
* stutter model option????? *Depending on your project you may decide which kind of stutter model to use, a nice guide can be found [here](https://hipstr-tool.github.io/HipSTR/#in-depth-usage). In this case we will use the default stutter model (`--def-stutter-model`) as our samples are from different sequencing projects and have different coverages.*
* output filtered haplotypes (in which format?)
