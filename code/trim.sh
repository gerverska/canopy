#!/bin/bash

# remove 5` gene primers and degenerate heterogeneity spacers from forward and reverse paired-end sequences
# using cutadapt v1.18 https://github.com/marcelm/cutadapt and SeqPurge 2019_05 (part of https://github.com/imgag/ngs-bits)

# define gene primers for trimming
fwd_primer="CCTCCGCTTATTGATATGCTTAART" # ITS4-fun
rc_fwd_primer="AYTTAAGCATATCAATAAGCGGAGG" # ITS4-fun, reverse complement
rev_primer="AACTTTYRRCAAYGGATCWCT" #5.8S-fun
rc_rev_primer="AGWGATCCRTTGYYRAAAGTT" # 5.8S-fun, reverse complement

# make directory for scratch data
scratch="output/scratch/"
mkdir -p $scratch

# make directory for trimming output
out="output/trim/"
mkdir -p $out

# make file for trimming summary output
log="${out}canopy.summary.tab"
echo -e "sample\traw\ttrim.1\ttrim.2" > ${log}

# loop over sample
for fwd in $( find output/demux/ -name "*R1.fastq.gz" | grep -v "undetermined"); do

    rev=`echo $fwd | sed 's/R1.fastq.gz/R2.fastq.gz/'`
    samp=`echo $fwd | awk -F'/' '{print $3}'| awk -F'.' '{print $1}'`
    
    # trim primers in paired-end mode
    fwd_trim1="${scratch}${samp}.R1.trim1.fq.gz"
    rev_trim1="${scratch}${samp}.R2.trim1.fq.gz"
    
    # important cutadapt setting: overlap = length of shorter primer, min/max = set to global expectations, e = error rate
    cutadapt --quiet -g $fwd_primer -G $rev_primer --discard-untrimmed --overlap 21 -e 0.2 --minimum-length 270 --maximum-length 277 -o $fwd_trim1 -p $rev_trim1 $fwd $rev
    
    # count seqs after first trim step (the "@M01498" should match the first 5+ characters of your fastq headers)
    if [[ -f $fwd_trim1 && -f $rev_trim1 ]]; then
        trim1sum=`gzip -cd $fwd_trim1 | grep -c '^@M01498'`
    else
        trim1sum=0
    fi
    
    # trim any remaining primer overhang
    if [ $trim1sum -gt 0 ]; then
        fwd_trim2="${scratch}${samp}.R1.trim2.fq.gz"
        rev_trim2="${scratch}${samp}.R2.trim2.fq.gz"
        SeqPurge -in1 $fwd_trim1 -in2 $rev_trim1 -out1 $fwd_trim2 -out2 $rev_trim2 -a1 $rc_rev_primer -a2 $rc_fwd_primer -qcut 0 -ncut 0 -min_len 150 -summary ${scratch}${samp}.log.txt
    fi
    
    # count seqs after final step (the "@M01498" should match the first 5+ characters of your fastq headers)
    if [[ -f $fwd_trim2 && -f $rev_trim2 ]]; then
        trim2sum=`gzip -cd $fwd_trim2 | grep -c '^@M01498'`
    else
        trim2sum=0
    fi
    
    # keep only samples that made it through the processing (also truncate reads)
    if [ $trim2sum -gt 0 ]; then
        fwdout="${out}${samp}.R1.fq.gz"
        revout="${out}${samp}.R2.fq.gz"
        cutadapt --quiet -g XX -G XX --length 270 -o $fwdout -p $revout $fwd_trim2 $rev_trim2
    fi
    
    raw=`gzip -cd $fwd | grep -c '^@M01498'`
    echo -e $samp"\t"$raw"\t"$trim1sum"\t"$trim2sum | cat >> ${log}
    #remove tmp files
    rm ${scratch}${samp}.*
done

rm -r ${scratch}