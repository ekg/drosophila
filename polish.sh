#!/bin/bash

# exit when any command fails
set -eo pipefail
# trace our commands
set -o xtrace

# first the alignment
#bwa index asm1/mercatorum_wildtype.ont.wtdbg2.asm1.cns.fa
#bwa mem -t 48 asm1/mercatorum_wildtype.ont.wtdbg2.asm1.cns.fa \
#    mercatorum_wildtype/illumina/SLX-18528.NEBNext04.000000000-K3MKR.s_1.r_1.fq.gz \
#    mercatorum_wildtype/illumina/SLX-18528.NEBNext04.000000000-K3MKR.s_1.r_2.fq.gz \
#    >asm1/mercatorum_wildtype.illumina.sam
#samtools view -b asm1/mercatorum_wildtype.illumina.sam >asm1/mercatorum_wildtype.illumina.raw.bam
#sambamba sort -t 48 -o asm1/mercatorum_wildtype.illumina.bam asm1/mercatorum_wildtype.illumina.raw.bam
#rm asm1/mercatorum_wildtype.illumina.sam asm1/mercatorum_wildtype.illumina.raw.bam

# then the variant calling
#mkdir -p calls
#cat asm1/mercatorum_wildtype.ont.wtdbg2.asm1.cns.fa.fai | cut -f 1 | parallel 'freebayes -f asm1/mercatorum_wildtype.ont.wtdbg2.asm1.cns.fa --region {} asm1/mercatorum_wildtype.illumina.bam | bcftools view --no-version -Ob -o - >calls/{}.polish1.bcf'

# finally concatenating the variant calls
#cat asm1/mercatorum_wildtype.ont.wtdbg2.asm1.cns.fa.fai | cut -f 1 | while read f; do echo calls/$f.polish1.bcf; done >concat_list.txt
#bcftools concat -nf concat_list.txt | bcftools view -Ou -e'type="ref"' | bcftools norm -Ob -f asm1/mercatorum_wildtype.ont.wtdbg2.asm1.cns.fa -o asm1/mercatorum_wildtype.illumina.bcf
bcftools index asm1/mercatorum_wildtype.illumina.bcf

# and applying the polishing
bcftools consensus -i'QUAL>1 && (GT="AA" || GT="Aa")' -Hla -f asm1/mercatorum_wildtype.ont.wtdbg2.asm1.cns.fa asm1/mercatorum_wildtype.illumina.bcf >asm1/mercatorum_wildtype.ont.wtdbg2.asm1.cns.polish.fa
