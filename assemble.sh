#!/bin/bash

wtdbg2 -x ont -g 200m -i ../mercatorum_wildtype.ont.fq.gz -t 24 -fo mercatorum_wildtype.ont.wtdbg2.asm1
minimap2 -t24 -ax map-ont -r2k mercatorum_wildtype.ont.wtdbg2.asm1.raw.fa ../mercatorum_wildtype.ont.fq.gz | samtools sort -@4 >mercatorum_wildtype.asm1.ont.bam && samtools view -F0x900 mercatorum_wildtype.asm1.ont.bam | wtpoa-cns -t 24 -d mercatorum_wildtype.ont.wtdbg2.asm1.raw.fa -i - -fo mercatorum_wildtype.ont.wtdbg2.asm1.cns.fa
