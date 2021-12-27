# mercatorum wildtype assembly and polishing

First we generate the assembly from the ONT data using wtdbg2.

```
mkdir asm1
cd asm1
wtdbg2 -x ont -g 200m -i ../mercatorum_wildtype.ont.fq.gz -t 24 -fo mercatorum_wildtype.ont.wtdbg2.asm1
wtpoa-cns -t 24 -i mercatorum_wildtype.ont.wtdbg2.asm1.ctg.lay.gz -fo mercatorum_wildtype.ont.wtdbg2.asm1.raw.fa
minimap2 -t24 -ax map-ont -r2k mercatorum_wildtype.ont.wtdbg2.asm1.raw.fa ../mercatorum_wildtype.ont.fq.gz \
    | samtools sort -@4 >mercatorum_wildtype.asm1.ont.bam && samtools view -F0x900 mercatorum_wildtype.asm1.ont.bam \
    | wtpoa-cns -t 24 -d mercatorum_wildtype.ont.wtdbg2.asm1.raw.fa -i - -fo mercatorum_wildtype.ont.wtdbg2.asm1.cns.fa
```

Now, we'll polish using some illumina data from the same isolate.

```
# first the alignment
bwa mem -t 48 asm1/mercatorum_wildtype.ont.wtdbg2.asm1.cns.fa mercatorum_wildtype/illumina/SLX-18528.NEBNext04.000000000-K3MKR.s_1.r_1.fq.gz mercatorum_wildtype/illumina/SLX-18528.NEBNext04.000000000-K3MKR.s_1.r_1.fq.gz >asm1/mercatorum_wildtype.illumina.sam
samtools view -b asm1/mercatorum_wildtype.illumina.sam >asm1/mercatorum_wildtype.illumina.raw.bam
sambamba sort -t 48 -o asm1/mercatorum_wildtype.illumina.bam asm1/mercatorum_wildtype.illumina.raw.bam
rm asm1/mercatorum_wildtype.illumina.sam asm1/mercatorum_wildtype.illumina.raw.bam
# then the variant calling
cat asm1/mercatorum_wildtype.ont.wtdbg2.asm1.cns.fa.fai | cut -f 1 | parallel 'freebayes -f asm1/mercatorum_wildtype.ont.wtdbg2.asm1.cns.fa --region {} asm1/mercatorum_wildtype.illumina.raw.bam | bcftools view --no-version -Ob -o - >calls/{}.polish1.bcf'
# finally concatenating the variant calls
cat asm1/mercatorum_wildtype.ont.wtdbg2.asm1.cns.fa.fai | cut -f 1 | while read f; do echo calls/$f.polish1.bcf; done >concat_list.txt
bcftools concat -nf concat_list.txt | bcftools view -Ou -e'type="ref"' | bcftools norm -Ob -f asm1/mercatorum_wildtype.ont.wtdbg2.asm1.cns.fa -o asm1/mercatorum_wildtype.illumina.bcf
# and applying the polishing
bcftools consensus -i'QUAL>1 && (GT="AA" || GT="Aa")' -Hla -f asm1/mercatorum_wildtype.ont.wtdbg2.asm1.cns.fa asm1/mercatorum_wildtype.illumina.bcf >asm1/mercatorum_wildtype.ont.wtdbg2.asm1.cns.polish.fa
```