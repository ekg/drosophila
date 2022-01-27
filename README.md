# mercatorum wildtype assembly and polishing

## assembly

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

This same process was repeated for both the wildtype and parthenogenic strains.

## alignments

We wanted to understand the relationship between the assemblies (`dmerc_wildtype` and `dmerc_partho`), and between both and the drosophila melanogaster reference genome (`dmel6`).

And make alignments:

```
fastix -p "dmerc_partho#" dmerc.ctg.polished1.fa | bgzip -@ 24 >dmerc_partho.fa.gz && samtools faidx dmerc_partho.fa.gz
fastix -p "dmerc_wildtype#" mercatorum_wildtype.ont.wtdbg2.asm1.cns.polish.fa | bgzip -@ 24 >dmerc_wildtype.fa.gz && samtools faidx dmerc_wildtype.fa.gz
fastix -p "dmel6#" GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna | bgzip -@ 24 >dmel6.fa.gz && samtools faidx dmel6.fa.gz
```

Then, with version `/gnu/store/d6ajy0ajxkbb22gkgi32g5c4jy8rs8bv-wfmash-0.6.0+948f168-21/bin/wfmash` of wfmash:

Align the two mercatorum strains to dmel6 reference.

```
wfmash -p 70 -w 17 -s 5000 -l 0 -O -t 48 dmel6.fa.gz dmerc_partho.fa.gz >dmel6_vs_dmerc_partho.17.paf

wfmash -p 70 -w 17 -s 5000 -l 0 -O -t 48 dmel6.fa.gz dmerc_wildtype.fa.gz >dmel6_vs_dmerc_wildtype.17.paf
```

The alignment between the two assemblies we made is easier because the divergence is much lower, so the settings are different:

```
wfmash -p 95 -t 48 dmerc_wildtype.fa.gz dmerc_partho.fa.gz >dmerc_wildtype_vs_dmerc_partho.p95.paf
```

The average identities shown here, and the total length in the alignment:

```
-> % < dmel6_vs_dmerc_partho.17.paf | sed 's/gi:f://' | awk '$13 > 0.55' | awk '{ sum += $13 * $11; tot += $11; } END { print sum/tot; print tot; }'
0.63003
50242551

-> % < dmel6_vs_dmerc_wildtype.17.paf | sed 's/gi:f://' | awk '$13 > 0.55' | awk '{ sum += $13 * $11; tot += $11; } END { print sum/tot; print tot; }'
0.63238
50207983

-> % < dmerc_wildtype_vs_dmerc_partho.p95.paf | sed 's/gi:f://' | awk '$13 > 0.55' | awk '{ sum += $13 * $11; tot += $11; } END { print sum/tot; print tot; }'
0.98591
146250790
```

## coverage analysis

Finally we completed a coverage analysis using the input ONT reads to measure the completeness of the assemblies.

First with mapping:

```
minimap2 -x map-ont -a -t 48 assemblies/dmerc_wildtype.fa.gz mercatorum_wildtype.ont.fq.gz >mercatorum_wildtype.ont.sam && samtools view -b mercatorum_wildtype.ont.sam >mercatorum_wildtype.ont.bam.1 && samtools sort -@ 24 mercatorum_wildtype.ont.bam.1 >mercatorum_wildtype.ont.bam && samtools index mercatorum_wildtype.ont.bam && rm -f mercatorum_wildtype.ont.bam.1 mercatorum_wildtype.ont.sam
minimap2 -x map-ont -a -t 48 assemblies/dmerc_partho.fa.gz mercatorum_partho.ont.fq.gz >mercatorum_partho.ont.sam && samtools view -b mercatorum_partho.ont.sam >mercatorum_partho.ont.bam.1 && samtools sort -@ 24 mercatorum_partho.ont.bam.1 >mercatorum_partho.ont.bam && samtools index mercatorum_partho.ont.bam && rm -f mercatorum_partho.ont.bam.1 mercatorum_partho.ont.sam
```

Then `samtools coverage` to get information about coverage of the contigs:

```
samtools coverage mercatorum_partho.ont.bam >mercatorum_partho.ont.bam.samtools_coverage.txt
samtools coverage mercatorum_wildtype.ont.bam >mercatorum_wildtype.ont.bam.samtools_coverage.txt
samtools coverage -A -w 160 mercatorum_partho.ont.bam >mercatorum_partho.ont.bam.samtools_coverage.asciiart.txt
samtools coverage -A -w 160 mercatorum_wildtype.ont.bam >mercatorum_wildtype.ont.bam.samtools_coverage.asciiart.txt
```

To obtain a quick summary of the coverage statistics, weighted on a per-basepair (not per-contig) level, we used this awk script.

These show the fraction of covered bases, which suggests a very low rate of gaps in the assemblies:

```
< mercatorum_wildtype.ont.bam.samtools_coverage.txt awk '{ len=$3-$2; tot += len; sum += len * $6; } END { print sum / tot; }'
96.9121
< mercatorum_partho.ont.bam.samtools_coverage.txt awk '{ len=$3-$2; tot += len; sum += len * $6; } END { print sum / tot; }'  
99.398
```

In the wildtype, 96.9% of the bases in the assembly are covered. (n.b. the wildtype should be more heterozygous.)
In the parthenogenic strain, 99.4% of the bases are covered.

And the depth of coverage is high in both assemblies:

```
< mercatorum_wildtype.ont.bam.samtools_coverage.txt awk '{ len=$3-$2; tot += len; sum += len * $7; } END { print sum / tot; }'
79.2321
< mercatorum_partho.ont.bam.samtools_coverage.txt awk '{ len=$3-$2; tot += len; sum += len * $7; } END { print sum / tot; }'
93.8656
```

So 79.2X on average in the wildtype and 93.9X in the parthenote.