# drosophila mercatorum genome assemblies

This summarizes steps required to build and evaluate a drosophila mercatorum genome assembly.
Two assemblies are made using the same process, one for drosophila mercatorum wildtype, and another for a parthenogenic laboratory strain.
These are compared to each other, and to the drosophila melanogaster reference genome.
A coverage analysis shows that both assemblies are of high quality.

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

Then, with version `/gnu/store/afyq274dvnz71y7mmrsp53kx6r0m71fv-wfmash-0.9.2+1722e9e-1/bin/wfmash` of wfmash:

Align the two mercatorum strains to dmel6 reference.

```
wfmash -p 70 -t 48 dmel6.fa.gz dmerc_wildtype.fa.gz | gzip >dmel6_vs_dmerc_wildtype.26.paf.gz
wfmash -p 70 -t 48 dmel6.fa.gz dmerc_partho.fa.gz | gzip >dmel6_vs_dmerc_partho.paf.26.gz
```

The alignment between the two assemblies we made is easier because the divergence is much lower, so we leave default settings in wfmash:

```
wfmash -t 48 dmerc_wildtype.fa.gz dmerc_partho.fa.gz | gzip >dmerc_wildtype_vs_dmerc_partho.25.paf.gz
```

We compute summary statistics on top of the outputs of the alignments.

The two mercatorum strains have about 1.33% divergence between them, when using [gap compressed identity](http://lh3.github.io/2018/11/25/on-the-definition-of-sequence-identity).

```
-> % f=dmerc_wildtype_vs_dmerc_partho.25.paf.gz ; zcat $f | sed s/gi:f:// | awk '{ sum += $13 * $10; tot += $10; block += $11; } END { print "matches",tot; print "block",block; print "gap.id", sum/tot; print "block.id",tot/block * 100; }'
matches 142245925
block 151145805
gap.id 98.6674
block.id 94.1117
```

In contrast, between the wildtype and melanogaster, we see 24.2% divergence.

```
-> % f=dmel6_vs_dmerc_wildtype.26.paf.gz ; zcat $f | sed s/gi:f:// | awk '{ sum += $13 * $10; tot += $10; block += $11; } END { print "matches",tot; print "block",block; print "gap.id", sum/tot; print "block.id",tot/block * 100; }'
matches 17679988
block 75282428
gap.id 75.8478
block.id 23.4849
```

And between the parthenote and melanogaster, we see 24.5% divergence.

```
-> % f=dmel6_vs_dmerc_partho.26.paf.gz ; zcat $f | sed s/gi:f:// | awk '{ sum += $13 * $10; tot += $10; block += $11; } END { print "matches",tot; print "block",block; print "gap.id", sum/tot; print "block.id",tot/block * 100; }'
matches 17696931
block 73444426
gap.id 75.548
block.id 24.0957
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


## self heterozygosity analysis

To establish estimates for pairwise heterozygosity in each strain, we aligned the illumina data used for polishing back to the assemblies.

```
bwa mem -t 24 assemblies/dmerc_wildtype.fa.gz mercatorum_wildtype/illumina/SLX-18528.NEBNext04.000000000-K3MKR.s_1.r_1.fq.gz mercatorum_wildtype/illumina/SLX-18528.NEBNext04.000000000-K3MKR.s_1.r_2.fq.gz >dmerc_wildtype.illumina.sam && samtools view  -b dmerc_wildtype.illumina.sam >dmerc_wildtype.illumina.bam && sambamba sort -t 24 dmerc_wildtype.illumina.bam
bwa mem -t 24 assemblies/dmerc_partho.fa.gz mercatorum_parthenote/illumina/SLX-18526.DNAA001.HGYGGDRXX.s_2.r_1.fq.gz mercatorum_parthenote/illumina/SLX-18526.DNAA001.HGYGGDRXX.s_2.r_2.fq.gz >dmerc_partho.illumina.sam && samtools view  -b dmerc_partho.illumina.sam >dmerc_partho.illumina.bam && sambamba sort -t 24 dmerc_partho.illumina.bam
```

We then call variants.

```
freebayes -f assemblies/dmerc_wildtype.fa dmerc_wildtype.illumina.sorted.bam >dmerc_wildtype.illumina.vcf
freebayes -f assemblies/dmerc_partho.fa dmerc_partho.illumina.sorted.bam >dmerc_partho.illumina.vcf
```

Normalize them using vcfwave.

```
<dmerc_partho.illumina.vcf awk '/^#/ || NF == 10' | vcffilter -f 'QUAL > 20' | vcfwave | vcffilter -f 'TYPE = snp' >dmerc_partho.illumina.q20wave.vcf
<dmerc_wildtype.illumina.vcf awk '/^#/ || NF == 10' | vcffilter -f 'QUAL > 20' | vcfwave | vcffilter -f 'TYPE = snp' >dmerc_wildtype.illumina.q20wave.vcf
```

And measure heterozygosity using bcftools.

```
-> % bcftools stats -s - dmerc_partho.illumina.q20wave.vcf | grep ^PSC | cut -f 6
16474
```

That's 16k heterozygous SNPs in the entire parthenote genome.
The genome is 161570079bp, suggesting a pairwise heterozygosity of 0.0102% _for SNPs_.

In contrast, the wildtype genome has higher heterozygosity.

```
-> % bcftools stats -s - dmerc_wildtype.illumina.q20wave.vcf | grep ^PSC | cut -f 6
109035
```

The genome is 171182504bp, suggesting a pairwise heterozygosity of 0.0637% _for SNPs_.

These are likely underestimates, but the basis of the assemblies in ONT means that we can't easily repeat this analysis for indels, which will have a much higher error rate.

## inter-strain divergence

To establish if the estimates for inter-strain divergence of around 1% that we established by alignment of the genomes were correct, we then repeated the above analysis, but for each strain we aligned it against the other mercatorum assembly, and called variants.

```
bwa mem -t 48 assemblies/dmerc_wildtype.fa.gz mercatorum_parthenote/illumina/SLX-18526.DNAA001.HGYGGDRXX.s_2.r_1.fq.gz mercatorum_parthenote/illumina/SLX-18526.DNAA001.HGYGGDRXX.s_2.r_2.fq.gz >dmerc_partho.vs_wildtype.illumina.sam && samtools view  -b dmerc_partho.vs_wildtype.illumina.sam >dmerc_partho.vs_wildtype.illumina.bam && sambamba sort -t 48 dmerc_partho.vs_wildtype.illumina.bam && freebayes -f assemblies/dmerc_wildtype.fa dmerc_partho.vs_wildtype.illumina.sorted.bam >dmerc_partho.vs_wildtype.illumina.vcf
bwa mem -t 48 assemblies/dmerc_partho.fa.gz mercatorum_wildtype/illumina/SLX-18528.NEBNext04.000000000-K3MKR.s_1.r_1.fq.gz mercatorum_wildtype/illumina/SLX-18528.NEBNext04.000000000-K3MKR.s_1.r_2.fq.gz >dmerc_wildtype.vs_partho.illumina.sam && samtools view  -b dmerc_wildtype.vs_partho.illumina.sam >dmerc_wildtype.vs_partho.illumina.bam && sambamba sort -t 48 dmerc_wildtype.vs_partho.illumina.bam && freebayes -f assemblies/dmerc_partho.fa dmerc_wildtype.vs_partho.illumina.sorted.bam >dmerc_wildtype.vs_partho.illumina.vcf
```

We normalized these as before.

```
<dmerc_wildtype.vs_partho.illumina.vcf awk '/^#/ || NF == 10' | vcffilter -f 'QUAL > 20' | vcfwave | vcffilter -f 'TYPE = snp' >dmerc_wildtype.vs_partho.illumina.q20wave.vcf
<dmerc_partho.vs_wildtype.illumina.vcf awk '/^#/ || NF == 10' | vcffilter -f 'QUAL > 20' | vcfwave | vcffilter -f 'TYPE = snp' >dmerc_partho.vs_wildtype.illumina.q20wave.vcf
```

And measured the number of non-reference homozogyous SNPs (not hets) for the parthenote:

```
-> % bcftools stats -s - dmerc_partho.vs_wildtype.illumina.q20wave.vcf | grep ^PSC | cut -f 5
1442703
```

This suggests a strain-wise divergence of 0.84%.

```
-> % bcftools stats -s - dmerc_wildtype.vs_partho.illumina.q20wave.vcf | grep ^PSC | cut -f 5  
1328517
```

We see the same for the other direction, with an estimate of 0.82%.

Due to our focus on SNPs, and limitation to non-reference homozygous cases, these should be strict underestimates of the pairwise divergence.
However, they support the idea that the mercatorum strains differ by around 1%.

## merqury kmer spectrum analysis

We have estimated a very low heterozygosity level for both mercatorum strains.
To further evaluate this, we explored the kmer copy spectrum of the assemblies relative to their respective illumina read sets using [merqury](https://github.com/marbl/merqury).

<img src="https://github.com/ekg/drosophila/raw/main/dmerc_partho.merqury.dmerc_partho.spectra-cn.hist.png">

<img src="https://github.com/ekg/drosophila/raw/main/dmerc_wildtype.merqury.dmerc_wildtype.spectra-cn.hist.png">