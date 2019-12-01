### Testing hisnpper

```

# Pick applicable filepath depending on the system
fasta="/data/aryee/caleb/mESC_allele/mm10-masked/genome.fa"
fasta="/Users/clareau/dat/genomes/GRCm38/GRCm38_hardmasked.fa"

hisnpper ase -b ds.bam -s CAST_129S_testSNPtable.tsv.gz -o out --fasta $fasta -bt XB
hisnpper haplotype -b ds.bam -s CAST_129S_testSNPtable.tsv.gz -o out --fasta $fasta -bt XB -z
hisnpper-edits -b ds.bam -o edits --fasta $fasta -e N_T -z

```

<br><br>
