### Testing hisnpper

```

# Pick applicable filepath depending on the system
fasta="/data/aryee/caleb/mESC_allele/mm10-masked/genome.fa"
fasta="/Users/clareau/dat/genomes/mm10/mm10-Nmask.fa"

hisnpper ase -b mESC_test.bam -s CAST_129S_testSNPtable.tsv.gz -o out --fasta $fasta -bt XB
hisnpper haplotype -b mESC_test.bam -s CAST_129S_testSNPtable.tsv.gz -o out --fasta $fasta -bt XB -z
hisnpper-edits -b mESC_test.bam -o edits --fasta $fasta -e N_T -z

```

<br><br>
