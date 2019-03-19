# hisnpper
Toolkit for extract reads overlapping SNPs ++ annotating bam files with haplotype SNVs


### Usage

```
hisnpper --help

Usage: hisnpper [OPTIONS] [ase|haplotype]

  hisnpper: Infer allele-specific attributes from sequence data

  Caleb Lareau, clareau <at> broadinstitute <dot> org

Options:
  --version                  Show the version and exit.
  -b, --bamfile TEXT         Bam file that will be parsed / annotated.
  -s, --snps TEXT            Input table containing SNP information.
  -f, --fasta TEXT           Input fasta file; needs to be hard-masked at
                             relevant base pairs
  -o, --output TEXT          Output directory for analysis; see documentation.
  -n, --name TEXT            Name for all of the output files (default: uses
                             the .bam prefix)
  -bt, --barcode-tag TEXT    Two letter tag that indicates the single-cell ID
  -ht, --haplotype-tag TEXT  Two letter tag for bam file for phased reads
  -ma, --min-aq INTEGER      Minimum alignment quality for read to be
                             considered.
  -c, --ncores INTEGER       Number of cores to be used in analysis
  -z, --keep-temp-files      Keep all intermediate files.
  --help                     Show this message and exit.

```

Description of essential parameters

#### mode
This parameter is either `ase` or `haplotype` by design. `haplotype` does everything `ase` does,
and more. One should only use `haplotype` if SNPs are phased. 

##### -b, --bamfile
Input bam file from CellRanger or bap. Nothing too crazy here. 

#### -f, fasta
**VERY IMPORTANT** The filepath supplied here must point to a fasta with
the SNPs (used in the parameter below) HARD-MASKED in the reference fasta. 

The backbone of the algorithm is to identify single bases with mismatches. 
If the reference genome isn't hard masked, then anything that matches the 
reference will not be reported. 

#### -s, --snps
This should be a four-column table with 

```
head CAST_129_finalSNPtable.chr.tsv 
chr	bp	X129S1	CASTEiJ
chr1	3000258	G	T
chr1	3001579	T	C
chr1	3003561	A	T
chr1	3005006	T	A
chr1	3007334	G	A
chr1	3007431	G	A
chr1	3007757	C	T
```

Note that for `ase`, only the first two columns matter... so this could technically be either a 2 column 
file or any other number of columns... only the first two will be considered

#### -o, --output 
Folder filepath for the output

#### -n, --name
Name of the output files to follow. By default, will utilize the base name of the corresponding .bam file.

#### --barcode-tag
Currently a required parameter; the expectation is that this is the single-cell barcode

#### --haplotype-tag
Only needed in `haplotype` mode. Two letters to specify the new sam tag associated with the assigned haplotype

#### --ncores
Necessary to set; defaults to 2 otherwise

<br><br>
