from vcf_parser import VCFParser
import warnings
warnings.filterwarnings("ignore")
vcff = "mgp.v3.snps.rsIDdbSNPv137.vcf.gz"
my_parser = VCFParser(infile=vcff, split_variants=True, check_info=True)

varfile = open("mm10_simple_snp.tsv", "w")
varfile.write("chr\tbp\tREF\tALT\tC57BL6NJ\tBALBcJ\t129S1\tCASTEiJ\tDBA2J\tPWKPhJ\n")
for variant in my_parser:
	pos_info = variant['CHROM'] + "\t" +  variant['POS'] + "\t" + variant['REF'] + "\t" + variant['ALT']
	geno = variant['genotypes']
	
	# Pull SNPs and annotate with whether or not they have a reference allele
	C57BL6NJ = "\t" + str(int(str(geno['C57BL6NJ']) == "1/1"))
	BALBcJ = "\t" + str(int(str(geno['BALBcJ']) == "1/1"))
	X129S1 = "\t" + str(int(str(geno['129S1']) == "1/1"))
	CASTEiJ = "\t" + str(int(str(geno['CASTEiJ']) == "1/1"))
	DBA2J = "\t" + str(int(str(geno['DBA2J']) == "1/1"))
	PWKPhJ = "\t" + str(int(str(geno['PWKPhJ']) == "1/1"))
	outstr = pos_info + C57BL6NJ + BALBcJ + X129S1 + CASTEiJ + DBA2J + PWKPhJ
	varfile.write(outstr + "\n")
varfile.close()
