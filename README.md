
Filter VCF file to identify germline de novo mutations.

The VCF should contain a mother, a father and their multiple offspring.

For a de novo mutation:
	
	The child with denovo should be heterozygous(0/1) for autosomes and homozygous reference(1/1) for male X.
	
	The denovo in the child should have supporing reads for both forward and reverse strans(ADF>0, ADR>0).

	The child with denovo should meet  `--allele-balance`, `--min-dp`, `--max-dp`, `--min-gq`. (for male X, depth halved and not check allelic balance).	

	The parents are both homozygous reference(0/0), no reads contain the mutation allele(AD alt = 0);

	The parents shoule meet  `--min-dp`, `--max-dp`, `--min-gq`.
	
	The denovo should not appear in other offspring genotype.

	The number of impurity samples should no greater than `--max-impurity-sample`. The impurity sample defined as other offspring samples(not the denovo offspring)
whose number of reads supporting the denovo larger than or equal `--ad-min-impurity` 




# Install

Require: BOOST library

```bash
cd src 
g++ --std=c++11 main.cpp -o findDeNovo -lboost_iostreams

```

# Input

- VCF file(gz or text format), contains a paternal sample, a maternal sample and multiple offspring samples;
	assume having GT, GQ, DP, ADF and ADR annotation for each sample.

- PED file, no header and no comment lines. 
	For ped file: FamilyID, ChildID, FatherID, MotherID, Sex, Phenotype
	Here we only use column: ChildID, FatherID, MotherID and Sex.     

# Output

```
CHROM  POS      REF  ALT  CHILD      CHILD.DP  CHILD.AD  CHILD.GQ  MOTHER     MOTHER.DP  MOTHER.AD  MOTHER.GQ  FATHER     FATHER.DP  FATHER.AD  FATHER.GQ
2L     1886521  C    T    sample_84  52        43,9      99        sample_32  42         42,0       99         sample_43  31         31,0       99
2L     1934106  G    C    sample_88  12        4,8       99        sample_32  30         30,0       97         sample_43  36         36,0       99
2L     1934106  G    T    sample_74  19        10,9      99        sample_32  30         30,0       97         sample_43  36         36,0       99
```



# Usage:

```
Options: 
        --vcf                 Input VCF file or gz VCF file, assume having GT,GQ,DP,AD,ADF,ADR
        --ped                 Input PED file, not comment lines allowed
        --out                 Output file.
        -h,--help             Show help message
        --min-dp              The minimum read depth. int [0,Inf); default 10
        --max-dp              The maximum read depth. int [0,Inf); default 150
        --min-gq              The minimum genotype quality. int [0,Inf); default 70
        --allele-balance      The allele balance for a heterozygote
                                allele balance less than this value or larger than 1-this value
                                will be filtered out. double [0,0.5); default 0.15
        --ad-min-impurity     For homoRef, the number of ALT read above which will be considered
                                as impurity sample. int [1,Inf); default 1 (most strigent)
        --max-impurity-sample The maximum number of impurity samples allowed, above which the
                                SNP will not be considered for mutations. int [0, Inf); default 0 (most strigent)
```


# Example

```bash
./findDeNovo --vcf ../testData/example.vcf.gz --ped ../testData/samples.ped --out example.out

```
