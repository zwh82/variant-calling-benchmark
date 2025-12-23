# variant-calling-benchmark


**Evaluation of variant detection for metagenomic long-read sequencing data (ONT R10.4). **

## benchmark

Although many variant calling tools already exist, no tools have been found that are suitable for both high-error-rate long reads and metagenomic data. Metagenomic data are highly complex, with multiple strains coexisting and uneven coverage depth.

```
mamba create -n vc
mamba activate vc
mamba install snakemake
```

```
cd workflow
snakemake -s variant_calling.smk --cores 5 --use-conda
```


## vcf eval

Having extensively utilized two prominent variant evaluation tools — `RTG vcfeval` and `hap.py (som.py)`—I've concluded that there remains no optimal solution for comprehensive variant assessment. Each tool exhibits distinct strengths and limitations that must be carefully considered for different applications.

### rtg vcfeval

+ Advantages

RTG excels at reconciling different representations of the same variant (indels or SVs) across variant calling tools. For instance, it correctly identifies that `G → GC` (GATK representation) and `GGT → GCGT` (Freebayes representation) represent the same underlying variant at identical genomic positions.

+ Limitations (I think)

ALT allele count insensitivity: RTG fails to distinguish between multi-allelic representations. For example, `G → A,C` and `G → A` are incorrectly considered identical variants, despite representing different allelic compositions.

GT field dependency: When format and sample columns are specified, RTG incorporates genotype (GT) information into matching decisions. This leads to counterintuitive mismatches—even when REF and ALT alleles are identical—if genotypes differ (e.g., truth=2/2 vs prediction=1/1).

### hap.py (som.py)
Stringent positional matching: hap.py requires exact matches of position, reference allele, and alternative allele(s). This simplistic approach fails to recognize equivalent variant representations that RTG successfully reconciles.


## The challenge and ideal vcf evaluation tools of metagenome
Standard variant calling tools typically predict genotypes as 0/0, 0/1, 1/1, or 0/2, this may be not suitable for metagenome. In microbial communities where multiple strains coexist, representing genetic diversity through binary diploid genotypes is biologically unrealistic .

A ideal evaluation tools:

+ performs variant comparison as RTG: can deal with differences in representation that commonly arise, particularly when dealing with complex variants or when comparing variants produced by different callers.

+ each ALT comparison: systematically evaluate each alternative allele independently without comparing GT.
