# follow https://fastly.jsdelivr.net/gh/RealTimeGenomics/rtg-tools@master/installer/resources/tools/RTGOperationsManual.pdf

sdf=../test/reference/ref.sdf

# no squash-ploidy
vcf=pred_alt2.vcf.gz
truth=truth_alt2.vcf.gz
output=eval_alt2_no_squash-ploidy
if [ ! -d $output ]; then 
    rtg vcfeval -c $vcf -b $truth -o $output -t $sdf 
fi 

# alt2
vcf=pred_alt2.vcf.gz
truth=truth_alt2.vcf.gz
output=eval_alt2
if [ ! -d $output ]; then 
    rtg vcfeval -c $vcf -b $truth -o $output -t $sdf --squash-ploidy
fi 

# alt2 predict missing1
vcf=pred_alt2_missing1.vcf.gz
truth=truth_alt2.vcf.gz
output=eval_alt2_missing1
if [ ! -d $output ]; then 
    rtg vcfeval -c $vcf -b $truth -o $output -t $sdf --squash-ploidy
fi 

# truth_alt2_GT2-2
vcf=pred_alt2.vcf.gz
truth=truth_alt2_GT2-2.vcf.gz
output=eval_alt2_truth_GT2-2
if [ ! -d $output ]; then 
    rtg vcfeval -c $vcf -b $truth -o $output -t $sdf --squash-ploidy
fi 

vcf=pred_alt2.vcf.gz
truth=truth_alt2_GT0-1-2.vcf.gz
output=eval_alt2_truth_GT0-1-2
if [ ! -d $output ]; then 
    rtg vcfeval -c $vcf -b $truth -o $output -t $sdf --squash-ploidy
fi 

# no sample
vcf=pred_alt2.vcf.gz
truth=truth_alt2_no_sample.vcf.gz
output=eval_alt2_truth_no_example
if [ ! -d $output ]; then 
    rtg vcfeval -c $vcf -b $truth -o $output -t $sdf --squash-ploidy --sample ALT
fi 

# alt2 predict missing1, turth no sample
vcf=pred_alt2_missing1_no_sample.vcf.gz
truth=truth_alt2_no_sample.vcf.gz
output=eval_alt2_missing1_truth_no_example
if [ ! -d $output ]; then 
    rtg vcfeval -c $vcf -b $truth -o $output -t $sdf --squash-ploidy --sample ALT
fi 