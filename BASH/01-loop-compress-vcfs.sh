

for i in $1/*.vcf; 
do bgzip -c $i > $i.gz
done
