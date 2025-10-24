#!/bin/bash


#Input files
#set id, file(perbase), file(depth) from coverageData

#Set variables
#strain id
id=$1
echo "id=$id"

#reference genome
reference=$2
echo "reference=$reference"

#baseDir for snpEff
baseDir=$3
echo "baseDir=$baseDir"

#annotation genome for snpeff
snpeff=$4
echo "snpEff=$snpeff"

fast=$5
echo "fast = $fast"


#import GATK filtering parameters
source ${baseDir}/configs/gatk_source.config

# Prefer running snpEff from the isolated Java-21 env `ardap_snpeff21`.
# If `SNPEFF_CMD` is provided in the environment it will be used (useful for
# containers). Otherwise prefer `mamba run -n ardap_snpeff21 snpEff`, fall back
# to `conda run -n ardap_snpeff21 snpEff`, and finally to plain `snpEff`.
if [ -z "${SNPEFF_CMD+x}" ]; then
  if command -v mamba >/dev/null 2>&1; then
    SNPEFF_CMD="mamba run -n ardap_snpeff21 snpEff"
  elif command -v conda >/dev/null 2>&1; then
    SNPEFF_CMD="conda run -n ardap_snpeff21 snpEff"
  else
    SNPEFF_CMD="snpEff"
  fi
fi

if [ "$fast" == "true" ]; then

  #test for interval file
  if [ ! -s ${baseDir}/Databases/${snpeff}/intervals.list ]; then


cat << _EOF_ > interval.coverage.query.txt
SELECT 
	Coverage.Gene
FROM 
	Coverage
_EOF_

cat << _EOF_ > interval.query.txt
SELECT 
	Variants_SNP_indel.Gene_name
FROM 
	Variants_SNP_indel
_EOF_

  #create interval file
  sqlite3 ${baseDir}/Databases/${snpeff}/${snpeff}.db < interval.query.txt > ${id}.gene.list
  sqlite3 ${baseDir}/Databases/${snpeff}/${snpeff}.db < interval.coverage.query.txt > ${id}.gene.list2

cat ${id}.gene.list ${id}.gene.list2 interval.coverage.query.txt > ${id}.gene.list.tmp
mv ${id}.gene.list.tmp ${id}.gene.list

  #create interval file

  uniq ${id}.gene.list > ${id}.gene.list.tmp
  mv ${id}.gene.list.tmp ${id}.gene.list
  ${SNPEFF_CMD} genes2bed ${snpeff} -dataDir ${baseDir}/resources/snpeff -f ${id}.gene.list > ${id}.intervals.list
  sed -i 's/\t/:/' ${id}.intervals.list
  sed -i 's/\t/-/' ${id}.intervals.list
  awk '{print $1}' ${id}.intervals.list | tail -n+2 > ${id}.intervals.list.tmp
  mv ${id}.intervals.list.tmp ${id}.intervals.list
  gatk HaplotypeCaller -R ${reference} --ploidy 1 --I ${id}.dedup.bam -O ${id}.raw.snps.indels.vcf -L ${id}.intervals.list --interval-padding 200
 else 
  gatk HaplotypeCaller -R ${reference} --ploidy 1 --I ${id}.dedup.bam -O ${id}.raw.snps.indels.vcf -L ${baseDir}/Databases/${snpeff}/intervals.list --interval-padding 200
 fi
else
#echo "intervals == $intervals"
  gatk HaplotypeCaller -R ${reference} --ploidy 1 --I ${id}.dedup.bam -O ${id}.raw.snps.indels.vcf 
fi



gatk SelectVariants -R ${reference} -V ${id}.raw.snps.indels.vcf -O ${id}.raw.snps.vcf -select-type SNP
gatk SelectVariants -R ${reference} -V ${id}.raw.snps.indels.vcf -O ${id}.raw.indels.vcf -select-type INDEL

gatk VariantFiltration -R ${reference} -O ${id}.filtered.snps.vcf -V ${id}.raw.snps.vcf \
-filter "MLEAF < ${MLEAF_SNP}" --filter-name "AFFilter" \
-filter "QD < ${QD_SNP}" --filter-name "QDFilter" \
-filter "MQ < ${MQ_SNP}" --filter-name "MQFilter" \
-filter "FS > ${FS_SNP}" --filter-name "FSFilter" \
-filter "QUAL < ${QUAL_SNP}" --filter-name "StandardFilters"

  header=`grep -a -n "#CHROM" ${id}.filtered.snps.vcf | cut -d':' -f 1`
  head -n $header ${id}.filtered.snps.vcf > ${id}.snp_head
  cat ${id}.filtered.snps.vcf | grep PASS | cat ${id}.snp_head - > ${id}.PASS.snps.vcf

gatk VariantFiltration -R ${reference} -O ${id}.failed.snps.vcf -V ${id}.raw.snps.vcf \
-filter "MLEAF < ${MLEAF_SNP}" --filter-name "FAIL" \
-filter "QD < ${QD_SNP}" --filter-name "FAIL1" \
-filter "MQ < ${MQ_SNP}" --filter-name "FAIL2" \
-filter "FS > ${FS_SNP}" --filter-name "FAIL3" \
-filter "QUAL < ${QUAL_SNP}" --filter-name "FAIL5"

  header=`grep -a -n "#CHROM" ${id}.failed.snps.vcf | cut -d':' -f 1`
  head -n $header ${id}.failed.snps.vcf > ${id}.snp_head
  cat ${id}.failed.snps.vcf | grep FAIL | cat ${id}.snp_head - > ${id}.FAIL.snps.vcf

gatk VariantFiltration -R ${reference} -O ${id}.filtered.indels.vcf -V ${id}.raw.indels.vcf \
-filter "MLEAF < ${MLEAF_INDEL}" --filter-name "AFFilter" \
-filter "QD < ${QD_INDEL}" --filter-name "QDFilter" \
-filter "FS > ${FS_INDEL}" --filter-name "FSFilter" \
-filter "QUAL < ${QUAL_INDEL}" --filter-name "QualFilter"

  header=`grep -a -n "#CHROM" ${id}.filtered.indels.vcf | cut -d':' -f 1`
  head -n $header ${id}.filtered.indels.vcf > ${id}.indel_head
  cat ${id}.filtered.indels.vcf | grep PASS | cat ${id}.indel_head - > ${id}.PASS.indels.vcf

gatk VariantFiltration -R  ${reference} -O ${id}.failed.indels.vcf -V ${id}.raw.indels.vcf \
-filter "MLEAF < ${MLEAF_INDEL}" --filter-name "FAIL" \
-filter "MQ < ${MQ_INDEL}" --filter-name "FAIL1" \
-filter "QD < ${QD_INDEL}" --filter-name "FAIL2" \
-filter "FS > ${FS_INDEL}" --filter-name "FAIL3" \
-filter "QUAL < ${QUAL_INDEL}" --filter-name "FAIL5"

  header=`grep -a -n "#CHROM" ${id}.failed.indels.vcf | cut -d':' -f 1`
  head -n $header ${id}.failed.indels.vcf > ${id}.indel_head
  cat ${id}.failed.indels.vcf | grep FAIL | cat ${id}.indel_head - > ${id}.FAIL.indels.vcf

  ${SNPEFF_CMD} eff -nodownload -no-downstream -no-intergenic -ud 100 -v -dataDir ${baseDir}/resources/snpeff ${snpeff} ${id}.PASS.snps.vcf > ${id}.PASS.snps.annotated.vcf

  ${SNPEFF_CMD} eff -nodownload -no-downstream -no-intergenic -ud 100 -v -dataDir ${baseDir}/resources/snpeff ${snpeff} ${id}.PASS.indels.vcf > ${id}.PASS.indels.annotated.vcf

echo "Calculating duplication and deletion events"

echo -e "Chromosome\tStart\tEnd\tInterval" > ${id}.tmp.header
zcat output.per-base.bed.gz | awk '$4 ~ /^0/ { print $1,$2,$3,$3-$2 }' > ${id}.del.summary.tmp
cat ${id}.tmp.header ${id}.del.summary.tmp > ${id}.deletion_summary.txt

covdep=$(head -n 1 ${id}.depth.txt)
DUP_CUTOFF=$(echo $covdep*3 | bc)
echo "dup cutoff is $DUP_CUTOFF"

zcat output.per-base.bed.gz | awk -v DUP_CUTOFF="$DUP_CUTOFF" '$4 >= DUP_CUTOFF { print $1,$2,$3,$3-$2 }' > ${id}.dup.summary.tmp

i=$(head -n1 ${id}.dup.summary.tmp | awk '{ print $2 }')
k=$(tail -n1 ${id}.dup.summary.tmp | awk '{ print $3 }')
chr=$(head -n1 ${id}.dup.summary.tmp | awk '{ print $1 }')

awk -v i="$i" -v k="$k" -v chr="$chr" 'BEGIN {printf "chromosome " chr " start " i " "; j=i} {if (i==$2 || i==$2-1 || i==$2-2 ) {
i=$3;
}
else {
  print "end "i " interval " i-j;
  j=$2;
  i=$3;
  printf "chromosome " $1 " start "j " ";
}} END {print "end "k " interval "k-j}' < ${id}.dup.summary.tmp > ${id}.dup.summary.tmp1

sed -i 's/chromosome\|start \|end \|interval //g' ${id}.dup.summary.tmp1
echo -e "Chromosome\tStart\tEnd\tInterval" > ${id}.dup.summary.tmp.header
cat ${id}.dup.summary.tmp.header ${id}.dup.summary.tmp1 > ${id}.duplication_summary.txt

awk '{
  if (match($0,"ANN=")){print substr($0,RSTART)}
  }' ${id}.PASS.indels.annotated.vcf > ${id}.indel.effects.tmp

awk -F "|" '{ print $4,$10,$11,$15 }' ${id}.indel.effects.tmp | sed 's/c\.//' | sed 's/p\.//' | sed 's/n\.//'> ${id}.annotated.indel.effects

awk '{
  if (match($0,"ANN=")){print substr($0,RSTART)}
  }' ${id}.PASS.snps.annotated.vcf > ${id}.snp.effects.tmp
awk -F "|" '{ print $4,$10,$11,$15 }' ${id}.snp.effects.tmp | sed 's/c\.//' | sed 's/p\.//' | sed 's/n\.//' > ${id}.annotated.snp.effects

echo 'Identifying high consequence mutations'

grep 'HIGH' ${id}.snp.effects.tmp  | awk -F"|" '{ print $4,$11 }' >> ${id}.Function_lost_list.txt
grep 'HIGH' ${id}.indel.effects.tmp | awk -F"|" '{ print $4,$11 }' >> ${id}.Function_lost_list.txt

sed -i 's/p\.//' ${id}.Function_lost_list.txt

delly call -q 5 -o ${id}.delly.bcf -g ${reference} ${id}.dedup.bam
bcftools view ${id}.delly.bcf > ${id}.delly.vcf
grep "#" ${id}.delly.vcf > ${id}.delly.header
grep "<INV>" ${id}.delly.vcf > ${id}.delly.inv.vcf
grep -v "LowQual" ${id}.delly.inv.vcf > ${id}.delly.inv.vcf.tmp
cat ${id}.delly.header ${id}.delly.inv.vcf.tmp > ${id}.delly.inv.vcf

  ${SNPEFF_CMD} eff -no-downstream -no-intergenic -ud 100 -v -dataDir ${baseDir}/resources/snpeff ${snpeff} ${id}.delly.inv.vcf > ${id}.delly.inv.annotated.vcf

if [ -s ${id}.delly.inv.vcf.tmp ]; then
  bcftools query -f '%CHROM %POS[\t%GT\t%GL]\n' ${id}.delly.inv.vcf > ${id}.likelihoods.delly
  while read line; do
    echo "$line" > ${id}.line.desc;
    awk '{print $4}' ${id}.line.desc > ${id}.geno.likelihoods;
    genotype_RR=$(awk -F"," '{print $1}' ${id}.geno.likelihoods);
    genotype_RA=$(awk -F"," '{print $2}' ${id}.geno.likelihoods);
    genotype_AA=$(awk -F"," '{print $3}' ${id}.geno.likelihoods);
    genotype=$(awk '{print $3}' ${id}.line.desc);
    if [ "$genotype" == "0/1"  ]; then
      if [ "$genotype_RR" == 0 ]; then
      echo "Genotype ignored";
      fi;
    if [ "$genotype_AA" == 0 ]; then
      echo "Genotype included";
      chromosome=$(awk '{print $1}' ${id}.line.desc);
      location=$(awk '{print $2}' ${id}.line.desc);
      echo -e "$chromosome\t$location" >> ${id}.filtered.inversions;
    fi;
    if [ "$genotype_RA" == 0 ]; then
      alt_ref_check=0;
      alt_ref_check=$(awk -v a="$genotype_RR" -v b="$genotype_AA" 'BEGIN {if (a < b) {print "1" }}');
      if [ "$alt_ref_check" == 1 ]; then
        echo "calculating log likelihood";
        log_genotype_AA=$(awk -v a="$genotype_AA" 'BEGIN {print (10^a)}');
        log_genotype_RA=$(awk -v a="$genotype_RA" 'BEGIN {print (10^a)}');
        log_genotype_RR=$(awk -v a="$genotype_RR" 'BEGIN {print (10^a)}');
        sum_AA_RR=$(awk -v a="$log_genotype_AA" -v b="$log_genotype_RR" 'BEGIN {print (a+b)}' );
        likelihood_ratio=$(awk -v a="$log_genotype_RA" -v b="$sum_AA_RR" 'BEGIN {print (a/b)}');
        echo -e "$log_genotype_AA\t$log_genotype_RA\t$log_genotype_RR" >> ${id}.likelihood.ratios.2
        echo -e "$likelihood_ratio\n" >> ${id}.likelihood.ratios.2
        likelihood_ratio_test=$(awk -v a="$likelihood_ratio" 'BEGIN {if (a < 100000) {print "1" }}')
        if [ "$likelihood_ratio_test" == 1 ]; then
          echo "changing genotype to 1/1";
          chromosome=$(awk '{print $1}' ${id}.line.desc);
          location=$(awk '{print $2}' ${id}.line.desc);
          echo -e "$chromosome\t$location" >> ${id}.filtered.inversions;
        else
          echo "Ignoring variant due to poor quality"
        fi;
      else
        echo "Ignoring variant due to likely reference allele"
      fi;
    fi;
  fi;
  if [ "$genotype" == "1/1"  ]; then
    echo "Genotype included";
    chromosome=$(awk '{print $1}' ${id}.line.desc);
    location=$(awk '{print $2}' ${id}.line.desc);
    echo -e "$chromosome\t$location" >> ${id}.filtered.inversions;
  fi;
  done < ${id}.likelihoods.delly
fi;

if [ -s ${id}.filtered.inversions ]; then
  while read line; do grep -w "$line" ${id}.delly.inv.annotated.vcf >> ${id}.delly.inv.annotated.vcf.tmp ; done < ${id}.filtered.inversions
  cat ${id}.delly.header ${id}.delly.inv.annotated.vcf.tmp > ${id}.delly.inv.annotated.vcf
  awk -F"|" '/HIGH/ {f=NR} f&&NR-1==f' RS="|" ${id}.delly.inv.annotated.vcf > ${id}.delly.tmp
  sed -i '/^\s*$/d' ${id}.delly.tmp
  cat ${id}.delly.tmp ${id}.Function_lost_list.txt > ${id}.Function_lost_list.txt.tmp
  mv ${id}.Function_lost_list.txt.tmp ${id}.Function_lost_list.txt
fi;