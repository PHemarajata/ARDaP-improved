#!/bin/bash

cpus=$1
echo "task.cpus=$cpus"

forward=$2
echo "forward=$forward"

reverse=$3
echo "reverse=$reverse"

id=$4
echo "id=$id"

#baseDir for snpEff
baseDir=$5
echo "baseDir=$baseDir"

#annotation genome for snpeff
snpeff=$6
echo "snpEff=$snpeff"

#create interval list
cat << _EOF_ > ${id}.interval.coverage.query.txt
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
cat << _EOF_ > ${id}.interval.query.txt
SELECT 
	Variants_SNP_indel.Gene_name
FROM 
	Variants_SNP_indel
_EOF_

sqlite3 ${baseDir}/Databases/${snpeff}/${snpeff}.db < ${id}.interval.query.txt > ${id}.gene.list
sqlite3 ${baseDir}/Databases/${snpeff}/${snpeff}.db < ${id}.interval.coverage.query.txt > ${id}.gene.list2
  
cat ${id}.gene.list ${id}.gene.list2 ${id}.interval.coverage.query.txt > ${id}.gene.list.tmp
mv ${id}.gene.list.tmp ${id}.gene.list 

#create interval file

uniq ${id}.gene.list > ${id}.gene.list.tmp
mv ${id}.gene.list.tmp ${id}.gene.list
${SNPEFF_CMD:-snpEff} genes2bed ${snpeff} -dataDir ${baseDir}/resources/snpeff -f ${id}.gene.list > ${id}.intervals.list
sed -i 's/\t/:/' ${id}.intervals.list
sed -i 's/\t/-/' ${id}.intervals.list
awk '{print $1}' ${id}.intervals.list | tail -n+2 > ${id}.intervals.list.tmp
mv ${id}.intervals.list.tmp ${id}.intervals.list

#bedtools maskfasta -fi ref.fasta -bed intervals.list -fo ref.intervals

#mv ref.intervals ref.fasta
