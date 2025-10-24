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

${SNPEFF_CMD} genes2bed ${snpeff} -dataDir ${baseDir}/resources/snpeff -f ${id}.gene.list > ${id}.intervals.list
sed -i 's/\t/:/' ${id}.intervals.list
sed -i 's/\t/-/' ${id}.intervals.list
awk '{print $1}' ${id}.intervals.list | tail -n+2 > ${id}.intervals.list.tmp
mv ${id}.intervals.list.tmp ${id}.intervals.list

#bedtools maskfasta -fi ref.fasta -bed intervals.list -fo ref.intervals

#mv ref.intervals ref.fasta
