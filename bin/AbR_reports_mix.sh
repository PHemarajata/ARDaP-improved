#!/bin/bash

seq=$1
RESISTANCE_DB=$2

#This variable restricts resfinder matches to "3" which means 100% nucleotide identity across the entire length of the gene. Set to "2" to reduce stringency 
stringency=3

echo -e "Importing isolate data"
#echo -e "ID,Barcode,LName,FName,DOB,Location,sampType,sampID,sampDate,sampSource,sampSeq,reportLab,reportDate,comments,organism,requestor,requestorContact,lineageNum,lineageName" > patientMetaData.csv
date=$(date +"%F")
echo -e "Looking for specific strain in the metadata file"
grep -w "$seq" patientMetaData.csv
status=$?
if [ $status == 0 ]; then
    echo "Found strain information"
	head -n1 patientMetaData.csv >> patientMetaData_"$seq".csv
	grep -w "$seq" patientMetaData.csv >> patientMetaData_"$seq".csv
else
    echo "Couldn't find strain specific data in file, reverting to default"
	head -n1 patientMetaData.csv >> patientMetaData_"$seq".csv
    echo -e "$seq,BARCODE,Smith,James,1/01/1990,Shanghai,Blood,$seq,$date,Blood,Cultured isolate,unspecified,$date,No words needed,Burkholderia pseudomallei,Dr. Requestor Name,req_contact@genome.com,XX,NA" >> patientMetaData_"$seq".csv	  
fi  

mv patientMetaData_"$seq".csv patientMetaData.csv

Report_structure () {

cat << _EOF_ >  ${seq}.Drug.table
.separator ","
SELECT Antibiotics."Drug_class",
Antibiotics.Antibiotic,
Antibiotics."Abbreviation"
FROM Antibiotics
ORDER BY Antibiotics."Drug_class"
_EOF_

sqlite3 "$RESISTANCE_DB" < ${seq}.Drug.table >> ${seq}.drug.table.txt

}

Report_structure

cp ${seq}.drug.table.txt ${seq}.drug.table.txt.backup

#sed manipulation of the drug table here is simplistic and relies on the order of the drugs in the drug.table.txt.back file
#TO DO -  replace with awk pattern matching is case users want to add custom drug classes

#format Resfinder output for report
while read line; do
    echo "$line" > ${seq}_resfinder_query_string.txt
    Antibiotic=$(awk -F "," '{print $2}' ${seq}_resfinder_query_string.txt)
    abbrev=$(awk -F "," '{print $3}' ${seq}_resfinder_query_string.txt)
  echo "Antibiotic = $Antibiotic"
  echo "abbrev=$abbrev"
    awk -F "\t" -v Ab="${Antibiotic}" 'BEGIN{IGNORECASE=1} $1==Ab' "${seq}"_resfinder.txt > ${seq}_resfinder_tmp.txt
    gene=$(awk -F "\t" '{ print $5 }' ${seq}_resfinder_tmp.txt)
    }
      for(j=1; j<=p; j++) {
    }
    for (i=1; i<=NF; i++) {
    }
    END {
Report_structure () {

cat << _EOF_ >  ${seq}.Drug.table
.separator ","
SELECT Antibiotics."Drug_class",
Antibiotics.Antibiotic,
Antibiotics."Abbreviation"
FROM Antibiotics
ORDER BY Antibiotics."Drug_class"
_EOF_

# write per-sample drug table (overwrite to avoid cross-run contamination)
sqlite3 "$RESISTANCE_DB" < ${seq}.Drug.table > ${seq}.drug.table.txt
# remove temporary SQL file
rm -f ${seq}.Drug.table

}

Report_structure

cp ${seq}.drug.table.txt ${seq}.drug.table.txt.backup

#format Resfinder output for report
while read line; do
    echo "$line" > ${seq}_resfinder_query_string.txt
    Antibiotic=$(awk -F "," '{print $2}' ${seq}_resfinder_query_string.txt)
    abbrev=$(awk -F "," '{print $3}' ${seq}_resfinder_query_string.txt)
  echo "Antibiotic = $Antibiotic"
  echo "abbrev=$abbrev"
    awk -F "\t" -v Ab="${Antibiotic}" 'BEGIN{IGNORECASE=1} $1==Ab' "${seq}"_resfinder.txt > ${seq}_resfinder_tmp.txt
    gene=$(awk -F "\t" '{ print $5 }' ${seq}_resfinder_tmp.txt)
  #apply resfinder stringency filter
  #awk -v stringency="${stringency}" -F "\t" '$4==stringency' resfinder_tmp.txt > resfinder_tmp2.txt
 
  #look for resistance
  awk -F "\t" '$3 ~ "Resistant" { exit 1 } ' resfinder_tmp.txt &> /dev/null
  status=$?
  #echo $status  
  if [ "$status" == 1 ]; then  #Resistant strain
    echo "Found resistance to ${Antibiotic} for ${seq}"
    echo -e "Resfinder|$gene||${abbrev}r|100" >> "${seq}"_resfinder_report.txt #If the resistance mechanism is found then report 100% for resfinder AMR
  else
    echo "Strain appears sensitive to $Antibiotic"
  fi
        }


# ensure AbR output files exist (use /dev/null fallback)
snp_file="${seq}.AbR_output_snp_indel_mix.txt"
if [ ! -s "$snp_file" ]; then
    snp_file="${seq}.AbR_output_snp_indel.txt"
fi
del_file="${seq}.AbR_output_del_dup_mix.txt"
if [ ! -s "$del_file" ]; then
    del_file="${seq}.AbR_output_del_dup.txt"
fi
if [ ! -s "$snp_file" ]; then
    snp_file="/dev/null"
fi
if [ ! -s "$del_file" ]; then
    del_file="/dev/null"
fi

cat "$snp_file" "$del_file" "${seq}"_resfinder_report.txt | tee ${seq}.AbR_output.txt ${seq}.AbR_output.final.txt


#Deduplicate any repetition in the resistance list
awk -F"|" '!seen[$1,$2,$3,$4,$5]++' ${seq}.AbR_output.final.txt > ${seq}.AbR_output.temp
mv ${seq}.AbR_output.temp ${seq}.AbR_output.final.txt

#Deduplicate any repetition in the resistance list
awk -F"|" '!seen[$1,$2,$3,$4,$5]++' ${seq}.AbR_output.txt > ${seq}.AbR_output.temp
mv ${seq}.AbR_output.temp ${seq}.AbR_output.txt

# Robust CSV generator: iterate per-sample drug table backup and inspect deduplicated AbR output
out_csv="${seq}_patientDrugSusceptibilityData.csv"
echo "ID,Class,Drug,Status,Level,Details" > "$out_csv"

while IFS=, read -r drug_class drug_name abbrev; do
    # skip header lines if present
    if [[ "$abbrev" == "Abbreviation" || -z "$abbrev" ]]; then
        continue
    fi
    # find matching lines in AbR output
    matches=$(awk -F"|" -v a="$abbrev" 'BEGIN{IGNORECASE=1} $4 ~ a"r" || $4 ~ a"i" || $4 ~ a"s" { print $0 }' ${seq}.AbR_output.final.txt)
    status_val="Sensitive"
    level_val=0
    details="No resistance detected"
    if [[ -n "$matches" ]]; then
        if echo "$matches" | grep -iq "${abbrev}r"; then
            status_val="Resistant"
        elif echo "$matches" | grep -iq "${abbrev}i"; then
            status_val="Intermediate"
        else
            status_val="Sensitive"
        fi
        level_sum=0
        details_list=
        while IFS= read -r m; do
            lvl=$(echo "$m" | awk -F"," '{print $NF}' | tr -d '\r')
            if [[ "$lvl" =~ ^[0-9]+$ ]]; then
                level_sum=$((level_sum + lvl))
            fi
            gene=$(echo "$m" | awk -F"|" '{print $2}' )
            mut=$(echo "$m" | awk -F"|" '{print $3}' )
            details_list+="${gene} ${mut};"
        done <<< "$matches"
        level_val=$level_sum
        details=$(echo "$details_list" | sed 's/;$/\n/' | tr -d '\n')
    fi
    printf '%s,%s,%s,%s,%s,%s\n' "$seq" "$drug_class" "$drug_name" "$status_val" "$level_val" "$details" >> "$out_csv"
done < ${seq}.drug.table.txt.backup

# move per-sample result into canonical name used by Report_html.sh
mv "$out_csv" patientDrugSusceptibilityData.csv

# also ensure per-sample AbR_output.final copy exists
cp ${seq}.AbR_output.final.txt "${seq}.AbR_output.final.txt"

exit 0
    }
    NF>p { p = NF }
    END {
      for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str"\t"a[i,j];
        }
        print str
       }
    }' | sed 's/\t/\n/g' | grep "$f" -A1 > "$f"s.level_calc
	
	cat "$f"{r,s,i}.level_calc > "$f".level_calc.tmp
	
	sum_resistance=$(awk '{ for(i=1; i<=NF; i++) { if($i ~ /^-?[0-9]+$/) { sum += $i } } } END { print sum }' "$f".level_calc.tmp)
	#echo "$sum_resistance"
	if [ ! -n "$sum_resistance" ]; then
	  sum_resistance=0
    fi
	
	grep -w "$f"r "$f"r.output &> /dev/null
	status=$?
	if [[ "$status" -eq 0 ]]; then
		echo "found mechanism for $f resistance"
		length=$(wc -l "$f"r.output | awk '{print $1}' )
		if [[ "$length" -gt 1 ]]; then
			echo "found multiple determinants for $f resistance"
			#cat "$f"s.output >> drug.table.tertiary.txt
			sed -i "${i}s/.*/&,Resistant,${sum_resistance},Multiple determinants/" drug.table.txt
			i=$((i+1))
		else
			echo "found single mechanism for $f resistance" 
			mech=$(awk -F "|" '{ print $1,$2,$3 }' "$f"r.output) #Prints gene name (column 2 from SQL query) and mutation (#col 3
			sed -i "${i}s/.*/&,Resistant,${sum_resistance},${mech}/" drug.table.txt
			#cat "$f"s.output >> drug.table.tertiary.txt
			i=$((i+1))
		fi
	else
		echo "no mechanism identified for $f resistance"
		sed -i "${i}s/.*/&,Sensitive,${sum_resistance},No resistance detected/" drug.table.txt
		i=$((i+1))
	fi
done < <(grep -E "tertiary|Tertiary" drug.table.txt.backup | awk -F "," '{ print $3 }')

#Looking for speciation and QC
echo "Looking for speciation markers"
awk -F"|" -v f="speciation" 'BEGIN{IGNORECASE=1} $4~ f' AbR_output.txt > speciation.output
grep -iw "speciation" speciation.output &> /dev/null
status=$?
if [[ "$status" -eq 0 ]]; then
	echo "Species specific marker is missing. Likely quality control issue" #TO DO This should go into final report
	echo "Interpret AMR profile with caution" #TO DO This should go into final report
fi	

# create patientDrugSusceptibilityData.csv
# ID refers to individual strains
sed -i "s/^/$seq,/" drug.table.txt
awk -v FS="," -v OFS="," '{print $1,$2,$3,$5,$6,$7 }' drug.table.txt > drug.table.txt.tmp

#debug line to check format
#cp drug.table.txt drug.table.txt.format_check

mv drug.table.txt.tmp drug.table.txt
    sed -i '1 i\ID,Class,Drug,Status,Level,Details' ${seq}.drug.table.txt 
    cp ${seq}.drug.table.txt ${seq}_patientDrugSusceptibilityData.csv

if [ -s ${seq}.AbR_output.final.txt ]; then
    cp ${seq}.AbR_output.final.txt "${seq}.AbR_output.final.txt"
else	
    echo "No antibiotic resistance identified in $seq" >> "${seq}.AbR_output.final.txt"
fi

# move per-sample patientDrugSusceptibilityData into the expected final filename
mv ${seq}_patientDrugSusceptibilityData.csv patientDrugSusceptibilityData.csv

exit 0
