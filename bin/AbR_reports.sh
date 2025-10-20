#!/bin/bash

seq=$1
RESISTANCE_DB=$2


#This variable restricts resfinder matches to "3" which means 100% nucleotide identity across the entire length of the gene. Set to "2" to reduce stringency 
stringency=3

echo -e "Importing isolate data"
date=$(date +"%F")
# work on a per-sample copy of the metadata to avoid clobbering the repo-level file
cp patientMetaData.csv ${seq}_patientMetaData.csv
echo -e "Looking for specific strain in the metadata file"
grep -w "$seq" ${seq}_patientMetaData.csv
status=$?
if [ $status == 0 ]; then
    echo "Found strain information"
    head -n1 ${seq}_patientMetaData.csv >> ${seq}_patientMetaData_${seq}.csv
    grep -w "$seq" ${seq}_patientMetaData.csv >> ${seq}_patientMetaData_${seq}.csv
else
    echo "Couldn't find strain specific data in file, reverting to default"
    head -n1 ${seq}_patientMetaData.csv >> ${seq}_patientMetaData_${seq}.csv
    echo -e "$seq,BARCODE,Smith,James,1/01/1990,Shanghai,Blood,$seq,$date,Blood,Cultured isolate,unspecified,$date,No words needed,Burkholderia pseudomallei,Dr. Requestor Name,req_contact@genome.com,XX,NA" >> ${seq}_patientMetaData_${seq}.csv
fi

# use the per-sample metadata file going forward
mv ${seq}_patientMetaData_${seq}.csv patientMetaData_${seq}.csv

Report_structure () {

cat << _EOF_ >  Drug.table
.separator ","
SELECT Antibiotics."Drug_class",
Antibiotics.Antibiotic,
Antibiotics."Abbreviation"
FROM Antibiotics
ORDER BY Antibiotics."Drug_class"
_EOF_


sqlite3 "$RESISTANCE_DB" < Drug.table > ${seq}.drug.table.txt
# remove temporary SQL file
rm -f Drug.table

}

Report_structure

cp ${seq}.drug.table.txt ${seq}.drug.table.txt.backup

# remove any header lines that might have been produced by sqlite
sed -i '/Drug_class\|Antibiotics\|Abbreviation/d' ${seq}.drug.table.txt
sed -i '/Drug_class\|Antibiotic\|Abbreviation/d' ${seq}.drug.table.txt.backup 2>/dev/null || true

#format Resfinder output for report
while read line; do
    echo "$line" > ${seq}_resfinder_query_string.txt
    Antibiotic=$(awk -F "," '{print $2}' ${seq}_resfinder_query_string.txt)
    abbrev=$(awk -F "," '{print $3}' ${seq}_resfinder_query_string.txt)
    # create a filesystem-safe abbreviation (replace non-alphanum with _)
    safe_abbrev=$(echo "${abbrev}" | sed 's/[^A-Za-z0-9]/_/g')
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
    echo -e "Resfinder|$gene||${safe_abbrev}r|100" >> "${seq}"_resfinder_report.txt #If the resistance mechanism is found then report 100% for resfinder AMR
  else
    echo "Strain appears sensitive to $Antibiotic"
  fi
done < ${seq}.drug.table.txt.backup


# Ensure placeholder files exist so concatenation doesn't fail when upstream produced no file
[ -f "${seq}.AbR_output_snp_indel.txt" ] || touch "${seq}.AbR_output_snp_indel.txt"
[ -f "${seq}.AbR_output_del_dup.txt" ] || touch "${seq}.AbR_output_del_dup.txt"
[ -f "${seq}_resfinder_report.txt" ] || touch "${seq}_resfinder_report.txt"

cat "${seq}".AbR_output_snp_indel.txt "${seq}".AbR_output_del_dup.txt "${seq}"_resfinder_report.txt | tee ${seq}.AbR_output.txt ${seq}.AbR_output.final.txt


#Deduplicate any repetition in the resistance list
awk -F"|" '!seen[$1,$2,$3,$4,$5]++' ${seq}.AbR_output.final.txt > ${seq}.AbR_output.temp
mv ${seq}.AbR_output.temp ${seq}.AbR_output.final.txt

#Deduplicate any repetition in the resistance list
awk -F"|" '!seen[$1,$2,$3,$4,$5]++' ${seq}.AbR_output.txt > ${seq}.AbR_output.temp
mv ${seq}.AbR_output.temp ${seq}.AbR_output.txt

# --- Generate patientDrugSusceptibilityData.csv now (moved up) ---
# create patientDrugSusceptibilityData.csv (robust generator)
out=${seq}.drug.table.txt
tmp=${seq}.drug.table.txt.tmp2
echo "ID,Class,Drug,Status,Level,Details" > "$tmp"

# If backup doesn't exist, fallback to the generated drug table
if [ ! -f ${seq}.drug.table.txt.backup ]; then
    cp ${seq}.drug.table.txt ${seq}.drug.table.txt.backup
fi

while IFS=',' read -r dclass ddrug dabbrev; do
    # skip empty lines or header-like lines
    if [[ -z "$dabbrev" || "$dclass" == "Drug_class" || "$dabbrev" == "Abbreviation" ]]; then
        continue
    fi
    abbrev=$(echo "$dabbrev" | sed 's/[^A-Za-z0-9]//g')
    # default values
    status="Sensitive"
    level=0
    details="No resistance detected"

    if [ -s ${seq}.AbR_output.final.txt ]; then
        # find any matching records (case-insensitive) where the 4th field contains the abbreviation + r/i/s
        # prefer Resistant > Intermediate > Sensitive
        r_match=$(awk -F"|" -v a="${abbrev}" 'BEGIN{IGNORECASE=1} tolower($4) ~ tolower(a)"r" {print; exit}' ${seq}.AbR_output.final.txt)
        i_match=$(awk -F"|" -v a="${abbrev}" 'BEGIN{IGNORECASE=1} tolower($4) ~ tolower(a)"i" {print; exit}' ${seq}.AbR_output.final.txt)
        s_match=$(awk -F"|" -v a="${abbrev}" 'BEGIN{IGNORECASE=1} tolower($4) ~ tolower(a)"s" {print; exit}' ${seq}.AbR_output.final.txt)

        if [ -n "$r_match" ]; then
            status="Resistant"
            # collect all mechanics for this abbrev
            details=$(awk -F"|" -v a="${abbrev}" 'BEGIN{IGNORECASE=1} tolower($4) ~ tolower(a)"r" {print $1" "$2" "$3}' ${seq}.AbR_output.final.txt | tr '\n' ';' | sed 's/;$//')
            level=$(awk -F"|" -v a="${abbrev}" 'BEGIN{IGNORECASE=1} { if(tolower($4) ~ tolower(a)"r" ){ g=$5; gsub(/[^0-9-]/,"",g); sum+=g } } END{print (sum+0)}' ${seq}.AbR_output.final.txt)
        elif [ -n "$i_match" ]; then
            status="Intermediate"
            details=$(awk -F"|" -v a="${abbrev}" 'BEGIN{IGNORECASE=1} tolower($4) ~ tolower(a)"i" {print $1" "$2" "$3}' ${seq}.AbR_output.final.txt | tr '\n' ';' | sed 's/;$//')
            level=$(awk -F"|" -v a="${abbrev}" 'BEGIN{IGNORECASE=1} { if(tolower($4) ~ tolower(a)"i" ){ g=$5; gsub(/[^0-9-]/,"",g); sum+=g } } END{print (sum+0)}' ${seq}.AbR_output.final.txt)
        elif [ -n "$s_match" ]; then
            status="Sensitive"
            details=$(awk -F"|" -v a="${abbrev}" 'BEGIN{IGNORECASE=1} tolower($4) ~ tolower(a)"s" {print $1" "$2" "$3}' ${seq}.AbR_output.final.txt | tr '\n' ';' | sed 's/;$//')
            level=$(awk -F"|" -v a="${abbrev}" 'BEGIN{IGNORECASE=1} { if(tolower($4) ~ tolower(a)"s" ){ g=$5; gsub(/[^0-9-]/,"",g); sum+=g } } END{print (sum+0)}' ${seq}.AbR_output.final.txt)
        fi
    fi

    # Write single line for this drug
    echo "${seq},${dclass},${ddrug},${status},${level},${details}" >> "$tmp"
done < <(awk -F"," 'NF>=3{print $1","$2","$3}' ${seq}.drug.table.txt.backup)

# move into place
mv "$tmp" ${seq}.drug.table.txt
cp ${seq}.drug.table.txt ${seq}_patientDrugSusceptibilityData.csv

if [ -s ${seq}.AbR_output.final.txt ]; then
    cp ${seq}.AbR_output.final.txt "${seq}.AbR_output.final.txt"
else
    echo "No antibiotic resistance identified in $seq" >> "${seq}.AbR_output.final.txt"
fi

# move per-sample patientDrugSusceptibilityData into the expected final filename
mv ${seq}_patientDrugSusceptibilityData.csv patientDrugSusceptibilityData.csv

# cleanup per-sample metadata copy
rm -f ${seq}_patientMetaData.csv

exit 0

#TO DO -  replace with awk pattern matching in case users want to add custom drug classes

i=1
while read f; do 
    # sanitize f into safe token (replace non-alphanum with underscore)
    safe_f=$(echo "$f" | sed 's/[^A-Za-z0-9]/_/g')
    awk -F"|" -v f="$f" '$4~ f"r"' ${seq}.AbR_output.txt > "${seq}_${safe_f}r.output"
    awk -F"|" -v f="$f" '$4~ f"i"' ${seq}.AbR_output.txt > "${seq}_${safe_f}i.output"
    awk -F"|" -v f="$f" '$4~ f"s"' ${seq}.AbR_output.txt > "${seq}_${safe_f}s.output"
	#awk -F"|" -v f="$f" '$4~ f"r"' "${seq}"_resfinder_report.txt >> "$f"r.output #already included above
	#awk -F"|" -v f="$f" '$4~ f"i"' "${seq}"_resfinder_report.txt >> "$f"i.output
	
	
	##calc level of resistance - messy code here. TO DO - neaten and reduce redundancy 
	
		awk -F"|" 'BEGIN { OFS = "\n" } {print $4,$5}' "$f"r.output | awk -F"," '{
    for (i=1; i<=NF; i++) {
        a[NR,i] = $i
        }
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
    }' | sed 's/\t/\n/g' | grep "$f" -A1 > "$f"r.level_calc
	
	awk -F"|" 'BEGIN { OFS = "\n" } {print $4,$5}' "$f"i.output | awk -F"," '{
    for (i=1; i<=NF; i++) {
        a[NR,i] = $i
        }
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
    }' | sed 's/\t/\n/g' | grep "$f" -A1 > "$f"i.level_calc
	
	awk -F"|" 'BEGIN { OFS = "\n" } {print $4,$5}' "$f"s.output | awk -F"," '{
    for (i=1; i<=NF; i++) {
        a[NR,i] = $i
        }
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
	
    cat "${seq}_${safe_f}"{r,s,i}.level_calc > "${seq}_${safe_f}".level_calc.tmp
    
    sum_resistance=$(awk '{ for(i=1; i<=NF; i++) { if($i ~ /^-?[0-9]+$/) { sum += $i } } } END { print sum }' "${seq}_${safe_f}".level_calc.tmp)
	
	#echo "$sum_resistance"
	if [ ! -n "$sum_resistance" ]; then
	  sum_resistance=0
    fi
	
    grep -w "$f"r "${seq}_${safe_f}"r.output &> /dev/null #looks for full resistance
	status=$?
	if [[ "$status" -eq 0 ]]; then 
		echo "found mechanism for $f resistance"
		length=$(wc -l "$f"r.output | awk '{print $1}' )
		if [[ "$length" -gt 1 ]]; then
			echo "found multiple determinants for $f resistance"
            sed -i "${i}s/.*/&,Resistant,${sum_resistance},Multiple determinants/" ${seq}.drug.table.txt
			i=$((i+1))
		else
			echo "found single mechanism for $f resistance" 
            mech=$(awk -F "|" '{ print $1,$2,$3 }' "${seq}_${safe_f}"r.output) #Prints gene name (column 2 from SQL query) and mutation (col 3
            sed -i "${i}s/.*/&,Resistant,${sum_resistance},${mech}/" ${seq}.drug.table.txt
			i=$((i+1))
		fi
	else
		echo "no mechanism identified for $f resistance, looking for intermediate resistance"
		grep -w "${f}"i "$f"i.output &> /dev/null
		status=$?
		if [[ "$status" -eq 0 ]]; then
			echo "found intermediate resistance mechanism for $f"
			length=$(wc -l "$f"i.output | awk '{print $1}' )
			if [[ "$length" -gt 1 ]]; then
				echo "found multiple determinants for intermediate $f resistance"
                sed -i "${i}s/.*/&,Intermediate,${sum_resistance},Multiple determinants/" ${seq}.drug.table.txt
				i=$((i+1))
			else
				echo "found single mechanism for intermediate $f resistance" 
				mech=$(awk -F "|" '{ print $1,$2,$3 }' "$f"i.output) #Prints gene name (column 2 from SQL query) and mutation (col 3
                sed -i "${i}s/.*/&,Intermediate,${sum_resistance},${mech}/" ${seq}.drug.table.txt
				i=$((i+1))
			fi
		else
			echo "no intermediate resistance found"
            sed -i "${i}s/.*/&,Sensitive,${sum_resistance},No resistance detected/" ${seq}.drug.table.txt
			i=$((i+1))
		fi
	fi
done < <(grep -E "First-line|first-line" ${seq}.drug.table.txt.backup | awk -F "," '{ print $3 }') #Three letter abbreviation of antibiotic

# --- CHANGED: operate on per-sample AbR output and drug table for sensitivity section ---
while read f; do
    safe_f=$(echo "$f" | sed 's/[^A-Za-z0-9]/_/g')
    awk -F"|" -v f="$f" '$4~ f"s" ' ${seq}.AbR_output.txt > "${seq}_${safe_f}s.output" 
	#awk -F"|" -v f="$f" '$4~ f"s" ' "${seq}"_resfinder_report.txt >> "$f"s.output
	
	##calc level of resistance - messy code here. TO DO - neaten and reduce redundancy 
	
	awk -F"|" 'BEGIN { OFS = "\n" } {print $4,$5}' "$f"r.output | awk -F"," '{
    for (i=1; i<=NF; i++) {
        a[NR,i] = $i
        }
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
    }' | sed 's/\t/\n/g' | grep "$f" -A1 > "$f"r.level_calc
	
	awk -F"|" 'BEGIN { OFS = "\n" } {print $4,$5}' "$f"i.output | awk -F"," '{
    for (i=1; i<=NF; i++) {
        a[NR,i] = $i
        }
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
    }' | sed 's/\t/\n/g' | grep "$f" -A1 > "$f"i.level_calc
	
	awk -F"|" 'BEGIN { OFS = "\n" } {print $4,$5}' "$f"s.output | awk -F"," '{
    for (i=1; i<=NF; i++) {
        a[NR,i] = $i
        }
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

    cat "${seq}_${safe_f}"{r,s,i}.level_calc > "${seq}_${safe_f}".level_calc.tmp

    sum_resistance=$(awk '{ for(i=1; i<=NF; i++) { if($i ~ /^-?[0-9]+$/) { sum += $i } } } END { print sum }' "${seq}_${safe_f}".level_calc.tmp)
	#echo "$sum_resistance"
	if [ ! -n "$sum_resistance" ]; then
	  sum_resistance=0
    fi
	
    grep -w "$f"s "${seq}_${safe_f}"s.output &> /dev/null
	status=$?
	if [[ "$status" -eq 0 ]]; then
		echo "found mechanism for $f sensitivity"
		length=$(wc -l "$f"s.output | awk '{print $1}' )
		if [[ "$length" -gt 1 ]]; then
			echo "found multiple determinants for $f sensitivity"
			#cat "$f"s.output >> ${seq}.drug.table.tertiary.txt
            sed -i "${i}s/.*/&,Sensitive,${sum_resistance},Multiple determinants/" ${seq}.drug.table.txt
			i=$((i+1))
		else
			echo "found single mechanism for $f sensitivity" 
            mech=$(awk -F "|" '{ print $1,$2,$3 }' "${seq}_${safe_f}"s.output) #Prints gene name (column 2 from SQL query) and mutation (#col 3)
            sed -i "${i}s/.*/&,Sensitive,${sum_resistance},${mech}/" ${seq}.drug.table.txt
			#cat "$f"s.output >> ${seq}.drug.table.tertiary.txt
			i=$((i+1))
		fi
	else
		echo "no mechanism identified for $f sensitivity"
		sed -i "${i}s/.*/&,Resistant,${sum_resistance},No sensitivity detected/" ${seq}.drug.table.txt
		i=$((i+1))
	fi
done < <(grep -E "intrinsic|Intrinsic" ${seq}.drug.table.txt.backup | awk -F "," '{ print $3 }')

# --- CHANGED: operate on per-sample AbR output and drug table for second-line section ---
while read f; do 
    safe_f=$(echo "$f" | sed 's/[^A-Za-z0-9]/_/g')
    awk -F"|" -v f="$f" '$4~ f"r"' ${seq}.AbR_output.txt > "${seq}_${safe_f}r.output"
    awk -F"|" -v f="$f" '$4~ f"i"' ${seq}.AbR_output.txt > "${seq}_${safe_f}i.output"
    awk -F"|" -v f="$f" '$4~ f"s"' ${seq}.AbR_output.txt > "${seq}_${safe_f}s.output"
	#awk -F"|" -v f="$f" '$4~ f"r"' "${seq}"_resfinder_report.txt >> "$f"r.output
	#awk -F"|" -v f="$f" '$4~ f"i"' "${seq}"_resfinder_report.txt >> "$f"i.output
	
	##calc level of resistance - messy code here. TO DO - neaten and reduce redundancy 
	
	awk -F"|" 'BEGIN { OFS = "\n" } {print $4,$5}' "$f"r.output | awk -F"," '{
    for (i=1; i<=NF; i++) {
        a[NR,i] = $i
        }
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
    }' | sed 's/\t/\n/g' | grep "$f" -A1 > "$f"r.level_calc
	
	awk -F"|" 'BEGIN { OFS = "\n" } {print $4,$5}' "$f"i.output | awk -F"," '{
    for (i=1; i<=NF; i++) {
        a[NR,i] = $i
        }
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
    }' | sed 's/\t/\n/g' | grep "$f" -A1 > "$f"i.level_calc
	
	awk -F"|" 'BEGIN { OFS = "\n" } {print $4,$5}' "$f"s.output | awk -F"," '{
    for (i=1; i<=NF; i++) {
        a[NR,i] = $i
        }
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

    cat "${seq}_${safe_f}"{r,s,i}.level_calc > "${seq}_${safe_f}".level_calc.tmp

    sum_resistance=$(awk '{ for(i=1; i<=NF; i++) { if($i ~ /^-?[0-9]+$/) { sum += $i } } } END { print sum }' "${seq}_${safe_f}".level_calc.tmp)
	#echo "$sum_resistance"
	if [ ! -n "$sum_resistance" ]; then
	  sum_resistance=0
    fi
	
    grep -w "$f"r "${seq}_${safe_f}"r.output &> /dev/null
	status=$?
	if [[ "$status" -eq 0 ]]; then
		echo "found mechanism for $f resistance"
		length=$(wc -l "$f"r.output | awk '{print $1}' )
		if [[ "$length" -gt 1 ]]; then
			echo "found multiple determinants for $f resistance"
			#cat "$f"s.output >> ${seq}.drug.table.tertiary.txt
			sed -i "${i}s/.*/&,Resistant,${sum_resistance},Multiple determinants/" ${seq}.drug.table.txt
			i=$((i+1))
		else
			echo "found single mechanism for $f resistance" 
            mech=$(awk -F "|" '{ print $1,$2,$3 }' "${seq}_${safe_f}"r.output) #Prints gene name (column 2 from SQL query) and mutation (#col 3
            sed -i "${i}s/.*/&,Resistant,${sum_resistance},${mech}/" ${seq}.drug.table.txt
			#cat "$f"s.output >> ${seq}.drug.table.tertiary.txt
			i=$((i+1))
		fi
	else
		echo "no mechanism identified for $f resistance"
    sed -i "${i}s/.*/&,Sensitive,${sum_resistance},No resistance detected/" ${seq}.drug.table.txt
		i=$((i+1))
	fi
 done < <(grep -E "Second-line|second-line" ${seq}.drug.table.txt.backup | awk -F "," '{ print $3 }') 

#Looking for resistance
while read f; do
	awk -F"|" -v f="$f" '$4~ f"r"' ${seq}.AbR_output.txt > "$f"r.output
	awk -F"|" -v f="$f" '$4~ f"i"' ${seq}.AbR_output.txt > "$f"i.output
	awk -F"|" -v f="$f" '$4~ f"s"' ${seq}.AbR_output.txt > "$f"s.output
	#awk -F"|" -v f="$f" '$4~ f"r"' "${seq}"_resfinder_report.txt >> "$f"r.output
	
		##calc level of resistance - messy code here. TO DO - neaten and reduce redundancy 
	
	awk -F"|" 'BEGIN { OFS = "\n" } {print $4,$5}' "$f"r.output | awk -F"," '{
    for (i=1; i<=NF; i++) {
        a[NR,i] = $i
        }
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
    }' | sed 's/\t/\n/g' | grep "$f" -A1 > "$f"r.level_calc
	
	awk -F"|" 'BEGIN { OFS = "\n" } {print $4,$5}' "$f"i.output | awk -F"," '{
    for (i=1; i<=NF; i++) {
        a[NR,i] = $i
        }
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
    }' | sed 's/\t/\n/g' | grep "$f" -A1 > "$f"i.level_calc
	
	awk -F"|" 'BEGIN { OFS = "\n" } {print $4,$5}' "$f"s.output | awk -F"," '{
    for (i=1; i<=NF; i++) {
        a[NR,i] = $i
        }
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
	
	grep -w "$f"r "$f"r.output &> /dev/null #looks for full resistance
	status=$?
	if [[ "$status" -eq 0 ]]; then
		echo "found mechanism for $f resistance"
		length=$(wc -l "$f"r.output | awk '{print $1}' )
		if [[ "$length" -gt 1 ]]; then
			echo "found multiple determinants for $f resistance"
			#cat "$f"s.output >> ${seq}.drug.table.tertiary.txt
			sed -i "${i}s/.*/&,Resistant,${sum_resistance},Multiple determinants/" ${seq}.drug.table.txt
			i=$((i+1))
		else
			echo "found single mechanism for $f resistance" 
			mech=$(awk -F "|" '{ print $1,$2,$3 }' "$f"r.output) #Prints gene name (column 2 from SQL query) and mutation (#col 3
			sed -i "${i}s/.*/&,Resistant,${sum_resistance},${mech}/" ${seq}.drug.table.txt
			#cat "$f"s.output >> ${seq}.drug.table.tertiary.txt
			i=$((i+1))
		fi
	else
		echo "no mechanism identified for $f resistance, looking for intermediate resistance"
		grep -w "${f}"i "$f"i.output &> /dev/null
		status=$?
		if [[ "$status" -eq 0 ]]; then
			echo "found intermediate resistance mechanism for $f"
			length=$(wc -l "$f"i.output | awk '{print $1}' )
			if [[ "$length" -gt 1 ]]; then
				echo "found multiple determinants for intermediate $f resistance"
				sed -i "${i}s/.*/&,Intermediate,${sum_resistance},Multiple determinants/" ${seq}.drug.table.txt
				i=$((i+1))
			else
				echo "found single mechanism for intermediate $f resistance" 
				mech=$(awk -F "|" '{ print $1,$2,$3 }' "$f"i.output) #Prints gene name (column 2 from SQL query) and mutation (col 3
				sed -i "${i}s/.*/&,Intermediate,${sum_resistance},${mech}/" ${seq}.drug.table.txt
				i=$((i+1))
			fi
		else
			echo "no intermediate resistance found"
			sed -i "${i}s/.*/&,Sensitive,${sum_resistance},No resistance detected/" ${seq}.drug.table.txt
			i=$((i+1))
		fi
	fi
done < <(grep -E "tertiary|Tertiary" ${seq}.drug.table.txt.backup | awk -F "," '{ print $3 }')

#Looking for speciation and QC
echo "Looking for speciation markers"
awk -F"|" -v f="speciation" 'BEGIN{IGNORECASE=1} $4~ f' ${seq}.AbR_output.txt > ${seq}.speciation.output
grep -iw "speciation" ${seq}.speciation.output &> /dev/null
status=$?
if [[ "$status" -eq 0 ]]; then
	echo "Species specific marker is missing. Likely quality control issue" #TO DO This should go into final report
	echo "Interpret AMR profile with caution" #TO DO This should go into final report
fi	

# create patientDrugSusceptibilityData.csv (robust generator)
# ID refers to individual strains
out=${seq}.drug.table.txt
tmp=${seq}.drug.table.txt.tmp2
echo "ID,Class,Drug,Status,Level,Details" > "$tmp"

# If backup doesn't exist, fallback to the generated drug table
if [ ! -f ${seq}.drug.table.txt.backup ]; then
    cp ${seq}.drug.table.txt ${seq}.drug.table.txt.backup
fi

while IFS=',' read -r dclass ddrug dabbrev; do
    # skip empty lines or header-like lines
    if [[ -z "$dabbrev" || "$dclass" == "Drug_class" || "$dabbrev" == "Abbreviation" ]]; then
        continue
    fi
    abbrev=$(echo "$dabbrev" | sed 's/[^A-Za-z0-9]//g')
    # default values
    status="Sensitive"
    level=0
    details="No resistance detected"

    if [ -s ${seq}.AbR_output.final.txt ]; then
        # find any matching records (case-insensitive) where the 4th field contains the abbreviation + r/i/s
        # prefer Resistant > Intermediate > Sensitive
        r_match=$(awk -F"|" -v a="${abbrev}" 'BEGIN{IGNORECASE=1} tolower($4) ~ tolower(a)"r" {print; exit}' ${seq}.AbR_output.final.txt)
        i_match=$(awk -F"|" -v a="${abbrev}" 'BEGIN{IGNORECASE=1} tolower($4) ~ tolower(a)"i" {print; exit}' ${seq}.AbR_output.final.txt)
        s_match=$(awk -F"|" -v a="${abbrev}" 'BEGIN{IGNORECASE=1} tolower($4) ~ tolower(a)"s" {print; exit}' ${seq}.AbR_output.final.txt)

        if [ -n "$r_match" ]; then
            status="Resistant"
            # collect all mechanics for this abbrev
            details=$(awk -F"|" -v a="${abbrev}" 'BEGIN{IGNORECASE=1} tolower($4) ~ tolower(a)"r" {print $1" "$2" "$3}' ${seq}.AbR_output.final.txt | tr '\n' ';' | sed 's/;$//')
            level=$(awk -F"|" -v a="${abbrev}" 'BEGIN{IGNORECASE=1} { if(tolower($4) ~ tolower(a)"r" ){ g=$5; gsub(/[^0-9-]/,"",g); sum+=g } } END{print (sum+0)}' ${seq}.AbR_output.final.txt)
        elif [ -n "$i_match" ]; then
            status="Intermediate"
            details=$(awk -F"|" -v a="${abbrev}" 'BEGIN{IGNORECASE=1} tolower($4) ~ tolower(a)"i" {print $1" "$2" "$3}' ${seq}.AbR_output.final.txt | tr '\n' ';' | sed 's/;$//')
            level=$(awk -F"|" -v a="${abbrev}" 'BEGIN{IGNORECASE=1} { if(tolower($4) ~ tolower(a)"i" ){ g=$5; gsub(/[^0-9-]/,"",g); sum+=g } } END{print (sum+0)}' ${seq}.AbR_output.final.txt)
        elif [ -n "$s_match" ]; then
            status="Sensitive"
            details=$(awk -F"|" -v a="${abbrev}" 'BEGIN{IGNORECASE=1} tolower($4) ~ tolower(a)"s" {print $1" "$2" "$3}' ${seq}.AbR_output.final.txt | tr '\n' ';' | sed 's/;$//')
            level=$(awk -F"|" -v a="${abbrev}" 'BEGIN{IGNORECASE=1} { if(tolower($4) ~ tolower(a)"s" ){ g=$5; gsub(/[^0-9-]/,"",g); sum+=g } } END{print (sum+0)}' ${seq}.AbR_output.final.txt)
        fi
    fi

    # Write single line for this drug
    echo "${seq},${dclass},${ddrug},${status},${level},${details}" >> "$tmp"
done < <(awk -F"," 'NF>=3{print $1","$2","$3}' ${seq}.drug.table.txt.backup)

# move into place
mv "$tmp" ${seq}.drug.table.txt
cp ${seq}.drug.table.txt ${seq}_patientDrugSusceptibilityData.csv

if [ -s ${seq}.AbR_output.final.txt ]; then
    cp ${seq}.AbR_output.final.txt "${seq}.AbR_output.final.txt"
else	
    echo "No antibiotic resistance identified in $seq" >> "${seq}.AbR_output.final.txt"
fi

# move per-sample patientDrugSusceptibilityData into the expected final filename
mv ${seq}_patientDrugSusceptibilityData.csv patientDrugSusceptibilityData.csv

# cleanup per-sample metadata copy
rm -f ${seq}_patientMetaData.csv

exit 0
