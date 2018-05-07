#!/bin/bash

reset

echo "****************************************"
echo "*    Pipeline to annotate variants     *"
echo "*                                      *"
echo "* Neurodevelopmental disorders project *"
echo "*                                      *"
echo "*             Author: Luigi Chiricosta *"
echo "****************************************"
echo
echo
read -r -p "Do you want run the complete pipeline? [Y/N] " completePipeline

if [ "$1" != "" ]; then
	project=$1
else
	project="ASD_Project"
	echo "I use the default project name."
fi

echo "Project name: $project"
echo

#initialize path variables
home="$HOME"
desktop="$home/Desktop"
annovarToolPath="$desktop/Software/annovar"
workingPath="$desktop/Pipeline_Luigi_Chiricosta"
#VCFInputPath="$workingPath/VCFfiles/12-13_04_2018"
VCFInputPath="$desktop/Exoms_data/SNP/Vcf"
path="$workingPath/$project"
annovarOutput="$path/annovarOutput"
mergeFileName="merged_output_vcf.vcf.gz"
outputAnnovarFileName="merg_out"
reducedGnomADRelease="$desktop/Pipeline_Luigi_Chiricosta/GnomAD releases/reduced_gnomad_release.txt" #extra file
#done

if [ $(dpkg -l | grep -E '^ii' | grep -w realpath | wc -l) -eq 0 ]; then
        echo "Realpath command is not in the system... installing..."
        sudo apt-get install realpath
fi

echo "> Project folder path: "$(realpath $path)

if [ "$completePipeline" = "N" ]; then
        echo
        echo
        read -r -p "Do you want to prepare the input vcf? [Y/N] " pipelineNextStep
else
	echo "Complete option confirmed."
	echo
fi

if [[ "$completePipeline" != "N" ]] || [[ "$pipelineNextStep" != "N" ]]; then
	if [ $(dpkg -l | grep -E '^ii' | grep -w tabix | wc -l) -eq 0 ]; then
		echo "Tabix command is not in the system... installing..."
		sudo apt-get install tabix
	fi

	if [[ -e $path ]] ; then
		i=1
	    	while [[ -e $path-$i ]] ; do
	        	let i++
		done
	
		path=$path-$i
	        echo "> Redefined project folder path: "$(realpath $path)
	fi
	
	mkdir "$path"
	
	annovarOutput="$path/annovarOutput" #rewrite annovarOutput if create a new project folder
        echo "> Annovar Output: $annovarOutput <"
        mkdir $annovarOutput

	#TODO use input folder for vcf files
	echo "> VCF input files in: "$(realpath $VCFInputPath)
	
	#if vcf does not have its tbi, do it
	for vcf in $VCFInputPath/*.vcf.gz; do
		name=$(basename "$vcf")
		if [[ ! -f "$vcf.tbi" ]]; then
			echo -e "$name is not alreay indexed, i am doing it... \c"
			tabix -p vcf "$vcf"
			echo "done."
		fi
	done
	
	echo
	echo -e "Merging the vcf files... \c"
	vcf-merge $VCFInputPath/*.vcf.gz 2>/dev/null | bgzip -c > $annovarOutput/$mergeFileName
	echo "done."
	echo "Annoting with annovar... "
	perl $annovarToolPath/table_annovar.pl $annovarOutput/$mergeFileName $annovarToolPath/humandb/ -buildver hg19 -vcfinput -out $annovarOutput/$outputAnnovarFileName -remove -protocol refGene,dbnsfp33a,intervar_20170202,clinvar_20170905,avsnp150,gnomad_genome,gnomad_exome,cosmic70 -operation g,f,f,f,f,f,f,f -nastring .
	echo
	echo "Completed."
	echo
fi
	
if [ "$completePipeline" = "N" ]; then
	echo
	echo
	read -r -p "Do you want to build the new table? [Y/N] " pipelineNextStep
fi

if [[ "$completePipeline" != "N" ]] || [[ "$pipelineNextStep" != "N" ]]; then
	echo "Building new table..."
	multiannoPath="$annovarOutput/$outputAnnovarFileName.hg19_multianno"
	
	#Catch the sort of the columns
	SPECIAL_FIELD="PATIENTS_LIST"
	vcfInfo=$(grep -w "#CHROM" "$multiannoPath.vcf")
	originalCompleteSort=$(head -n1 "$multiannoPath.txt" | sed "s/Otherinfo/Otherinfo1\tOtherinfo2\tOtherinfo3\t$vcfInfo\t/g")
	newSort=$(echo -e "Chr\tStart\tEnd\tRef\tAlt\tFunc.refGene\tGene.refGene\tGeneDetail.refGene\tExonicFunc.refGene\tAAChange.refGene\tills_in_db\tpatients_id_list\tavsnp150\tcosmic70\tgnomAD_exome_ALL\tgnomAD_genome_ALL\tAC_gnomAD\tHom_gnomAD\tHemi_gnomAD\tInterVar(automated)\tCLINSIG\tCLNACC\tSIFT_pred\tPolyphen2_HDIV_pred\tPolyphen2_HVAR_pred\tLRT_pred\tMutationTaster_pred\tMutationAssessor_pred\tFATHMM_pred\tPROVEAN_pred\tMetaSVM_pred\tMetaLR_pred\tM-CAP_pred\tCADD_phred\tfathmm-MKL_coding_pred\tGERP++_RS\tInterpro_domain\tINFO\tFORMAT\t$SPECIAL_FIELD") #insert column PATIENTS_LIST to have the list of patients starting from that point
	firstColumn="Chr\tPOS\tREF\tALT";
	
	#Retrevial the patients
	patientsList=$(echo -e "$vcfInfo" | awk -F '\t' '{for(i=10;i<=NF;i++) printf ("%s%s", $i, (i!=NF) ? "\t" : "")}')
	patients=$(echo -e "$patientsList" | awk -F '\t' '{print split($0, patients, "\t")}')
	echo "Found $patients patients."
	
	#Adjust the list if you want the patients too
	patientStartingList=${newSort#*$SPECIAL_FIELD}
	patientStartingPosition=$(((${#newSort}-${#patientStartingList}-${#SPECIAL_FIELD})+1))
	if [[ $patientStartingPosition -gt -1 ]]; then
		startingFieldPatientsList=$(echo -e "${newSort:0:$patientStartingPosition}" | awk -F "\t" '{print NF}')
		newSort=$(echo -e "$newSort" | sed "s/$SPECIAL_FIELD/$patientsList/g")
	fi
	
	columnFORMAT=$(grep -e ^[^#] "$multiannoPath.vcf" | awk -F '\t' '{valuesFORMATLength=split($9, valuesFORMAT, ":"); for(i=1;i<=valuesFORMATLength;i++) print valuesFORMAT[i]}' | sort | uniq | awk 'BEGIN{FORMAT=""}{if(FORMAT=="")FORMAT=$0; else FORMAT=FORMAT":"$0;}END{print FORMAT}')
	
	cp "$multiannoPath.txt" "$multiannoPath.complete.txt"
	
	#add Extra columns
	if [[ "$newSort" == *"AC_gnomAD"* || "$newSort" == *"Hom_gnomAD"* || "$newSort" == *"Hemi_gnomAD"* ]]; then
		echo "Found extra GnomAD columns"
	
		tmp="$multiannoPath.txt.tmp"
	
		awk -F '\t' -v originalCompleteSort="$originalCompleteSort" -v firstColumn="$firstColumn" '
		BEGIN{
	        	#elaborate the names row for each column
		        idxLength=split(firstColumn, idxArray, "\t");
		
		        #elaborate the idx column (first one)
		        for(i=1; i<=idxLength; i++)
		                printf("%s%s", idxArray[i], (i!=idxLength) ? ":" : "\t");
	        	
		        originalSortSplittedLength = split(originalCompleteSort, originalSortSplitted, "\t");
	        	for(i=1;i<=originalSortSplittedLength; i++)
	                	originalSortDictionary[originalSortSplitted[i]]=i
		}
		{
		        if(NR==1)
		                print $0
	        
		        #starting with a new row, reset the variables
		        if(NR!=1)
		        {
		                #elaborate idx column (first one)
		                for(i=1; i<=idxLength; i++)
		                        printf("%s%s", $originalSortDictionary[idxArray[i]], (i!=idxLength) ? ":" : "\t"$0"\n");
		
		        }
		}' "$multiannoPath.txt" > "$tmp"
		rm "$multiannoPath.complete.txt"
		mv "$tmp" "$multiannoPath.complete.txt"
		
		format="0"
	        var1=$(echo -e "$originalCompleteSort" | awk -F "\t" '{print NF+1}'); 
	        for i in $(seq 2 "$var1"); do
	                format=$format",1.$i"; 
	        done
	
	        var2=$(head -n1 "$reducedGnomADRelease" | awk -F "\t" '{print NF}'); 
	        for i in $(seq 2 "$var2"); do
	                format=$format",2.$i";
	        done
		
		echo -e "I am joining it... \c"
		
		header1=$(head -n1 "$multiannoPath.complete.txt")
		header2=$(head -n1 "$reducedGnomADRelease" | cut -f2- -d$'\t')
		echo -e "$header1\t$header2" | awk -F '\t' '{for(i=2; i<=NF; i++) printf("%s%s", $i, (i!=NF) ? "\t" : "\n");}' > "$tmp"
		
		join -t $'\t' -a1 -e "-" -o "$format" <(awk -F '\t' '{if(NR!=1)print $0}' "$multiannoPath.complete.txt" | sort -k1,1) <(awk -F '\t' '{if(NR!=1){print $0}}' "$reducedGnomADRelease" | sort -k1,1) | awk -F '\t' '{for(i=2; i<=NF; i++) printf("%s%s", $i, (i!=NF) ? "\t" : "\n");}' | vcf-sort -c >> "$tmp" 2>/dev/null
		
	
		if [[ "$newSort" == *"AC_gnomAD"* ]]; then
			originalCompleteSort="$originalCompleteSort\tAC_gnomAD"
		fi
		if [[ "$newSort" == *"Hom_gnomAD"* ]]; then
	                originalCompleteSort="$originalCompleteSort\tHom_gnomAD"
	        fi
		if [[ "$newSort" == *"Hemi_gnomAD"* ]]; then
	                originalCompleteSort="$originalCompleteSort\tHemi_gnomAD"
	        fi
		echo "done."
		
		rm "$multiannoPath.complete.txt"
		mv "$tmp" "$multiannoPath.complete.txt"
	fi
	
	echo "Building the new table... \c"
	
	awk -F '\t' -v originalCompleteSort="$originalCompleteSort" -v newSort="$newSort" -v firstColumn="$firstColumn" -v patients="$patients" -v firstPatientsField="$startingFieldPatientsList" -v columnFORMAT="$columnFORMAT" 'BEGIN{
		#elaborate the names row for each column
		#elaborate the idx column (first one)
	        idxLength=split(firstColumn, idxArray, "\t");
	        for(i=1; i<=idxLength; i++)
	                printf("%s%s", idxArray[i], (i!=idxLength) ? ":" : "\t");
	
		originalSortSplittedLength = split(originalCompleteSort, originalSortSplitted, "\t");
		for(i=1;i<=originalSortSplittedLength; i++)
			originalSortDictionary[originalSortSplitted[i]]=i
		newSortSplittedLength = split(newSort, newSortSplitted, "\t");
		for(i=1; i<=newSortSplittedLength; i++)
	        	printf("%s%s", newSortSplitted[i], (i!=newSortSplittedLength) ? "\t" : "\n");
		lastPatientsField=firstPatientsField+patients
	}
	{
		#starting with a new row, reset the variables
		if(NR!=1)
		{	
			ills=0
			illsList=""
	
			#elaborate idx column (first one)
	                for(i=1; i<=idxLength; i++)
	                        printf("%s%s", $originalSortDictionary[idxArray[i]], (i!=idxLength) ? ":" : "\t");
	
			#elaborate all the other columns
			finalArray=""
			for(i=1; i<=newSortSplittedLength; i++)
		        {
				#check if the actual value is a special column or not
		                if(originalSortDictionary[newSortSplitted[i]] == "") #the actual column is special
				{
					if(newSortSplitted[i] == "ills_in_db")
						finalArray=finalArray"ills_in_db\t"
					else if(newSortSplitted[i] == "patients_id_list")
	                                        finalArray=finalArray"patients_id_list\t"
					
					else
					{
						finalArray=finalArray"COLUMN_ERROR\t"
			                        printf "COLUMN_ERROR("newSortSplitted[i]")\t"
					}
		                }
				else #the actual column is not special
				{
					if(firstPatientsField == -1 || i<firstPatientsField || i>=lastPatientsField)
					{
						if(newSortSplitted[i] == "FORMAT")
							finalArray=finalArray""columnFORMAT"\t"
						else
							finalArray=finalArray""$originalSortDictionary[newSortSplitted[i]]"\t"
					}
					else
					{
						#check if is a patient column and its status
						if(firstPatientsField != -1) #no one patient could be present
						{
							if(i>=firstPatientsField && i<lastPatientsField)
							{
								#print "patient ("i"): "$originalSortDictionary[newSortSplitted[i]]
								if($originalSortDictionary[newSortSplitted[i]]!=".")
								{
									#print $originalSortDictionary["FORMAT"]
									columnFORMATSplittedLength = split($originalSortDictionary["FORMAT"], columnFORMATSplitted, ":")
								        for(j=1;j<=columnFORMATSplittedLength; j++)
										columnFORMATDictionary[columnFORMATSplitted[j]]=j
									split($originalSortDictionary[newSortSplitted[i]], valuesFORMATSplitted, ":")
								
									#recostruct the new FORMAT column sort
									newFORMATColumn = ""
									columnFORMATSplittedLength = split(columnFORMAT, columnFORMATSplitted, ":")
									for(j=1; j<=columnFORMATSplittedLength;j++)
									{
										newFORMATColumn=newFORMATColumn""valuesFORMATSplitted[columnFORMATDictionary[columnFORMATSplitted[j]]]
										if(j<columnFORMATSplittedLength)
											newFORMATColumn=newFORMATColumn":"
									}
									finalArray=finalArray""newFORMATColumn"\t"
	
									#specific format details for patients with variant
									GQ = valuesFORMATSplitted[columnFORMATDictionary["GQ"]]
									DP = valuesFORMATSplitted[columnFORMATDictionary["DP"]]
									AF = valuesFORMATSplitted[columnFORMATDictionary["AF"]]
									#print newSortSplitted[i]" : "GQ"("columnFORMATDictionary["GQ"]") - "DP"("columnFORMATDictionary["DP"]") -"AF"("columnFORMATDictionary["AF"]")"
			
	        		                                        ills=ills+1
	                		                                illsList=illsList""newSortSplitted[i]"(GQ: "GQ", DP: "DP", AF: "AF"); "
								}
								else
									finalArray=finalArray".\t"
							}
			
							#filter to not visualize the patients of common variants
							if(i==lastPatientsField-1 && ills>5)
								illsList="-"
						}
					}
				}
			}
			
			#replace the special columns
			gsub("ills_in_db", ills"/"patients, finalArray)
			gsub("patients_id_list", illsList, finalArray)
			
			#print the final record
			print finalArray
		}
	}' "$multiannoPath.complete.txt" | awk -F "\t" '{if(NR==1){gsub("CLNACC", "ClinVar_ID", $0); gsub("AC_gnomAD", "AC", $0); gsub("Hom_gnomAD", "Hom", $0); gsub("Hemi_gnomAD", "Hem", $0)}; print $0}' | uniq > "$path/$project.xls" #uniq command is necessary because in txt annovar output some rows are duplicated maybe due to more than one original reference allele

	echo "done."
	echo
	echo "> Final output in: $path/$project.xls"
	echo
	echo "Completed."
fi

if [ "$completePipeline" = "N" ]; then
        echo
        echo
        read -r -p "Do you want to attach the transcript information? [Y/N] " pipelineNextStep
fi

if [[ "$completePipeline" != "N" ]] || [[ "$pipelineNextStep" != "N" ]]; then
        echo "Attaching transcript information..."
	
	if [ $(dpkg -l | grep -E '^ii' | grep -w curl | wc -l) -eq 0 ]; then
        	echo "Realpath command is not in the system... installing..."
       		sudo apt-get install curl
	fi
	
	tmp="$path/$project.plus.xls"; 
	> $tmp;
 	
	lines=$(cat "$path/$project.xls" | wc -l)
	actualLine=0
	while read in; do 
	        if [[ "$in" == *"Chr:POS:REF:ALT"* ]]; then 
			echo -e "Mutalyzer\t$in" >> "$tmp"
		else
			actualLine=$((actualLine+1))
			id=$(awk -F '\t' '{split($1, id, ":"); multiplesLength=split(id[4], multiples, ","); if(multiplesLength==1) print id[1]":g."id[2]""id[3]">"id[4]; else print ""}' <<< "$in"); 
			if [ "$id" = "" ]; then
				transcript="-"
			else
		                transcript=$(curl 'https://mutalyzer.nl/json/numberConversion?build=hg19&variant='$id'' 2>/dev/null | sed 's/\[//g' | sed 's/\]//g' | sed 's/\"//g')
		fi

	                echo -e "$transcript\t$in" >> "$tmp"; 
			echo -ne "Progress: $actualLine/$lines\r"
	        fi; 
	done < "$path/$project.xls"
	
	echo
	echo "done."
fi

echo
echo
echo "Pipeline finished."
