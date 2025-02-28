#!/bin/bash

function p04_krona() {
  if [ -f "$path4/.done" ]; then
    echo "Krona process already finished, skipped"
    return
  fi

  local BLAST=$path3/${OUT}_blast_hit.txt
  local TR=`cat $path1/${OUT}_log_p1.txt | awk 'NR==3 {print $4}' | tr -d ","`
  local OUTR=$path4/${OUT}_report.csv
  local OUTR2=$path4/${OUT}_bacteria_report.csv

  ##Basic data build
  mkdir -p $path4/${OUT}_CSV
  mkdir -p $path4/${OUT}_CSV2
  mkdir -p $path4/${OUT}_fRep
  mkdir -p $path4/${OUT}_bRep

  ktClassifyBLAST -p -s ${BLAST} -o $path4/${OUT}_CSV/${OUT}_blast.csv
  cut -f 1,2 $path4/${OUT}_CSV/${OUT}_blast.csv > $path4/${OUT}_CSV/${OUT}_merge.csv

  ktImportTaxonomy \
    -m 1 -s 3 \
    $path4/${OUT}_CSV/${OUT}_merge.csv \
    -o $path4/${OUT}_pathogen_final.html 

  cat $path4/${OUT}_CSV/${OUT}_merge.csv \
    | sed -e '1d' \
    | awk '{ tax[$2] += $1; } END { for (i in tax) print tax[i] "\t" i; }' \
    | sort -k 1nr,1 \
    > $path4/${OUT}_CSV/${OUT}_merge_taxid.csv

  taxonkit lineage -i 2 $path4/${OUT}_CSV/${OUT}_merge.csv \
    > $path4/${OUT}_CSV/${OUT}_merge_tax.csv

  cat $path4/${OUT}_CSV/${OUT}_merge_tax.csv \
    | awk -v "r=${TR}" '{print $1/r*1000000"\t"$0}' \
    | cut -f 1,3- \
    > $path4/${OUT}_CSV/${OUT}_merge_tax2.csv

  #Create CSV file for each taxonomic rank
  #Kingdom
  rank=k
  taxonkit reformat -j 8 -i 3 -f {$rank} $path4/${OUT}_CSV/${OUT}_merge_tax2.csv \
    | cut -f 1,4- | sed 's/ /_/g' | awk '{ tax[$2] += $1; } END { for (i in tax) print i "\t" tax[i]; }' \
    | sort -k 2nr,2 \
    | (s=$(cat)&&awk 'NR==FNR{a=a+$2;next}{p=$2/a;print $0"\t"p}' <(echo -e "$s") <(echo -e "$s")) \
    > $path4/${OUT}_CSV2/${OUT}_mf${rank^^}.csv

  taxonkit reformat -j 8 -i 3 -f {$rank} $path4/${OUT}_CSV/${OUT}_merge_tax.csv \
    | cut -f 1,4- | sed 's/ /_/g' | awk '{ tax[$2] += $1; } END { for (i in tax) print i "\t" tax[i]; }' \
    | sort -k 2nr,2 \
    > $path4/${OUT}_CSV2/${OUT}_rf${rank^^}.csv

  paste $path4/${OUT}_CSV2/${OUT}_rf${rank^^}.csv $path4/${OUT}_CSV2/${OUT}_mf${rank^^}.csv \
    | cut -f 1,2,4,5 \
    > $path4/${OUT}_fRep/${OUT}_f${rank^^}1.csv

  #Phylum, Family, Genus, Species
  for rank in p f g s ; do
    taxonkit reformat -j 8 -i 3 -f {$rank} $path4/${OUT}_CSV/${OUT}_merge_tax2.csv \
      | cut -f 1,4- | sed 's/ /_/g' | awk '{ tax[$2] += $1; } END { for (i in tax) print i "\t" tax[i]; }' \
      | sort -k 2nr,2 \
      | grep -v -e '^\s' \
      | (s=$(cat)&&awk 'NR==FNR{a=a+$2;next}{p=$2/a;print $0"\t"p}' <(echo -e "$s") <(echo -e "$s")) \
      > $path4/${OUT}_CSV2/${OUT}_mf${rank^^}.csv

    taxonkit reformat -j 8 -i 3 -f {$rank} $path4/${OUT}_CSV/${OUT}_merge_tax.csv \
      | cut -f 1,4- | sed 's/ /_/g' | awk '{ tax[$2] += $1; } END { for (i in tax) print i "\t" tax[i]; }' \
      | sort -k 2nr,2 \
      | grep -v -e '^\s' \
      > $path4/${OUT}_CSV2/${OUT}_rf${rank^^}.csv

    paste $path4/${OUT}_CSV2/${OUT}_rf${rank^^}.csv $path4/${OUT}_CSV2/${OUT}_mf${rank^^}.csv \
      | cut -f 1,2,4,5 \
      > $path4/${OUT}_fRep/${OUT}_f${rank^^}1.csv
  done

  ##Bacteria, Viral & Fungal read extraction
  #Bacteria
  grep Bacteria $path4/${OUT}_CSV/${OUT}_merge_tax2.csv > $path4/${OUT}_CSV/${OUT}_merge_taxb.csv
  grep Bacteria $path4/${OUT}_CSV/${OUT}_merge_tax.csv > $path4/${OUT}_CSV/${OUT}_merge_taxrb.csv

  #Bac_species, Bac_genus, Bac_family
  for rank in s g f ; do
    taxonkit reformat -j 8 -i 3 -f {$rank} $path4/${OUT}_CSV/${OUT}_merge_taxb.csv \
      | cut -f 1,4- | sed 's/ /_/g' | awk '{ tax[$2] += $1; } END { for (i in tax) print i "\t" tax[i]; }' \
      | sort -k 2nr,2 \
      | grep -v -e '^\s' \
      | (s=$(cat)&&awk 'NR==FNR{a=a+$2;next}{p=$2/a;print $0"\t"p}' <(echo -e "$s") <(echo -e "$s")) \
      > $path4/${OUT}_CSV2/${OUT}_mBacteria${rank^^}.csv

    taxonkit reformat -j 8 -i 3 -f {$rank} $path4/${OUT}_CSV/${OUT}_merge_taxrb.csv \
      | cut -f 1,4- | sed 's/ /_/g' | awk '{ tax[$2] += $1; } END { for (i in tax) print i "\t" tax[i]; }' \
      | sort -k 2nr,2 \
      | grep -v -e '^\s' \
      > $path4/${OUT}_CSV2/${OUT}_rBacteria${rank^^}.csv

    paste $path4/${OUT}_CSV2/${OUT}_rBacteria${rank^^}.csv $path4/${OUT}_CSV2/${OUT}_mBacteria${rank^^}.csv \
      | cut -f 1,2,4,5 \
      > $path4/${OUT}_fRep/${OUT}_Bacteria${rank^^}1.csv
  done

  #Virus
  grep Virus $path4/${OUT}_CSV/${OUT}_merge_tax2.csv > $path4/${OUT}_CSV/${OUT}_merge_taxv.csv
  grep Virus $path4/${OUT}_CSV/${OUT}_merge_tax.csv > $path4/${OUT}_CSV/${OUT}_merge_taxrv.csv

  taxonkit reformat -j 8 -i 3 -f {s} $path4/${OUT}_CSV/${OUT}_merge_taxv.csv \
    | cut -f 1,4- | sed 's/ /_/g' | awk '{ tax[$2] += $1; } END { for (i in tax) print i "\t" tax[i]; }' \
    | sort -k 2nr,2 \
    | grep -v -e '^\s' \
    | (s=$(cat)&&awk 'NR==FNR{a=a+$2;next}{p=$2/a;print $0"\t"p}' <(echo -e "$s") <(echo -e "$s")) \
    > $path4/${OUT}_CSV2/${OUT}_mVirus.csv

  taxonkit reformat -j 8 -i 3 -f {s} $path4/${OUT}_CSV/${OUT}_merge_taxrv.csv \
    | cut -f 1,4- | sed 's/ /_/g' | awk '{ tax[$2] += $1; } END { for (i in tax) print i "\t" tax[i]; }' \
    | sort -k 2nr,2 \
    | grep -v -e '^\s' \
    > $path4/${OUT}_CSV2/${OUT}_rVirus.csv

  paste $path4/${OUT}_CSV2/${OUT}_rVirus.csv $path4/${OUT}_CSV2/${OUT}_mVirus.csv \
    | cut -f 1,2,4,5 \
    > $path4/${OUT}_fRep/${OUT}_Virus1.csv

  #Fungi
  grep Fungi $path4/${OUT}_CSV/${OUT}_merge_tax2.csv > $path4/${OUT}_CSV/${OUT}_merge_taxf.csv
  grep Fungi $path4/${OUT}_CSV/${OUT}_merge_tax.csv > $path4/${OUT}_CSV/${OUT}_merge_taxrf.csv

  taxonkit reformat -j 8 -i 3 -f {s} $path4/${OUT}_CSV/${OUT}_merge_taxf.csv \
    | cut -f 1,4- | sed 's/ /_/g' | awk '{ tax[$2] += $1; } END { for (i in tax) print i "\t" tax[i]; }' \
    | sort -k 2nr,2 \
    | grep -v -e '^\s' \
    | (s=$(cat)&&awk 'NR==FNR{a=a+$2;next}{p=$2/a;print $0"\t"p}' <(echo -e "$s") <(echo -e "$s")) \
    > $path4/${OUT}_CSV2/${OUT}_mFungi.csv

  taxonkit reformat -j 8 -i 3 -f {s} $path4/${OUT}_CSV/${OUT}_merge_taxrf.csv \
    | cut -f 1,4- | sed 's/ /_/g' | awk '{ tax[$2] += $1; } END { for (i in tax) print i "\t" tax[i]; }' \
    | sort -k 2nr,2 \
    | grep -v -e '^\s' \
    > $path4/${OUT}_CSV2/${OUT}_rFungi.csv

  paste $path4/${OUT}_CSV2/${OUT}_rFungi.csv $path4/${OUT}_CSV2/${OUT}_mFungi.csv \
    | cut -f 1,2,4,5 \
    > $path4/${OUT}_fRep/${OUT}_Fungi1.csv

  #Calculate parameters(RIV,DiversityIndex D,H)
  PR2=` awk '{s += $2} END {print s}' $path4/${OUT}_fRep/${OUT}_fF1.csv`
  PR3=` awk '{s += $2} END {print s}' $path4/${OUT}_fRep/${OUT}_BacteriaF1.csv`

  PRstd2=` awk '{s += $2} END {print s}' $path4/${OUT}_CSV2/${OUT}_mfK.csv`
  PRstd3=` awk '{s += $2} END {print s}' $path4/${OUT}_CSV2/${OUT}_mBacteriaF.csv`

  FF1=`head -n3 $path4/${OUT}_fRep/${OUT}_fF1.csv`
  FG1=`head -n3 $path4/${OUT}_fRep/${OUT}_fG1.csv`
  FS1=`head -n3 $path4/${OUT}_fRep/${OUT}_fS1.csv`

  FSV=`head -n3 $path4/${OUT}_fRep/${OUT}_Virus1.csv`
  FSF=`head -n3 $path4/${OUT}_fRep/${OUT}_Fungi1.csv`

  FSB=`head -n3 $path4/${OUT}_fRep/${OUT}_BacteriaS1.csv`
  FSBG=`head -n3 $path4/${OUT}_fRep/${OUT}_BacteriaG1.csv`
  FSBF=`head -n3 $path4/${OUT}_fRep/${OUT}_BacteriaF1.csv`
  BF1=`head -n3 $path4/${OUT}_fRep/${OUT}_BacteriaF1.csv`
  BG1=`head -n3 $path4/${OUT}_fRep/${OUT}_BacteriaG1.csv`
  BS1=`head -n3 $path4/${OUT}_fRep/${OUT}_BacteriaS1.csv`

  #total Diversity
  FDF=`cut -f3 $path4/${OUT}_CSV2/${OUT}_mfF.csv | perl -nle '$D += $_*$_; print 1-$D if(eof);'`
  FDG=`cut -f3 $path4/${OUT}_CSV2/${OUT}_mfG.csv | perl -nle '$D += $_*$_; print 1-$D if(eof);'`
  FDS=`cut -f3 $path4/${OUT}_CSV2/${OUT}_mfS.csv | perl -nle '$D += $_*$_; print 1-$D if(eof);'`

  FHF=`cut -f3 $path4/${OUT}_CSV2/${OUT}_mfF.csv | perl -nle '$H += $_*log($_); print $H*(-1) if(eof);'`
  FHG=`cut -f3 $path4/${OUT}_CSV2/${OUT}_mfG.csv | perl -nle '$H += $_*log($_); print $H*(-1) if(eof);'`
  FHS=`cut -f3 $path4/${OUT}_CSV2/${OUT}_mfS.csv | perl -nle '$H += $_*log($_); print $H*(-1) if(eof);'`

  #bacteria Diversity
  BDF=`cut -f3 $path4/${OUT}_CSV2/${OUT}_mBacteriaF.csv | perl -nle '$D += $_*$_; print 1-$D if(eof);'`
  BDG=`cut -f3 $path4/${OUT}_CSV2/${OUT}_mBacteriaG.csv | perl -nle '$D += $_*$_; print 1-$D if(eof);'`
  BDS=`cut -f3 $path4/${OUT}_CSV2/${OUT}_mBacteriaS.csv | perl -nle '$D += $_*$_; print 1-$D if(eof);'`

  BHF=`cut -f3 $path4/${OUT}_CSV2/${OUT}_mBacteriaF.csv | perl -nle '$H += $_*log($_); print $H*(-1) if(eof);'`
  BHG=`cut -f3 $path4/${OUT}_CSV2/${OUT}_mBacteriaG.csv | perl -nle '$H += $_*log($_); print $H*(-1) if(eof);'`
  BHS=`cut -f3 $path4/${OUT}_CSV2/${OUT}_mBacteriaS.csv | perl -nle '$H += $_*log($_); print $H*(-1) if(eof);'`

  #total report###
  echo -e "[${OUT}_REPORT]" >> $OUTR
  echo -e "Total_read(reads)\t${TR}" >> $OUTR
  echo -e "Pathogen_read(reads)\t${PR2}" >> $OUTR
  echo -e "Pathogen_read(RPM)\t${PRstd2}" >> $OUTR

  echo -e "\n>Diversity_Index" >> $OUTR
  echo -e "Simpson_Index(family)\t${FDF}" >> $OUTR
  echo -e "Simpson_Index(genus)\t${FDG}" >> $OUTR
  echo -e "Simpson_Index(species)\t${FDS}" >> $OUTR
  echo -e "\nShannon_Index(family)\t${FHF}" >> $OUTR
  echo -e "Shannon_Index(genus)\t${FHG}" >> $OUTR
  echo -e "Shannon_Index(species)\t${FHS}" >> $OUTR
  echo -e "\n>Top_3hit_pathogen(family) Name,reads,RPM,RA\n${FF1}" >> $OUTR
  echo -e "\n>Top_3hit_pathogen(genus) Name,reads,RPM,RA\n${FG1}" >> $OUTR
  echo -e "\n>Top_3hit_pathogen(species) Name,reads,RPM,RA\n${FS1}" >> $OUTR

  #Bacterial report
  echo -e "\n>Top_3hit_Bacteria(species) Name,reads,RPM,RA" >> $OUTR
  VC=`cat $path4/${OUT}_fRep/${OUT}_BacteriaS1.csv | wc -l`
  if test ${VC} -eq 0 ; then
     echo -e "No Bacteria could be found." >> $OUTR
  else
     echo -e "${FSB}" >> $OUTR
  fi

  #Viral report
  echo -e "\n>Top_3hit_Virus(species) Name,reads,RPM,RA" >> $OUTR
  VC=`cat $path4/${OUT}_fRep/${OUT}_Virus1.csv | wc -l`
  if test ${VC} -eq 0 ; then
     echo -e "No Virus could be found." >> $OUTR
  else
     echo -e "${FSV}" >> $OUTR
  fi

  #Fungi report
  echo -e "\n>Top_3hit_Fungi(species) Name,reads,RPM,RA" >> $OUTR
  VC=`cat $path4/${OUT}_fRep/${OUT}_Fungi1.csv | wc -l`
  if test ${VC} -eq 0 ; then
     echo -e "No Fungi could be found." >> $OUTR
  else
     echo -e "${FSF}" >> $OUTR
  fi

  #Bacterial mini report
  echo -e "[${OUT}_BACTERIA_REPORT]" >> $OUTR2
  echo -e "Total_read(reads)\t${TR}" >> $OUTR2
  echo -e "Bacterial_read(reads)\t${PR3}" >> $OUTR2
  echo -e "Bacterial_read(RPM)\t${PRstd3}" >> $OUTR2

  echo -e "\n>Diversity_Index" >> $OUTR2
  echo -e "Simpson_Index(family)\t${BDF}" >> $OUTR2
  echo -e "Simpson_Index(genus)\t${BDG}" >> $OUTR2
  echo -e "Simpson_Index(species)\t${BDS}" >> $OUTR2
  echo -e "\nShannon_Index(family)\t${BHF}" >> $OUTR2
  echo -e "Shannon_Index(genus)\t${BHG}" >> $OUTR2
  echo -e "Shannon_Index(species)\t${BHS}" >> $OUTR2
  echo -e "\n>Top_3hit_bacteria(family) Name,reads,RPM,RA\n${BF1}" >> $OUTR2
  echo -e "\n>Top_3hit_bacteria(genus) Name,reads,RPM,RA\n${BG1}" >> $OUTR2
  echo -e "\n>Top_3hit_bacteria(species) Name,reads,RPM,RA\n${BS1}" >> $OUTR2

  for i in K P F G S ; do
    echo "Microorganisms,reads,RPM,RA" > $path4/${OUT}_fRep/${OUT}_f$i.csv
    cat $path4/${OUT}_fRep/${OUT}_f${i}1.csv >> $path4/${OUT}_fRep/${OUT}_f$i.csv
    rm -f $path4/${OUT}_fRep/${OUT}_f${i}1.csv
  done
  for i in F G S ; do
    echo "Microorganisms,reads,RPM,RA" > $path4/${OUT}_bRep/${OUT}_Bacteria_f$i.csv
    cat $path4/${OUT}_fRep/${OUT}_Bacteria${i}1.csv >> $path4/${OUT}_bRep/${OUT}_Bacteria_f$i.csv
    rm -f $path4/${OUT}_fRep/${OUT}_Bacteria${i}1.csv
  done
  echo "Microorganisms,reads,RPM,RA" > $path4/${OUT}_fRep/${OUT}_Virus.csv
  cat $path4/${OUT}_fRep/${OUT}_Virus1.csv >> $path4/${OUT}_fRep/${OUT}_Virus.csv
  rm -f $path4/${OUT}_fRep/${OUT}_Virus1.csv
  echo "Microorganisms,reads,RPM,RA" > $path4/${OUT}_fRep/${OUT}_Fungi.csv
  cat $path4/${OUT}_fRep/${OUT}_Fungi1.csv >> $path4/${OUT}_fRep/${OUT}_Fungi.csv
  rm -f $path4/${OUT}_fRep/${OUT}_Fungi1.csv

  echo "PATHDET_KRONA_FINISHED"

  touch $path4/.done
}
