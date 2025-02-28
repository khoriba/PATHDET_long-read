#!/bin/bash

set -u

cd Data
THREADS=64
SPLIT=50000
LEN=1000
DEBUG=False
MODE=long

while [ $# -gt 0 ] ; do
  case "$1" in
    -i             ) FQ1=$2
                     shift
                     shift ;;
    -o | --out     ) OUT=$2
                     shift
                     shift ;;
    -t | --threads ) THREADS=$2
                     shift
                     shift ;;
    -s | --split   ) SPLIT=$2
                     shift
                     shift ;;
    -l | --length   ) LEN=$2
                     shift
                     shift ;;
    --debug        ) DEBUG=True
                     shift ;;
                 * ) echo "ERROR : illegal option"
                     echo "you need to supply path to fastq using -1 and -2 option, and output name using -o/--out option"
                     echo "you could overwrite thread number using -t/--threads option, 64 by default"
                     echo "you could overwrite split number using -s/--split option, 50000 by default"
                     exit ;;
  esac
done

if [ -z "$FQ1" ] || [ -z "$OUT" ]; then
  echo "you need to supply path to fastq using -1 and -2 option, and output name using -o option"
  exit
elif [ ! -f "$FQ1" ]; then
  echo "fastq file does NOT exists"
  exit
fi

##PATH
path0=output/${OUT}
mkdir -p $path0

path1=$path0/${OUT}_qc
path2=$path0/${OUT}_hgs
pathr=$path0/${OUT}_rapid
path3=$path0/${OUT}_blast
path4=$path0/${OUT}_tbl
path5=$path0/${OUT}_map

##DataBase
datadir="/path/to/datadir"
GRCH38="${datadir}/kraken2/human"
HUMAN="${datadir}/hg38.fa"
NT=${datadir}/blast/nt/nt
HOST2=${datadir}/blast/t2t/GCF_009914755.1_T2T-CHM13v2.0_genomic
PATHDB1=${datadir}/kraken2/k2_pluspf_20241228
path_ng=${datadir}/ncbi_genome
MASHDB=${datadir}/mash/RefSeq88n.msh

spark=$CONDA_PREFIX/opt/git/SparK/SparK.py
export TAXONKIT_DB=${datadir}/taxonkit/taxdump

function p01_qc () {
  ## ===========================
  ## pathdet_p1L.sh
  ## ===========================
  
  ##LOG
  local START=`date '+%y%m%d%H%M%S'`
  local TIMEA=`date +%s`
  local OUTL=${OUT}_log_p1.txt
  if [ ! -f ${path1}/.done ] ; then

    mkdir -p $path1
  
    echo "${START} ${OUT} raw FQ stats" >> ${path1}/${OUTL}
    seqkit stats -j $THREADS ${FQ1} >> ${path1}/${OUTL}
  
    ##QC_1st
    NanoPlot --fastq ${FQ1} --loglength -t $THREADS -o ${path1}/${OUT}_qc1
  
    ##Trimming, default_length=1,000
    gzip -dc ${FQ1} > ${path1}/${OUT}_L.fastq
    NanoFilt ${path1}/${OUT}_L.fastq -q 10 --headcrop 50 -l ${LEN} > ${path1}/${OUT}_trimmed_L.fastq
  
    ##QC_2nd
    NanoPlot --fastq ${path1}/${OUT}_trimmed_L.fastq --loglength -t $THREADS -o ${path1}/${OUT}_qc2
  
    QC=`date '+%y%m%d%H%M%S'`
    echo "${QC} ${OUT} trimmed FQ stats" >> ${path1}/${OUTL}
    seqkit stats -j $THREADS ${path1}/${OUT}_trimmed_L.fastq >> ${path1}/${OUTL}

    touch $path1/.done
  fi
}

function p02_hgs () {
  ## ===========================
  ## pathdet_p2L.sh
  ## ===========================

  local FQ1=${path1}/${OUT}_trimmed_L.fastq

  ##LOG
  local START=`date '+%y%m%d%H%M%S'`
  local TIMEA=`date +%s`
  local OUTL=${OUT}_log_p2.txt

  if [ ! -f ${path2}/.done ] ; then
    ##Human Genome Subtraction
    mkdir -p ${path2}/HGS
    if [ ! -f ${path2}/.${OUT}_hgs1.done ]; then
      kraken2 \
        --threads $THREADS \
        --db ${GRCH38} \
        ${FQ1} \
        --unclassified-out ${path2}/HGS/${OUT}_unclass1.fastq \
        > ${path2}/HGS/${OUT}.kraken.txt
      if [ $? -eq 0 ]; then touch ${path2}/.${OUT}_hgs1.done; fi
    fi
    if [ ! -f ${path2}/.${OUT}_hgs2.done ]; then
      if  [ -f ${path2}/.${OUT}_hgs1.done ]; then
        minimap2 -t $THREADS -ax map-ont ${HUMAN} ${path2}/HGS/${OUT}_unclass1.fastq \
        > ${path2}/HGS/${OUT}_humanj.sam
        samtools view -@ $THREADS -bS ${path2}/HGS/${OUT}_humanj.sam | samtools sort | samtools view -f 4 > ${path2}/HGS/${OUT}_unmap.sam
        samtools fastq -@ $THREADS ${path2}/HGS/${OUT}_unmap.sam > ${path2}/HGS/${OUT}_hgs.fastq
        if [ $? -eq 0 ]; then touch ${path2}/.${OUT}_hgs2.done; fi
      else
        echo An unexpected error.
        exit 1
      fi
    fi
    
    #Host Genome Subtraction plus
    if [ ! -f ${path2}/.${OUT}_hgs_plus.done ]; then
      if  [ -f ${path2}/.${OUT}_hgs1.done ] && [ -f ${path2}/.${OUT}_hgs2.done ]; then
        #FQ2FA
        seqkit fq2fa ${path2}/HGS/${OUT}_hgs.fastq -o ${path2}/HGS/${OUT}_hgs.fasta
        blastn \
          -db ${HOST2} \
          -query ${path2}/HGS/${OUT}_hgs.fasta \
          -evalue 1.0e-10 \
          -max_target_seqs 1 -max_hsps 5 \
          -outfmt "7 std" -num_threads $THREADS > ${path2}/HGS/${OUT}_hgs.txt
        cat ${path2}/HGS/${OUT}_hgs.txt \
          | grep -B 2 "# 0 hits found" \
          | grep "Query" \
          | sed -e 's/# Query: //g' \
          > ${path2}/HGS/${OUT}_nohitID.txt
  
        seqkit grep -w0 -f ${path2}/HGS/${OUT}_nohitID.txt ${path2}/HGS/${OUT}_hgs.fasta > ${path2}/HGS/${OUT}_nohit.fasta
  
        touch ${path2}/.${OUT}_hgs_plus.done 
      else
        echo An unexpected error.
        exit 1
      fi
    fi

    if [ -f ${path2}/.${OUT}_hgs_plus.done ]; then
      local QC=`date '+%y%m%d%H%M%S'`
      echo "${QC} ${OUT} host genome subtracted FQ stats" >> ${path2}/${OUTL}
      seqkit stats ${path2}/HGS/${OUT}_nohit.fasta >> ${path2}/${OUTL}
      touch $path2/.done
    fi
  fi
}

function rapid () {
  ## ===================
  ## Rapid (pathdet3_v6.0.sh)
  ## ===================
  
  ## Report
  local OUTR=${OUT}_report.csv
  local OUTP=${OUT}_pre-report.csv

  ##LOG
  local START=`date '+%y%m%d%H%M%S'`
  local TIMEA=`date +%s`
  local OUTL=${OUT}_log_rapid.txt
  if [ ! -f ${pathr}/.done ] ; then
  
    ##PATHOGEN detection :MASH screen
    
    mkdir -p ${pathr}/MS
    
    if [ ! -f ${pathr}/.${OUT}_mash_screen.done ]; then
      mash screen \
        -p $THREADS -w \
        $MASHDB \
        ${path2}/HGS/${OUT}_nohit.fasta \
        > ${pathr}/MS/${OUT}_screen.tab
      if [ $? -eq 0 ]; then touch ${pathr}/.${OUT}_mash_screen.done; fi
      
      echo -e "\n[MASH_SCREEN_REPORT]" >> ${pathr}/${OUTR}
      echo -e "identity, shared-hashes, median-multiplicity, p-value, query-ID, query-comment" >> ${pathr}/${OUTR}
      sort -gr ${pathr}/MS/${OUT}_screen.tab | head >> ${pathr}/${OUTR}
      
      local MASH=`date '+%y%m%d%H%M%S'`
      echo "${MASH} pathogen screening(RefSeq release 88) finished" >> ${pathr}/${OUTL}
    
    fi

    local HG2=`seqkit stat -T ${path2}/HGS/${OUT}_nohit.fasta | tail -n1 | cut -f4`

    local WCA1=`grep 'Number of reads' ${path1}/${OUT}_qc1/NanoStats.txt | sed -e 's/^Number of reads:\s\+\([0-9,]\+\)\.0$/\1/' | sed -e 's/,//g'`
    
    ##PATHOGEN detection :Kraken2
    mkdir -p ${pathr}/KK2
    
    if [ ! -f ${pathr}/.${OUT}_kraken_p1.done ]; then
      kraken2 \
        --threads $THREADS \
        --report ${pathr}/KK2/${OUT}_p1.kreport \
        --unclassified-out ${pathr}/KK2/${OUT}_unclass2.fasta \
        --classified-out ${pathr}/KK2/${OUT}_class2.fasta \
        --db ${PATHDB1} ${path2}/HGS/${OUT}_nohit.fasta \
          > ${pathr}/KK2/${OUT}_kraken_p1.txt
      if [ $? -eq 0 ]; then touch ${pathr}/.${OUT}_kraken_p1.done; fi
        #--confidence 0.5 > ${pathr}/KK2/${OUT}_kraken_p1.txt
      ## 一旦confidence指定なしで進める方針
    
      local KK2S=`date '+%y%m%d%H%M%S'`
      local PKK2=`seqkit stat -T ${pathr}/KK2/${OUT}_class2.fasta | tail -n1 | cut -f4`
      local UKK2=`seqkit stat -T ${pathr}/KK2/${OUT}_unclass2.fasta | tail -n1 | cut -f4`
      local KKAR=`echo "scale=4; ${PKK2} / ${HG2} * 100" | bc`
      local UCR=`echo "scale=4; ${UKK2} / ${WCA1} * 100" | bc`
      echo "${KK2S} kraken2 finished: ${PKK2}reads (${KKAR}% of karaken2 target)" >> ${pathr}/${OUTL}
      echo -e "\t* un-classified read: ${UKK2}reads (${UCR}%)" >> ${pathr}/${OUTL}
    fi
  
    if [ ! -f ${pathr}/.${OUT}_kraken_prompt_report.done ]; then
      if [ -f ${pathr}/.${OUT}_kraken_p1.done ]; then
    
        local PRstd=`echo "scale=0; ${HG2} * 1000000 / ${WCA1}" | bc`
    
        local RMDA2=`grep 'Number of reads' ${path1}/${OUT}_qc2/NanoStats.txt | sed -e 's/^Number of reads:\s\+\([0-9,]\+\)\.0$/\1/' | sed -e 's/,//g'`
    
        ##Prompt report
        mkdir -p ${pathr}/${OUT}_pRep
        echo -e "Sample_ID\t${OUT}" >> ${pathr}/${OUT}_pRep/${OUTP}
        echo -e "Total_read(reads)\t${WCA1}" >> ${pathr}/${OUT}_pRep/${OUTP}
        echo -e "Available_read(reads)\t${RMDA2}" >> ${pathr}/${OUT}_pRep/${OUTP}
        echo -e "Non-human_read(reads)\t${HG2}" >> ${pathr}/${OUT}_pRep/${OUTP}
        echo -e "Non-human_read(RPM)\t${PRstd}" >> ${pathr}/${OUT}_pRep/${OUTP}
        
        ktImportTaxonomy \
          -q 2 -t 3 -s 4 \
          ${pathr}/KK2/${OUT}_kraken_p1.txt \
          -o ${pathr}/${OUT}_pRep/${OUT}.plot.html
          grep -v 'unclassified' ${pathr}/KK2/${OUT}_p1.kreport | cut -f3,5 | /bin/awk -F, '$1 > 1{ print $0 }' >  ${pathr}/KK2/${OUT}_p2.kreport
    
        cat ${pathr}/KK2/${OUT}_p2.kreport | \
          /bin/awk '{ tax[$2] += $1; } END { for (i in tax) print tax[i] "\t" i; }' | sort -k 1nr,1 | \
        taxonkit lineage -i 2 \
          > ${pathr}/KK2/${OUT}_tax.kreport
        cat ${pathr}/KK2/${OUT}_tax.kreport | /bin/awk -v "r=${RMDA2}" '{print $1/r*1000000"\t"$0}' | cut -f 1,3- \
          > ${pathr}/KK2/${OUT}_tax2.kreport
        
        #Remove C.acnes
        #grep -v acnes ${pathr}/KK2/${OUT}_tax2.kreport > ${pathr}/KK2/${OUT}_tax3.kreport
        #grep -v acnes ${pathr}/KK2/${OUT}_tax.kreport > ${pathr}/KK2/${OUT}_tax3r.kreport
    
        #Create CSV file for each taxonomic rank
        mkdir -p ${pathr}/CSV
        
        #Kingdom
        taxonkit reformat -j 8 -i 3 -f {k} ${pathr}/KK2/${OUT}_tax2.kreport\
         | cut -f 1,4- | sed 's/ /_/g' | /bin/awk '{ tax[$2] += $1; } END { for (i in tax) print i "\t" tax[i]; }'\
         | sort -k 2nr,2\
         | (s=$(cat)&&/bin/awk 'NR==FNR{a=a+$2;next}{p=$2/a;print $0"\t"p}' <(echo -e "$s") <(echo -e "$s"))\
         > ${pathr}/CSV/${OUT}_mpK.csv
        
        taxonkit reformat -j 8 -i 3 -f {k} ${pathr}/KK2/${OUT}_tax.kreport\
         | cut -f 1,4- | sed 's/ /_/g' | /bin/awk '{ tax[$2] += $1; } END { for (i in tax) print i "\t" tax[i]; }'\
         | sort -k 2nr,2\
         > ${pathr}/CSV/${OUT}_rpK.csv
        
        paste ${pathr}/CSV/${OUT}_rpK.csv ${pathr}/CSV/${OUT}_mpK.csv | cut -f 1,2,4,5\
         > ${pathr}/${OUT}_pRep/${OUT}_pK1.csv
        
        #Phylum
        taxonkit reformat -j 8 -i 3 -f {p} ${pathr}/KK2/${OUT}_tax2.kreport\
         | cut -f 1,4- | sed 's/ /_/g' | /bin/awk '{ tax[$2] += $1; } END { for (i in tax) print i "\t" tax[i]; }'\
         | sort -k 2nr,2 | grep -v -e '^\s'\
         | (s=$(cat)&&/bin/awk 'NR==FNR{a=a+$2;next}{p=$2/a;print $0"\t"p}' <(echo -e "$s") <(echo -e "$s"))\
         > ${pathr}/CSV/${OUT}_mpP.csv
        
        taxonkit reformat -j 8 -i 3 -f {p} ${pathr}/KK2/${OUT}_tax.kreport\
         | cut -f 1,4- | sed 's/ /_/g' | /bin/awk '{ tax[$2] += $1; } END { for (i in tax) print i "\t" tax[i]; }'\
         | sort -k 2nr,2 | grep -v -e '^\s'\
         > ${pathr}/CSV/${OUT}_rpP.csv
        
        paste ${pathr}/CSV/${OUT}_rpP.csv ${pathr}/CSV/${OUT}_mpP.csv | cut -f 1,2,4,5\
         > ${pathr}/${OUT}_pRep/${OUT}_pP1.csv
        
        #Family
        taxonkit reformat -j 8 -i 3 -f {f} ${pathr}/KK2/${OUT}_tax2.kreport\
         | cut -f 1,4- | sed 's/ /_/g' | /bin/awk '{ tax[$2] += $1; } END { for (i in tax) print i "\t" tax[i]; }'\
         | sort -k 2nr,2 | grep -v -e '^\s'\
         | (s=$(cat)&&/bin/awk 'NR==FNR{a=a+$2;next}{p=$2/a;print $0"\t"p}' <(echo -e "$s") <(echo -e "$s"))\
         > ${pathr}/CSV/${OUT}_mpF.csv
        
        taxonkit reformat -j 8 -i 3 -f {f} ${pathr}/KK2/${OUT}_tax.kreport\
         | cut -f 1,4- | sed 's/ /_/g' | /bin/awk '{ tax[$2] += $1; } END { for (i in tax) print i "\t" tax[i]; }'\
         | sort -k 2nr,2 | grep -v -e '^\s'\
         > ${pathr}/CSV/${OUT}_rpF.csv
        
        paste ${pathr}/CSV/${OUT}_rpF.csv ${pathr}/CSV/${OUT}_mpF.csv | cut -f 1,2,4,5\
         > ${pathr}/${OUT}_pRep/${OUT}_pF1.csv
    
        #Genus
        taxonkit reformat -j 8 -i 3 -f {g} ${pathr}/KK2/${OUT}_tax2.kreport\
         | cut -f 1,4- | sed 's/ /_/g' | /bin/awk '{ tax[$2] += $1; } END { for (i in tax) print i "\t" tax[i]; }'\
         | sort -k 2nr,2 | grep -v -e '^\s'\
         | (s=$(cat)&&/bin/awk 'NR==FNR{a=a+$2;next}{p=$2/a;print $0"\t"p}' <(echo -e "$s") <(echo -e "$s"))\
         > ${pathr}/CSV/${OUT}_mpG.csv
        
        taxonkit reformat -j 8 -i 3 -f {g} ${pathr}/KK2/${OUT}_tax.kreport\
         | cut -f 1,4- | sed 's/ /_/g' | /bin/awk '{ tax[$2] += $1; } END { for (i in tax) print i "\t" tax[i]; }'\
         | sort -k 2nr,2 | grep -v -e '^\s'\
         > ${pathr}/CSV/${OUT}_rpG.csv
        
        paste ${pathr}/CSV/${OUT}_rpG.csv ${pathr}/CSV/${OUT}_mpG.csv | cut -f 1,2,4,5\
         > ${pathr}/${OUT}_pRep/${OUT}_pG1.csv
        
        #Species
        taxonkit reformat -j 8 -i 3 -f {s} ${pathr}/KK2/${OUT}_tax2.kreport\
         | cut -f 1,4- | sed 's/ /_/g' | /bin/awk '{ tax[$2] += $1; } END { for (i in tax) print i "\t" tax[i]; }'\
         | sort -k 2nr,2 | grep -v -e '^\s'\
         | (s=$(cat)&&/bin/awk 'NR==FNR{a=a+$2;next}{p=$2/a;print $0"\t"p}' <(echo -e "$s") <(echo -e "$s"))\
         > ${pathr}/CSV/${OUT}_mpS.csv
        
        taxonkit reformat -j 8 -i 3 -f {s} ${pathr}/KK2/${OUT}_tax.kreport\
         | cut -f 1,4- | sed 's/ /_/g' | /bin/awk '{ tax[$2] += $1; } END { for (i in tax) print i "\t" tax[i]; }'\
         | sort -k 2nr,2 | grep -v -e '^\s'\
         > ${pathr}/CSV/${OUT}_rpS.csv
        
        paste ${pathr}/CSV/${OUT}_rpS.csv ${pathr}/CSV/${OUT}_mpS.csv | cut -f 1,2,4,5\
         > ${pathr}/${OUT}_pRep/${OUT}_pS1.csv
        
        #Calculate parameters(RIV,DiversityIndex D,H)
        PF1=`head -n3 ${pathr}/${OUT}_pRep/${OUT}_pF1.csv`
        PG1=`head -n3 ${pathr}/${OUT}_pRep/${OUT}_pG1.csv`
        PS1=`head -n3 ${pathr}/${OUT}_pRep/${OUT}_pS1.csv`
        
        PDF=`cut -f3 ${pathr}/CSV/${OUT}_mpF.csv | perl -nle '$D += $_*$_; print 1-$D if(eof);'`
        PDG=`cut -f3 ${pathr}/CSV/${OUT}_mpG.csv | perl -nle '$D += $_*$_; print 1-$D if(eof);'`
        PDS=`cut -f3 ${pathr}/CSV/${OUT}_mpS.csv | perl -nle '$D += $_*$_; print 1-$D if(eof);'`
        
        PHF=`cut -f3 ${pathr}/CSV/${OUT}_mpF.csv | perl -nle '$H += $_*log($_); print $H*(-1) if(eof);'`
        PHG=`cut -f3 ${pathr}/CSV/${OUT}_mpG.csv | perl -nle '$H += $_*log($_); print $H*(-1) if(eof);'`
        PHS=`cut -f3 ${pathr}/CSV/${OUT}_mpS.csv | perl -nle '$H += $_*log($_); print $H*(-1) if(eof);'`
        
        echo -e "\n[PATHDET_REPORT]" >> ${pathr}/${OUT}_pRep/${OUTP}
        echo -e ">Diversity_Index" >> ${pathr}/${OUT}_pRep/${OUTP}
        echo -e "Simpson_Index(family)\t${PDF}" >> ${pathr}/${OUT}_pRep/${OUTP}
        echo -e "Simpson_Index(genus)\t${PDG}" >> ${pathr}/${OUT}_pRep/${OUTP}
        echo -e "Simpson_Index(species)\t${PDS}" >> ${pathr}/${OUT}_pRep/${OUTP}
        echo -e "\nShannon_Index(family)\t${PHF}" >> ${pathr}/${OUT}_pRep/${OUTP}
        echo -e "Shannon_Index(genus)\t${PHG}" >> ${pathr}/${OUT}_pRep/${OUTP}
        echo -e "Shannon_Index(species)\t${PHS}" >> ${pathr}/${OUT}_pRep/${OUTP}
        echo -e "\n>Top_3hit_pathogen(family) Name,reads,RPM,RA\n${PF1}" >> ${pathr}/${OUT}_pRep/${OUTP}
        echo -e "\n>Top_3hit_pathogen(genus) Name,reads,RPM,RA\n${PG1}" >> ${pathr}/${OUT}_pRep/${OUTP}
        echo -e "\n>Top_3hit_pathogen(species) Name,reads,RPM,RA\n${PS1}" >> ${pathr}/${OUT}_pRep/${OUTP}
        
        echo -e "Microorganisms,reads,RPM,RA" >  ${pathr}/${OUT}_pRep/${OUT}_pK.csv
        cat ${pathr}/${OUT}_pRep/${OUT}_pK1.csv  >> ${pathr}/${OUT}_pRep/${OUT}_pK.csv
        echo -e "Microorganisms,reads,RPM,RA" >  ${pathr}/${OUT}_pRep/${OUT}_pP.csv
        cat ${pathr}/${OUT}_pRep/${OUT}_pP1.csv  >> ${pathr}/${OUT}_pRep/${OUT}_pP.csv
        echo -e "Microorganisms,reads,RPM,RA" >  ${pathr}/${OUT}_pRep/${OUT}_pF.csv
        cat ${pathr}/${OUT}_pRep/${OUT}_pF1.csv  >> ${pathr}/${OUT}_pRep/${OUT}_pF.csv
        echo -e "Microorganisms,reads,RPM,RA" >  ${pathr}/${OUT}_pRep/${OUT}_pG.csv
        cat ${pathr}/${OUT}_pRep/${OUT}_pG1.csv  >> ${pathr}/${OUT}_pRep/${OUT}_pG.csv
        echo -e "Microorganisms,reads,RPM,RA" >  ${pathr}/${OUT}_pRep/${OUT}_pS.csv
        cat ${pathr}/${OUT}_pRep/${OUT}_pS1.csv  >> ${pathr}/${OUT}_pRep/${OUT}_pS.csv
        rm  ${pathr}/${OUT}_pRep/*1.csv
        
        TIMEC=`date +%s`
        PT1=`expr ${TIMEC} \- ${TIMEA}`
        H1=`expr ${PT1} \/ 3600`
        PT1=`expr ${PT1} \% 3600`
        M1=`expr ${PT1} \/ 60`
        S1=`expr ${PT1} \% 60`
        BRKS=`date '+%y%m%d%H%M%S'`
        echo "${BRKS} Kraken REPORT finished (${H1}:${M1}:${S1})" >> ${pathr}/${OUTL}
    
        touch ${pathr}/.${OUT}_kraken_prompt_report.done
      fi
    else
      echo An unexpected error.
      exit 1
    fi
    if [ -f ${pathr}/.${OUT}_kraken_prompt_report.done ]; then
      touch $pathr/.done
    fi
  fi
}

function tidyup() {
  intermediates=(
    $path1/${OUT}_L.fastq
    $path1/${OUT}_trimmed_L.fastq
    $path2/HGS
    $path3/${OUT}_FA
    $path3/${OUT}_BLAST
    $path3/${OUT}_hgs.fasta
    $path4/${OUT}_CSV
    $path4/${OUT}_CSV2
  )
  if [ -f "$path5/.done" ] && [ "$DEBUG" == "False" ]; then
    rm -fr ${intermediates[@]}
  fi
}

. ${CONDA_PREFIX}/bin/p03_blast.sh
. ${CONDA_PREFIX}/bin/p04_krona.sh
. ${CONDA_PREFIX}/bin/p05_map.sh

functions=(p01_qc p02_hgs rapid p03_blast p04_krona p05_map)
for func in ${functions[@]} ; do $func ; done

tidyup
