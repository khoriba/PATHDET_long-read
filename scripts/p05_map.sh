#!/bin/bash

function fetch() {
  table=$1
  pattern=$2
  for i in `cat $table | grep -E "$pattern" | cut -f 3` ; do
    i=`echo $i | sed "s,db/ncbi_genome/,,g"`
    i="$path_ng/$i"
    fasta_ae+=($i)
    i=`basename $i`
    echo "$i" >> $log
  done
}

function get_fasta() {
  fasta_ae=()
  csv=$path4/${OUT}_CSV/${OUT}_merge_taxid.csv

  num=$1
  if [ "$num" == 1 ]; then ordinal=1st ; elif [ "$num" == 2 ]; then ordinal=2nd ; else ordinal=3rd ; fi
  THT=`cat $csv | sed -n ${num}p | cut -f 2`
  log=$path5/ref_list${num}.txt
  echo -e ">${ordinal}_Hit microbe reference genome" > $log
  if printf '%s\n' "${id_ref1[@]}" | grep -qx "$THT" ; then
    fetch $table_ref "\s$THT\s"
  else
    if printf '%s\n' "${id_ref2[@]}" | grep -qx "$THT" ; then
      fetch $table_ref "(^|,)$THT(,|\s)"
    else
      echo -e "Not applicable" >> $log
      echo -e ">${ordinal}_Hit microbe registered genome" >> $log
      if printf '%s\n' "${id_all1[@]}" | grep -qx "$THT" ; then
        fetch $table_all "\s$THT\s"
      else
        if printf '%s\n' "${id_all2[@]}" | grep -qx "$THT" ; then
          fetch $table_all "(^|,)$THT(,|\s)"
        else
          echo -e "Not applicable" >> $log
          echo -e "Automatic acquisition of $ordinal microbe complete genome failed." >> $log
        fi
      fi
    fi
  fi
  if [ "${#fasta_ae[@]}" -gt 0 ]; then
    if [ ! -d "$path5/MAP" ]; then mkdir -p $path5/MAP ; fi
    for i in ${fasta_ae[@]} ; do
      zcat $i >> $path5/MAP/REF.fasta
    done
  fi
}

function p05_map() {
  if [ -f "$path5/.done" ]; then
    echo "Map   process already finished, skipped"
    return
  fi

  mkdir -p $path5
  OUTL5=$path5/${OUT}_log_p5.txt

  table_ref=$path_ng/reference/refseq_lineage_taxid_fasta.tsv
  table_all=$path_ng/all/refseq_lineage_taxid_fasta.tsv
  id_ref1=(`cat $table_ref | sed 1d | cut -f 2`)
  id_ref2=(`cat $table_ref | sed 1d | cut -f 1 | sed "s/,/\n/g" | sort -u`)
  id_all1=(`cat $table_all | sed 1d | cut -f 2 | uniq`)
  id_all2=(`cat $table_all | sed 1d | cut -f 1 | sed "s/,/\n/g" | sort -u`)
  FA=$path5/MAP/REF.fasta
  if [ ! -f "$FA" ]; then
    for num in 1 2 3 ; do
      get_fasta $num
    done
  fi

  SRG=`date '+%y%m%d%H%M%S'`
  echo -e "${SRG} reference genome search finished." >> $OUTL5
  cat $path5/ref_list1.txt $path5/ref_list2.txt $path5/ref_list3.txt >> $OUTL5
  mv $path5/ref_list*.txt $path5/MAP/

  if [ -f "$FA" ] ; then
    samtools faidx $FA
  else
    echo "$FA does NOT exist"
    echo "Map process finished"
    touch $path5/.done
    return
  fi

  # PATH setting
  mkdir -p $path5/${OUT}_map $path5/${OUT}_con $path5/${OUT}_bw

  # Mapping
  bam=$path5/${OUT}_map/${OUT}.sort.bam
  if [ "$MODE" == "short" ]; then
    # input: bbmap output
    in1=$path1/${OUT}_trimmed2_R1.fastq.gz
    in2=$path1/${OUT}_trimmed2_R2.fastq.gz
    minimap2 -t $THREADS -ax sr $FA $in1 $in2 \
      | samtools sort -o $bam
    samtools index $bam
  elif [ "$MODE" == "long" ]; then
    in1=${path1}/${OUT}_trimmed_L.fastq
    minimap2 -t $THREADS -ax map-ont $FA $in1 \
      | samtools sort -o $bam
    samtools index $bam
  fi
  pileup.sh in=$bam ref=$FA out=$path5/${OUT}_map/${OUT}_map.csv 

  # Variant call
  bcftools mpileup -Ou -f $FA $bam \
    | bcftools call -mv -Oz -o $path5/${OUT}_con/${OUT}.vcf.gz 
  bcftools index $path5/${OUT}_con/${OUT}.vcf.gz

  bcftools norm -f ${FA} $path5/${OUT}_con/${OUT}.vcf.gz \
    -Ob -o $path5/${OUT}_con/${OUT}.norm.bcf 
  bcftools filter --IndelGap 5 $path5/${OUT}_con/${OUT}.norm.bcf \
    -Oz -o $path5/${OUT}_con/${OUT}.norm.flt-idls.vcf.gz 
  bcftools index $path5/${OUT}_con/${OUT}.norm.flt-idls.vcf.gz

  cat ${FA} \
    | bcftools consensus $path5/${OUT}_con/${OUT}.norm.flt-idls.vcf.gz \
    > $path5/${OUT}_con/${OUT}_consensus.fa 

  #BigWig making
  bamCoverage -b $bam -o $path5/${OUT}_bw/${OUT}.bw -of bigwig \
    --binSize=10 --minMappingQuality 10 --numberOfProcessors=max 
  bdg=$path5/${OUT}_bw/${OUT}.bdg
  bamCoverage -b $bam -o $bdg -bs 1 -of bedgraph 
  bgzip $bdg
  tabix -p bed $bdg.gz

  #Drawing Mappung Graph (SVG)
  cp $FA $path5/${OUT}_bw/${OUT}_ref.fa
  while read line ; do
    contig=`echo $line | cut -d " " -f 1`
    len=`echo $line | cut -d " " -f 2`
    FAH="$contig:1-$len"
    python $spark -pr $FAH \
      -cf $bdg.gz -o $path5/${OUT}_bw/spark_$contig 
  done < $FA.fai

  echo "PATHDET_MAP_FINISHED"

  touch $path5/.done
}
