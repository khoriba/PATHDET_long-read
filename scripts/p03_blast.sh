

split_fasta=${CONDA_PREFIX}/bin/split_fasta.py

#blast_threads=8
#parallel_jobs=`expr $THREADS \/ $blast_threads`
#parallel_jobs=8

## blastのスレッド数
blast_threads=20
## parallelのジョブ数
parallel_jobs=1
## 理論上blast_threads * parallel_jobsだけスレッドを使用
## 20スレッドの場合は効率的に並列化できないと思われるため、現状parallelは1並列設定

function p03_blast () {
  if [ ! -f ${path3}/.done ] ; then

    in1=${path2}/HGS/${OUT}_nohit.fasta

    #split FA
    mkdir -p ${path3}/${OUT}_FA
    mkdir -p ${path3}/${OUT}_BLAST
    python $split_fasta $in1 ${SPLIT} ${path3}/${OUT}_FA/${OUT}_

    # parallel
    if [ ! -f ${path3}/.${OUT}_blast.done ]; then
      ls ${path3}/${OUT}_FA/${OUT}_*.fasta \
        | parallel --no-notice --jobs $parallel_jobs \
          "[ ! -f ${path3}/${OUT}_BLAST/.{/.}_blast.done ] \
            && blastn \
              -db "${NT}" \
              -query {} \
              -out ${path3}/${OUT}_BLAST/{/.}_blast.txt \
              -evalue 1.0e-10 \
              -max_target_seqs 1 -max_hsps 5 \
              -outfmt \"7 std\" -num_threads $blast_threads \
            && touch ${path3}/${OUT}_BLAST/.{/.}_blast.done"
      if [ $? -eq 0 ]; then touch ${path3}/.${OUT}_blast.done; fi
    fi

    cat ${path3}/${OUT}_BLAST/*_blast.txt > ${path3}/${OUT}_blast_fmt7.txt

    cat ${path3}/${OUT}_blast_fmt7.txt \
      | awk '/hits found/{getline;print}' \
      | grep -v "#" > ${path3}/${OUT}_blast_hit.txt

    mv ${path3}/${OUT}_FA/${OUT}_summary.txt ${path3}/

    #OUTPUT Krona fmt (combination HTML excluding "no hit")
    ktClassifyBLAST -p -s ${path3}/${OUT}_blast_hit.txt -o ${path3}/${OUT}_BLAST/${OUT}_blast_hit.csv
    cut -f 1,2 ${path3}/${OUT}_BLAST/${OUT}_blast_hit.csv > ${path3}/${OUT}_BLAST/${OUT}_blast_hit2.csv

    ktImportTaxonomy \
      -m 1 -s 3 \
      ${path3}/${OUT}_BLAST/${OUT}_blast_hit2.csv \
      -o ${path3}/${OUT}_pathogen_final.html

    echo "PATHDET_BLAST_FINISHED"

    if [ -f ${path3}/.${OUT}_blast.done ] ; then
      touch $path3/.done
    fi
  else
    echo "Blast process already finished, skipped"
  fi
}

