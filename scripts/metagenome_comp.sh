#!/bin/bash
## METAGENOME_COMP
## FOR COMPARISON OF MICROBIAL FLORA
## Please specify the taxonomy to compare

START=`date '+%y%m%d%H%M%S'`
TIMEA=`date +%s`

##Argument and Variable
##TAX='K','F','G','S','Fungi','Virus'
TAX="$1"
OUT="$2"
OUTI=${OUT}_int
OUTL=${OUT}_log.txt

##PATH
PATH0=output
PATH1=output_comp/${START}_${OUT}
PATH2=$PATH1/${OUTI}

##Output directory
mkdir -p $PATH2

##get metagenome_files
echo "${START} ${OUT} processing start" >> ${PATH1}/${OUTL}
echo "${START} column extraction" >> ${PATH1}/${OUTL}
mkdir -p ${PATH2}/COL1
mkdir -p ${PATH2}/COL2
mkdir -p ${PATH2}/COL3

exec 2> ${PATH2}/er.log

dirs=`find $PATH0/*/*_tbl/*_fRep/ -size +1c -type f -name *${TAX}.csv | sort`

for dir in ${dirs} ;
do
    fname_ext="${dir##*/}"
    fname="${fname_ext%.*}"
    cut -f 1 ${dir} > ${PATH2}/COL1/${fname}_c1.csv
    cut -f 1,2 ${dir} | sort -k 1,1 > ${PATH2}/COL2/${fname}_rpm.csv
    cut -f 1,3 ${dir} | sort -k 1,1 > ${PATH2}/COL3/${fname}_riv.csv
done

##merge metagenome filles column1
MER1=`date '+%y%m%d%H%M%S'`
echo "${MER1} RPM comparison started" >> ${PATH1}/${OUTL}
mkdir -p ${PATH2}/MERGE
cat ${PATH2}/COL1/*_c1.csv | sort -u > ${PATH1}/pathogen_list.csv
rpms=`find ${PATH2}/COL2/ -type f -name *.csv | sort`
mcsv=${PATH1}/pathogen_list.csv
echo "pathogen_name" > ${PATH2}/head.csv
headt=${PATH2}/head.csv

for rpm in ${rpms} ;
do
    fname_ext="${rpm##*/}"
    join -1 1 -2 1 -a 1 -t "$(printf '\011')" -o auto -e "0" ${mcsv} ${rpm} > ${PATH2}/MERGE/merge.csv
    mv ${PATH2}/MERGE/merge.csv ${PATH2}/MERGE/merge1.csv
    mcsv=${PATH2}/MERGE/merge1.csv
    echo "${fname_ext%.*}" > ${PATH2}/MERGE/head1.csv
    paste ${headt} ${PATH2}/MERGE/head1.csv > ${PATH2}/MERGE/head.csv
    mv ${PATH2}/MERGE/head.csv ${PATH2}/MERGE/headt.csv
    headt=${PATH2}/MERGE/headt.csv
done
cat ${headt} ${mcsv} > ${PATH1}/${OUT}_${TAX}_mRPM.csv
cat ${PATH1}/${OUT}_${TAX}_mRPM.csv | perl -anF'\t|\n' -e'$n=@F-1if!$n;for(0..$n){push@{$$m[$_]},$F[$_]}''END{print map{join"\t",@$_,"\n"}@$m}' > ${PATH1}/${OUT}_${TAX}_rRPM.csv

##merge metagenome filles column2
MER2=`date '+%y%m%d%H%M%S'`
echo "${MER2} RIV comparison started" >> ${PATH1}/${OUTL}
rivs=`find ${PATH2}/COL3/ -type f -name *.csv | sort`
mcsv=${PATH1}/pathogen_list.csv
headt=${PATH2}/head.csv

for riv in ${rivs} ;
do
    fname_ext="${riv##*/}"
    join -1 1 -2 1 -a 1 -t "$(printf '\011')" -o auto -e "0" ${mcsv} ${riv} > ${PATH2}/MERGE/merge.csv
    mv ${PATH2}/MERGE/merge.csv ${PATH2}/MERGE/merge2.csv
    mcsv=${PATH2}/MERGE/merge2.csv
    echo "${fname_ext%.*}" > ${PATH2}/MERGE/head2.csv
    paste ${headt} ${PATH2}/MERGE/head2.csv > ${PATH2}/MERGE/head.csv
    mv ${PATH2}/MERGE/head.csv ${PATH2}/MERGE/headt.csv
    headt=${PATH2}/MERGE/headt.csv
done
cat ${headt} ${mcsv} > ${PATH1}/${OUT}_${TAX}_mRIV.csv
cat ${PATH1}/${OUT}_${TAX}_mRIV.csv | perl -anF'\t|\n' -e'$n=@F-1if!$n;for(0..$n){push@{$$m[$_]},$F[$_]}''END{print map{join"\t",@$_,"\n"}@$m}' > ${PATH1}/${OUT}_${TAX}_rRIV.csv

rm -rf ${PATH2}

TIMEB=`date +%s`
PT=`expr ${TIMEB} \- ${TIMEA}`
H=`expr ${PT} \/ 3600`
PT=`expr ${PT} \% 3600`
M=`expr ${PT} \/ 60`
S=`expr ${PT} \% 60`
echo "total processing time ${H}:${M}:${S}" >> ${PATH1}/${OUTL}
echo "metagenome_comp_v1.sh finished." >> ${PATH1}/${OUTL}
