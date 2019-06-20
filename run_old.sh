ANADIR="/dati/Data_CGEM_IHEP_Integration_2019/raw_root/$1"
HERE=$PWD

if [ ! -d ${ANADIR} ]
then
    mkdir ${ANADIR}
else
    rm -f ${ANADIR}/Spill*.root
    rm -f ${ANADIR}/event.root
    rm -f ${ANADIR}/ana.root
    rm -f ${ANADIR}/log
fi

#if [ ! -f "${DATADIR}/file_list" ]; then
#    echo "cannot find \"file_list\" in ${DATADIR}"
#    return;
#fi

rm Spill_*GEMROC*.* -f
r=1
#ls ${DATADIR}
if [[ $1 -lt 88 ]]; then
    DATADIR="/dati/Data_CGEM_IHEP_Integration_2019/raw_dat/$1"
    for ROC in `ls ${DATADIR}/Spill* | sed 's/^.\{,65\}//' | sed 's/.$//' | sed 's/.$//'| sed 's/.$//' | sed 's/.$//'   `
    do
	#echo ROC: ${ROC}
	#echo FILE: ${DATADIR}/Spill_$1_GEMROC_*.dat
	cp ${DATADIR}/Spill_$1_GEMROC_*.dat .
	
	#echo NOME: Spill_$1_GEMROC_${ROC}.dat
	if [ -f "Spill_$1_GEMROC_${ROC}.dat" ]; then
	    python Decode.py Spill_$1_GEMROC_${ROC}.dat ${ROC} 1 $r # triggermatch
	    #python Decode.py Spill_$1_GEMROC_${ROC}.dat ${ROC} 0 $r #triggerless
	fi
	
	#   ls
    done  
    hadd -f decode.root Spill_$1_GEMROC*root
    ./bin/ana
    ./bin/event
    ./bin/post_event
 
    mv decode.root ${ANADIR}/Spill_$1_decode.root
    mv ana.root ${ANADIR}/Spill_$1_ana.root
    mv event.root ${ANADIR}/Spill_$1_event.root
    mv channel_chip.pdf ${ANADIR}/.
    echo "RunNo. $r\t${ROC}*" >> ${ANADIR}/log

    rm Spill_*GEMROC*.* -f

fi
if [[ $1 -gt 88 ]]; then
    DATADIR="/dati/Data_CGEM_IHEP_Integration_2019/raw_dat/RUN_$1"
    #ls ${DATADIR}/SubRUN*
    #for ROC in `ls ${DATADIR}/SubRUN* | sed 's/^.\{,69\}//' | sed 's/.$//' | sed 's/.$//'| sed 's/.$//' | sed 's/.$//'  | sed 's/.$//' | sed 's/.$//' | sed 's/.$//' `
    for i in {0..1000}
    do
	for ROC in {0..20}
	do
            echo ROC: ${ROC}  
	    echo sub: ${i}
            #echo FILE: ${DATADIR}/SubRUN_${i}_GEMROC_*.dat  
            #echo NOME: SubRUN_${i}_GEMROC_${ROC}_TM.dat 
            if [ -f "${DATADIR}/SubRUN_${i}_GEMROC_${ROC}_TM.dat" ]; then
		#echo ok
		#echo decode SubRUN_${i}_GEMROC_${ROC}_TM.dat
		#cp ${DATADIR}/SubRUN_${i}_GEMROC_${ROC}_TM.dat .
		python Decode.py ${DATADIR}/SubRUN_${i}_GEMROC_${ROC}_TM.dat ${ROC} 1 $r # triggermatch 
                #python Decode.py ${DATADIR}/SubRUN_${i}_GEMROC_${ROC}_TL.dat ${ROC} 1 $r # triggerless     
	    fi 
	done
	echo "SubRUN_${i}_GEMROC_*_TM.root"
	if [ -f "${DATADIR}/SubRUN_${i}_GEMROC_0_TM.root" ]; then
	    echo "Hello"
	    hadd -f ${DATADIR}/decode.root ${DATADIR}/SubRUN_*_GEMROC*root
	    echo "Begin ana"
            ./bin/ana ${DATADIR}
	    echo "Begin evt"
            ./bin/event ${DATADIR}
	    ./bin/post_event ${DATADIR}
            mv ${DATADIR}/decode.root      ${ANADIR}/Sub_RUN_dec_${i}.root
            mv ${DATADIR}/ana.root         ${ANADIR}/Sub_RUN_ana_${i}.root
            mv ${DATADIR}/event.root       ${ANADIR}/Sub_RUN_event_${i}.root
	    mv ${DATADIR}/post_event.root  ${ANADIR}/Sub_RUN_post_event_${i}.root
            mv channel_chip.pdf ${ANADIR}/.
            mv ${DATADIR}/SubRUN_${i}*root ${ANADIR}/.
            echo "RunNo. $r\t${ROC}*" >> ${ANADIR}/log
            #rm Spill_*GEMROC*.* -f
            #rm SubRUN_* -f
	fi
    done
    cd ${ANADIR}
    rm ${DATADIR}/event.root
    hadd -f ${ANADIR}/event.root ${ANADIR}/Sub_RUN_post_event*root
    cd -
    rm *.dat
fi

echo "Finish"

cd $HERE
