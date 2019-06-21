ANADIR="/dati/Data_CGEM_IHEP_Integration_2019/raw_root/$1"
HERE=$PWD
NROC=11
NSUB=10

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
    count=0
    for i in $(seq 0 $NSUB);
    do
	# Begin the loop on hte subrun
	for ROC in $(seq 0 $NROC);
	do
	    # Begin the loop on the ROC
            echo ROC: ${ROC}  
	    echo sub: ${i}
            if [ -f "${DATADIR}/SubRUN_${i}_GEMROC_${ROC}_TM.dat" ]; then
		ts bash -c "python Decode.py ${DATADIR}/SubRUN_${i}_GEMROC_${ROC}_TM.dat ${ROC} 1 $r" # triggermatch 
                #python Decode.py ${DATADIR}/SubRUN_${i}_GEMROC_${ROC}_TL.dat ${ROC} 1 $r # triggerless     
		count=`expr $count + 1`
	    fi
	done
    done
    ts -S 20
    ts -N 20 bash -c "sleep 0.0001"
    ts -w
    for i in $(seq 0 $NSUB);
    do
	echo "SubRUN_${i}_GEMROC_*_TM.root"
	if [ -f "${DATADIR}/SubRUN_${i}_GEMROC_0_TM.root" ]; then
	#if [[ $count -gt 0 ]]; then
	    ts bash -c "
	    echo 'Hello'
 	    hadd -f ${ANADIR}/Sub_RUN_dec_${i}.root ${DATADIR}/SubRUN_${i}_GEMROC*root
	    echo 'Begin ana'
	    ./bin/ana $1 $i
	    echo 'Begin evt'
	    ./bin/event $1 $i
	    ./bin/post_event $1 $i
            mv -f ${DATADIR}/SubRUN_${i}*root ${ANADIR}/. "
	    count=`expr $count + 1`
	fi
    done
    ts -S 20
    ts -N 20 bash -c "sleep 0.0001"
    ts -w
    hadd -f ${ANADIR}/event.root ${ANADIR}/Sub_RUN_post_event*root
fi

echo "Finish"

cd $HERE
