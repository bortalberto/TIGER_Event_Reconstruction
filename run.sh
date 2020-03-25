ANADIR="/dati/Data_CGEM_IHEP_Integration_2019/raw_root/$1"
HERE=$PWD
NROC=12
NSUB=2000

if [ ! -d ${ANADIR} ]
then
    mkdir ${ANADIR}
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
    #Decode
    ts -f bash -c "source run_decode.sh $1"
    #Merge decode
    ts -fd bash -c "source run_dec_merge.sh $1"
    #Reconstruction
    ts -fd bash -c "source run_recon.sh $1"
    #Merge recon
    #ts -fd $exe_ter -C $1
fi

echo "Finish"

cd $HERE
