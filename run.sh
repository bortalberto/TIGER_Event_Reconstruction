#DATADIR="$HOME/data/raw_dat/$1"
#ANADIR="$HOME/data/raw_root/$1"
DATADIR="/dati/Data_CGEM_IHEP_Integration_2019/raw_dat/$1"
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

mv decode.root ${ANADIR}/Spill_$1_decode.root
mv ana.root ${ANADIR}/Spill_$1_ana.root
mv event.root ${ANADIR}/Spill_$1_event.root
mv channel_chip.pdf ${ANADIR}/.
echo "RunNo. $r\t${ROC}*" >> ${ANADIR}/log

rm Spill_*GEMROC*.* -f


cd ${ANADIR}
hadd -f ana.root *ana.root
hadd -f event.root *event.root

#cp $HERE/check*cxx .
#root -b -q check_L1.cxx
#root -b -q check_L2.cxx
rm -f check*cxx

echo "Finish"

cd $HERE
