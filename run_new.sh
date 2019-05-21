DATADIR="$HOME/data/raw_dat/$1"
ANADIR="$HOME/data/raw_root/$1"
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

if [ ! -f "${DATADIR}/file_list" ]; then
    echo "cannot find \"file_list\" in ${DATADIR}"
    return;
fi

rm Spill_*GEMROC*_TM.* -f
r=1
for ROC in `ls ${DATADIR} | sed 's/^.\{,33\}//' | sed 's/.$//' | sed 's/.$//'| sed 's/.$//' | sed 's/.$//'  | sed 's/.$//' | sed 's/.$//' | sed 's/.$//'`
do
    #echo ${DATADIR}/Spill_$1_GEMROC_*_TM.dat
    cp ${DATADIR}/Spill_$1_GEMROC_*_TM.dat .

    #echo Spill_$1_GEMROC_${ROC}_TM.dat
    if [ -f "Spill_$1_GEMROC_${ROC}_TM.dat" ]; then
	python Decode.py Spill_$1_GEMROC_${ROC}_TM.dat ${ROC} 1 $r
    fi

    #ls

    hadd -f decode.root Spill_$1_GEMROC*root

    ./bin/ana
    ./bin/event

    mv decode.root ${ANADIR}/Spill_$1_decode.root
    mv ana.root ${ANADIR}/Spill_$1_ana.root
    mv event.root ${ANADIR}/Spill_$1_event.root
    echo "RunNo. $r\t${ROC}*" >> ${ANADIR}/log

    rm Spill_*GEMROC*_TM.* -f
    let r+=1
done

cd ${ANADIR}
hadd -f ana.root *ana.root
hadd -f event.root *event.root

cp $HERE/check*cxx .
root -b -q check_L1.cxx
root -b -q check_L2.cxx
rm -f check*cxx

echo "Finish"

cd $HERE
