DATADIR="~/data/raw_dat/$1"
ANADIR="~/data/raw_root/$1"
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
for DATA in `cat ${DATADIR}/file_list`
do
    cp ${DATADIR}/${DATA}GEMROC_*_TM.dat .

    for ((g=0; g<=11; g++))
    #for ((g=0; g<=3; g++))
    #for ((g=4; g<=11; g++))
    do
	if [ -f "${DATA}GEMROC_${g}_TM.dat" ]; then
	   python Decode.py ${DATA}GEMROC_${g}_TM.dat $g 1 $r
	fi
    done

    hadd decode.root ${DATA}GEMROC*root
    #root -b -q ana.cxx
    ./bin/ana
    #root -b -q event.cxx
    ./bin/event
    
    mv decode.root ${ANADIR}/${DATA}decode.root
    mv ana.root ${ANADIR}/${DATA}ana.root
    mv event.root ${ANADIR}/${DATA}event.root
    echo "RunNo. $r\t${DATA}*" >> ${ANADIR}/log

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
pwd
ls

#mv event.root $1_event.root
cd $HERE
#cd ..
