ANADIR="/dati/Data_CGEM_IHEP_Integration_2019/raw_root/$1"
HERE=$PWD
NROC=14
NSUB=2000

if [ ! -d ${ANADIR} ]
then
    mkdir ${ANADIR}
fi


if [[ $1 -gt 88 ]]; then
    DATADIR="/dati/Data_CGEM_IHEP_Integration_2019/raw_dat/RUN_$1"
    ts -S 50
    #Decode
    for i in $(seq 0 $NSUB);
    do
        ts bash -c "$exe_ter -D $1 $i"
    done
    ts -N 50 sleep 0.01
    ts -df sleep 0.01
    echo "Terminated decoding"
    source run_2.sh
fi
cd $HERE

