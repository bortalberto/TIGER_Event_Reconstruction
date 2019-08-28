ANADIR="/dati/Data_CGEM_IHEP_Integration_2019/raw_root/$1"
HERE=$PWD
NROC=12
NSUB=1000

if [ ! -d ${ANADIR} ]
then
    mkdir ${ANADIR}
fi


if [[ $1 -gt 88 ]]; then
    DATADIR="/dati/Data_CGEM_IHEP_Integration_2019/raw_dat/RUN_$1"
    ts -S 20
    for i in $(seq 0 20);
    do
    ts sleep 0.01
    done
    ts -N 20 sleep 0.01
    ts -df sleep 0.01
    #Merge
    for i in $(seq 0 $NSUB);
    do
        ts bash -c "$exe_ter -V $1 $i"
    done
    ts -N 20 sleep 0.01
    ts -df sleep 0.01
    echo "Terminated merging decoded"
    source run_3.sh

fi
cd $HERE

