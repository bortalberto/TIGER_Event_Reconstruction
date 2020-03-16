HERE=$PWD
NROC=$(tail $TER/setting/NROC)
NSUB=$(tail $TER/setting/NSUB)
MAXJOB=$(tail $TER/setting/MAXJOB)
TS_MAXCONN=100
export TS_MAXCONN=100

if [ ! -d ${ANADIR} ]
then
    mkdir ${ANADIR}
fi


if [[ $1 -gt 88 ]]; then
    #Init
    echo "Start merging decode"
    DATADIR="/dati/Data_CGEM_IHEP_Integration_2019/raw_dat/RUN_$1"
    ANADIR="/dati/Data_CGEM_IHEP_Integration_2019/raw_root/$1"
    ts -S $MAXJOB
    #Merge
    for i in $(seq 0 $NSUB);
    do
        ts bash -c "$exe_ter -V $1 $i"
    done
    #Wait
    ts -N $MAXJOB sleep 0.01
    ts -df sleep 0.01
    echo "Terminate merging decoded"
    rm -f /tmp/ihep_data/ts-out.*
    source run_3.sh
fi
cd $HERE

