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
    echo "Start Decode"
    #Init
    DATADIR="/dati/Data_CGEM_IHEP_Integration_2019/raw_dat/RUN_$1"
    ANADIR="/dati/Data_CGEM_IHEP_Integration_2019/raw_root/$1"
    ts -S $MAXJOB
    #Clean 
    cd /dati/Data_CGEM_IHEP_Integration_2019/raw_dat
    python purger.py
    rm -f /dati/Data_CGEM_IHEP_Integration_2019/raw_root/$1/*
    rm -f /dati/Data_CGEM_IHEP_Integration_2019/raw_graal/run716392$1.root
    rm -f /dati/Data_CGEM_IHEP_Integration_2019/rec_graal/rec_run716392$1.root
    rm -f /dati/Data_CGEM_IHEP_Integration_2019/evt_graal/evt_run716392$1.root
    #Decode
    cd $TER
    for i in $(seq 0 $NSUB);
    do
        ts bash -c "$exe_ter -D $1 $i"
    done
    #Wait
    ts -N $MAXJOB sleep 0.01
    ts -df sleep 0.01
    echo "Terminate decode"
    rm -f /tmp/ihep_data/ts-out.*
    source run_2.sh
fi
cd $HERE

