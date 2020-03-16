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
    echo "Start reconstruction"
    DATADIR="/dati/Data_CGEM_IHEP_Integration_2019/raw_dat/RUN_$1"
    ANADIR="/dati/Data_CGEM_IHEP_Integration_2019/raw_root/$1"
    ts -S $MAXJOB
    ts -N $MAXJOB sleep 0.01
    ts -df sleep 0.01
    #Reconstruction
    for i in $(seq 0 $NSUB);
    do
	ts bash -c "$exe_ter -P $1 $i"
    done
    ts -N $MAXJOB sleep 0.01
    ts -df sleep 0.01
    echo "Terminate reconstruction"
    rm -f /tmp/ihep_data/ts-out.*
    source run_4.sh

fi
cd $HERE

