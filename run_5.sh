ANADIR="/dati/Data_CGEM_IHEP_Integration_2019/raw_root/$1"
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

echo "Start clean bad subruns"
ts -S $MAXJOB
#Clean
ts -df bash -c "$exe_ter -c $1"
ts -N $MAXJOB sleep 0.01
ts -df sleep 0.01
echo "Terminate clean bad subruns"
rm -f /tmp/ihep_data/ts-out.*
source run_6.sh

