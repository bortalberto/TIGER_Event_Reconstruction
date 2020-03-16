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
echo "Start merge good subruns"
ts -S $MAXJOB
#Merge
ts -df bash -c "$exe_ter -C $1"
echo "Terminate merge good subruns"
