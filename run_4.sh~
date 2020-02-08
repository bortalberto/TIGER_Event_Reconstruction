ANADIR="/dati/Data_CGEM_IHEP_Integration_2019/raw_root/$1"
HERE=$PWD
NROC=12
NSUB=1000

if [ ! -d ${ANADIR} ]
then
    mkdir ${ANADIR}
fi

for i in $(seq 0 20);
do
    ts sleep 0.01
done
ts -N 50
ts -df sleep 0.01
    
#Merge Recon
ts -df bash -c "$exe_ter -C $1"
echo "Terminated all"

for i in $(seq 0 50);
do
    ts sleep 0.01
    done
ts -N 25
ts -df sleep 0.01
rm /tmp/ihep_data/ts-out*
cd $HERE

