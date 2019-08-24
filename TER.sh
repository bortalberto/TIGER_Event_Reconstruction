#!/bin/sh
#
# Examlple of using options in scripts
#
QUI=$PWD
NROC=12

if [ $# -eq 0 ]
then
    echo "Missing options!"
    echo "(run $0 -h for help)"
    echo ""
    exit 0
fi

run_number=$2
subrun_number=$3
roc_number=$4
DATADIR="/dati/Data_CGEM_IHEP_Integration_2019/raw_dat/RUN_$run_number"
ANADIR="/dati/Data_CGEM_IHEP_Integration_2019/raw_root/$run_number"

DATA_RAW_TIGER=$ANADIR
DATA_RAW_GRAAL="/dati/Data_CGEM_IHEP_Integration_2019/raw_graal"

echo RUN:    $run_number
echo SUBRUN: $subrun_number
echo ROC:    $roc_number
#echo $DATADIR

ECHO="false"
OPT_DEC="false"
OPT_DEC_I="false"
OPT_DEC_MERGE="false"
OPT_ANA="false"
OPT_EVT="false"
OPT_POS="false"
OPT_ALL="false"
OPT_all="false"
OPT_MAKE_ALL="false"
OPT_MAKE="false"
OPT_ROOT_FIRST="false"
OPT_ROOT_DEC="false"
OPT_ROOT_ANA="false"
OPT_ROOT_EVT="false"
OPT_ROOT_POS="false"
OPT_ROOT="false"
OPT_COPY="false"
while getopts "DAEPGMPQFmdawphegfCV" OPTION; do
    case $OPTION in

	w)
	    ECHO="true"
	    ;;

	h)
	    echo "Usage:"
	    echo ""
	    echo "   -w                  to execute echo \"hello world\""
	    echo "   -h                  help (this output)"
	    echo "   -F RUN SUBRUN ROC   Decoding one ROC  -> run Decode.py"
	    echo "   -D RUN SUBRUN       Decoding all ROCs -> run Decode.py"
	    echo "   -V RUN SUBRUN       Merge thr decoded files"
	    echo "   -A RUN SUBRUN       Calibration       -> run ana.C"
	    echo "   -E RUN SUBRUN       Create event      -> run event.C"
	    echo "   -Q RUN SUBRUN       Associate the time of each TIGER to its TP -> run post_event.C"
	    #echo "   -G RUN SUBRUN       Run decode ana event post_event"
	    echo "   -P RUN SUBRUN       RUN ana event"
	    echo "   -M                  make clean all"
	    echo "   -m                  make"
	    echo "   -f RUN SUBRUN ROC   open the decoded      root file for the run and subrun and roc given"
	    echo "   -d RUN SUBRUN       open the decoded      root file for the run and subrun given"
	    echo "   -a RUN SUBRUN       open the ana          root file for the run and subrun given"
	    echo "   -e RUN SUBRUN       open the event        root file for the run and subrun given"
	    echo "   -p RUN SUBRUN       open the post_event   root file for the run and subrun given"
	    echo "   -g RUN              open the merged event root file for the run"
	    echo "   -C RUN              copy the run into thr GRAAL folder"
	    exit 0
	    ;;

	E)
		OPT_EVT="true"
		;;

	D)
		OPT_DEC="true"
		;;

        V)
                OPT_DEC_MERGE="true"
                ;;

	F)
	        OPT_DEC_I="true"
		;;

	A)
		OPT_ANA="true"
		;;

        P)
                OPT_all="true"
                ;;

	G)
		OPT_ALL="true"
		;;		

	M)
		OPT_MAKE_ALL="true"
		;;		

	m)
		OPT_MAKE="true"
		;;		

	Q)
		OPT_POS="true"
		;;	

	d)
		OPT_ROOT_DEC="true"
		;;			

	e)
		OPT_ROOT_EVT="true"
		;;		

	a)
		OPT_ROOT_ANA="true"
		;;			

	p)
		OPT_ROOT_POS="true"
		;;	

        g)
                OPT_ROOT="true"
                ;;

	f)
	        OPT_ROOT_FIRST="true"
	       ;;

	C)
	        OPT_COPY="true"
	       ;;
    esac
done

#MAKE 
if [ $OPT_MAKE_ALL = "true" ]
then
    echo "make clean all";
    cd $TER;
    make clean all;
    cd $QUI;
fi 
#MAKE ALL
if [ $OPT_MAKE = "true" ]
then
    echo "make";
    cd $TER;
    make;
    cd $QUI;
fi
#Dummy 
if [ $ECHO = "true" ]
then
    echo "Hello world";
fi
#Decoding one roc
if [ $OPT_DEC_I = "true" ]
then
    if [ -z $roc_number ]; then echo "Use the command 'TER -F RUN SUBRUN ROC'"; exit; fi
    echo "Decoder $run_number $subrun_number $roc_number";
    cd $TER;
    python Decode.py ${DATADIR}/SubRUN_${subrun_number}_GEMROC_${roc_number}_TM.dat ${roc_number} 1 1
    mv ${DATADIR}/SubRUN_${subrun_number}_GEMROC_${roc_number}_TM.root ${ANADIR}/.
    cd $QUI;
fi
#Decoding all rocs
if [ $OPT_DEC = "true" ]
then
    if [ -z $subrun_number ]; then echo "Use the command 'TER -D RUN SUBRUN ROC'"; exit; fi
    for iroc in $(seq 0 $NROC);
    do
	if [ -f "${DATADIR}/SubRUN_${subrun_number}_GEMROC_${iroc}_TM.dat" ]; then
	    $exe_ter -F $run_number $subrun_number $iroc
	fi
    done
fi
#Merge the decoded files
if [ $OPT_DEC_MERGE = "true" ]
then
    if [ -z $subrun_number ]; then echo "Use the command 'TER -V RUN SUBRUN ROC'"; exit; fi
    if [ -f "${DATADIR}/SubRUN_${subrun_number}_GEMROC_0_TM.dat" ]; then
	hadd -f ${ANADIR}/Sub_RUN_dec_${subrun_number}.root ${ANADIR}/SubRUN_${subrun_number}_GEMROC*root
    fi
fi

#Event
if [ $OPT_EVT = "true" ]
then
    if [ -z $subrun_number ]; then echo "Use the command 'TER -E RUN SUBRUN'"; exit; fi
    echo "./bin/event $run_number $subrun_number";
    cd $TER;
    ./bin/event $run_number $subrun_number;
    cd $QUI;
fi
#Ana
if [ $OPT_ANA = "true" ]
then
    if [ -z $subrun_number ]; then echo "Use the command 'TER -A RUN SUBRUN'"; exit; fi
    echo "./bin/ana $run_number $subrun_number";
    cd $TER;
    ./bin/ana $run_number $subrun_number
    cd $QUI;
fi
#Post_Event
if [ $OPT_POS = "true" ]
then
    if [ -z $subrun_number ]; then echo "Use the command 'TER -P RUN SUBRUN'"; exit; fi
    echo "./bin/post_event $run_number $subrun_number";
    cd $TER;
    ./bin/post_event $run_number $subrun_number
     cd $QUI;
fi
#Dec+Ana+Evt
if [ $OPT_ALL = "true" ]
then
    if [ -z $roc_number ]; then echo "Use the command 'TER -G RUN SUBRUN ROC'"; exit; fi
    cd $TER;
    echo "Decoder $run_number $subrun_number $roc_number";
    python Decode.py ${DATADIR}/SubRUN_${$subrun_number}_GEMROC_${$roc_number}_TM.dat ${$roc_number} 1 1
    echo "./bin/ana $run_number $subrun_number";
    ./bin/ana $run_number $subrun_number
    echo "./bin/event $run_number $subrun_number";
    ./bin/event $run_number $subrun_number;
    #echo "./bin/post_event $run_number $subrun_number";
    #./bin/post_event $run_number $subrun_number
    cd $QUI;
fi

#Ana+Evt+Pos                                                                                                                                                                                                           
if [ $OPT_all = "true" ]
then
    if [ -z $subrun_number ]; then echo "Use the command 'TER -G RUN SUBRUN'"; exit; fi
    cd $TER;
    echo "./bin/ana $run_number $subrun_number";
    ./bin/ana $run_number $subrun_number
    echo "./bin/event $run_number $subrun_number";
    ./bin/event $run_number $subrun_number;
    #echo "./bin/post_event $run_number $subrun_number";
    #./bin/post_event $run_number $subrun_number
    cd $QUI;
fi

#dec roc root file open
if [ $OPT_ROOT_FIRST = "true" ]
then
    if [ -z $roc_number ]; then echo "Use the command 'TER -f RUN SUBRUN ROC'"; exit; fi
    root -l /home/ihep_data/data/raw_root/$run_number/SubRUN_${subrun_number}_GEMROC_${roc_number}_TM.root
fi


#dec root file open
if [ $OPT_ROOT_DEC = "true" ]
then
    if [ -z $subrun_number ]; then echo "Use the command 'TER -d RUN SUBRUN'"; exit; fi
    root -l /home/ihep_data/data/raw_root/$run_number/Sub_RUN_dec_$subrun_number.root
fi
#ana root file open
if [ $OPT_ROOT_ANA = "true" ]
then
    if [ -z $subrun_number ]; then echo "Use the command 'TER -a RUN SUBRUN'"; exit; fi
    root -l /home/ihep_data/data/raw_root/$run_number/Sub_RUN_ana_$subrun_number.root
fi
#evt root file open 
if [ $OPT_ROOT_EVT = "true" ]
then
    if [ -z $subrun_number ]; then echo "Use the command 'TER -e RUN SUBRUN'"; exit; fi
    root -l /home/ihep_data/data/raw_root/$run_number/Sub_RUN_event_$subrun_number.root
fi
#post_event root file open  
if [ $OPT_ROOT_POS = "true" ]
then
    if [ -z $subrun_number ]; then echo "Use the command 'TER -p RUN SUBRUN'"; exit; fi
    root -l /home/ihep_data/data/raw_root/$run_number/Sub_RUN_post_event_$subrun_number.root
fi
#merged event root file open
if [ $OPT_ROOT = "true" ]
then
    if [ -z $run_number ]; then echo "Use the command 'TER -g RUN'"; exit; fi
    root -l /home/ihep_data/data/raw_root/$run_number/event.root
fi
#copy the event foot into graal raw data folder
if [ $OPT_COPY = "true" ]
then
    if [ -z $run_number ]; then echo "Use the command 'TER -C RUN'"; exit; fi
    hadd -f ${ANADIR}/event.root ${ANADIR}/Sub_RUN_event*root
    if [ $run_number -lt 100 ]
    then
        echo $DATA_RAW_TIGER/event.root $DATA_RAW_GRAAL/run7163920$run_number.root
    fi
    if [ $run_number -gt 99 ]
        then
        echo $DATA_RAW_TIGER/event.root $DATA_RAW_GRAAL/run716392$run_number.root
	cp $DATA_RAW_TIGER/event.root $DATA_RAW_GRAAL/run716392$run_number.root
    fi

fi
