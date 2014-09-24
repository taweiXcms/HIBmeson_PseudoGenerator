###Plain root .C to be run
CONFIGFILE="fitB.C"

###All the header/related files needed
TRANSFERFILE="fitB.C,utilities.h"

###Output file location
DESTINATION=/net/hisrv0001/home/tawei/HeavyFlavor_20131030/Common/Pseudo20140320/PseudoBp_y/fout

###Log file location and it's name
LOGDIR=/net/hisrv0001/home/tawei/HeavyFlavor_20131030/Common/Pseudo20140320/PseudoBp_y/logout
LOGNAME=pseudo

########################## Create subfile ###############################
dateTime=$(date +%Y%m%d%H%M)
fileCounter=0
INFILE=""
mkdir -p $DESTINATION
mkdir -p $LOGDIR

for i in {1..50}
do
#    fileCounter=$((fileCounter+1))
    fileCounter=$i

# make the condor file
cat > subfile <<EOF

Universe = vanilla
Initialdir = .
Executable = exec_condor.sh
+AccountingGroup = "group_cmshi.$(whoami)"
Arguments =  $CONFIGFILE $DESTINATION ${fileCounter}
GetEnv       = True
Input = /dev/null

# log files
Output       = $LOGDIR/$LOGNAME-${fileCounter}.out
Error        = $LOGDIR/$LOGNAME-${fileCounter}.err
Log          = $LOGDIR/$LOGNAME-${fileCounter}.log

should_transfer_files = YES
when_to_transfer_output = ON_EXIT
transfer_input_files = $TRANSFERFILE

Queue
EOF

############################ Submit ###############################

#cat subfile
condor_submit subfile
#mv subfile $LOGDIR/$LOGNAME-$dateTime-$fileCounter.subfile
mv subfile $LOGDIR/$LOGNAME-$fileCounter.subfile
done
echo "Submitted $fileCounter jobs to Condor."
