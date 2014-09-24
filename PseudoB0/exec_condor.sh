CONFIGFILE=$1
DESTINATION=$2
SEED=$3

source /osg/app/cmssoft/cms/cmsset_default.sh
export SCRAM_ARCH=slc5_amd64_gcc462
cd /cvmfs/cms.cern.ch/slc5_amd64_gcc462/cms/cmssw-patch/CMSSW_5_3_2_patch4/src
#cd /cvmfs/cms.cern.ch/slc5_amd64_gcc434/cms/cmssw-patch/CMSSW_5_0_1_patch3/src
#eval `scramv1 runtime -sh` 
cmsenv
cd -

gcc --version
root-config --version
cat /proc/version                       

mkdir ResultsBzero
root -l -b -q  $CONFIGFILE\($SEED\)
#root -l -b -q  $CONFIGFILE++\($SEED\)
rm -r ${DESTINATION}/ResultsBzero-Pseudo-${SEED}
mv ResultsBzero ${DESTINATION}/ResultsBzero-Pseudo-${SEED}
cp ${DESTINATION}/ResultsBzero-Pseudo-${SEED}/SigmaBzero.root ${DESTINATION}/../ToMerge/SigmaBzero_${SEED}.root
