#!/bin/bash
#SBATCH --partition=unlimited
#SBATCH --time=30-00:00:00
cd /work-zfs/lhc/lkang/CMSSW_10_2_5/src
eval `scramv1 runtime -sh`
cd /work-zfs/lhc/lkang/JHUGenerator.v7.4.0/JHUGenerator
eval $(../JHUGenMELA/MELA/setup.sh env)
./JHUGen Process=66 DecayMode1=0 DecayMode2=0 Interf=0 ReweightInterf=0 MPhotonCutoff=4 VegasNc0=10000000 DataFile=EW_offshell/decay_same_sign_Cza_Caa/output/EW_offshell_run MReso=125 GaReso=0.00407 mJJcut=30 m4l_min=124.98 m4l_max=125.02 detajetcut=0 JetsOppositeEta=0 etajetcut=6.5 FacScheme=1 MuFacMultiplier=1 RenScheme=1 MuRenMultiplier=1 ghz1=0,0 ghgsgs2=0.051797691,0 ghzgs2=0.047047878,0 LHAPDF=NNPDF30_nlo_as_0118/NNPDF30_nlo_as_0118.info pTjetcut=15 deltaRcut=0.3 VBFoffsh_run=$SLURM_ARRAY_TASK_ID

# For the second step, simply add "VegasNc2=100000 ReadCSmax to the JHUGen inputs"
