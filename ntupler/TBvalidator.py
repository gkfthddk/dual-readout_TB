import pydrcTB
import re
import os
import argparse
import ROOT
from ROOT import TH1F
from glob import glob


def atoi(text) :
    return int(text) if text.isdigit() else text

def natural_keys(text) :
    return [ atoi(c) for c in re.split(r'(\d+)', text) ]

parser = argparse.ArgumentParser()

parser.add_argument('--runNumber', '-rn', action='store', type=int, required=True, help="Run number")
parser.add_argument('--inputNtuple', '-in', action='store', type=str, default="./",help="Path to input ntuples")
parser.add_argument('--doFast', '-f', action='store_true', default=False, help='If true, validate in fastmode. Else, validate in waveform mode')
parser.add_argument('--output', '-o', action='store', type=str, default="./", help="Path to store output png files")

args = parser.parse_args()
rn = args.runNumber
mode = "Wave" if not args.doFast else "Fast"
print(f"Validating in {mode} mode with run number {rn}...")

DataFileNumbers = []
for mid in range(1,16):
    FullDataPath = f"/gatbawi/dream/data/HDD_Run_{rn}_validated/Run_{rn}_{mode}/Run_{rn}_{mode}_MID_{mid}/*.dat"
    FullData_list = [files for files in glob(FullDataPath) if not os.stat(files).st_size == 0]
    DataFileNumbers.append( len(FullData_list) )
maxDataFileNum = min(DataFileNumbers) # just in case if MIDs have different # of fast .dat files... which should not happen
DataFiles=[]
for fn in range(maxDataFileNum) :
    DataFilesPerFN = []
    for mid in range(1,16) :
        fileName= f"/gatbawi/dream/data/HDD_Run_{rn}_validated/Run_{rn}_{mode}/Run_{rn}_{mode}_MID_{mid}/Run_{rn}_{mode}_MID_{mid}_FILE_{fn}.dat"
        if( os.path.isfile(fileName) and (os.stat(fileName).st_size != 0) ) :
            DataFilesPerFN.append(fileName)
        else:
            print(fileName,"not found or has size of zero")
    DataFiles.append(DataFilesPerFN)

if not args.inputNtuple.endswith("/") :
    inputNtuplePath = args.inputNtuple + f"/Run_{rn}/{mode}/"
else :
    inputNtuplePath = args.inputNtuple + f"Run_{rn}/{mode}/"
Ntuples = []
for fn in glob(inputNtuplePath + f"ntuple_Run_{rn}_{mode}_*.root") :
    Ntuples.append(fn)
Ntuples.sort(key = natural_keys)

if args.output.endswith("/") :
    outDir = args.output + f"valid_plots_Run_{rn}/{mode}/"
else :
    outDir = args.output + f"/valid_plots_Run_{rn}/{mode}/"
if not os.path.exists(outDir) : 
    os.makedirs(outDir)

print("Data files : ", DataFiles)
print("Ntuples : ", Ntuples)

# Auxiliary detector channel IDs
cid_mc = pydrcTB.TBcid(1, 16)
cid_tc = pydrcTB.TBcid(1, 14)
cid_ps = pydrcTB.TBcid(1, 12)
cid_dwc1_r = pydrcTB.TBcid(1, 17)
cid_dwc1_l = pydrcTB.TBcid(1, 19)
cid_dwc1_u = pydrcTB.TBcid(1, 21)
cid_dwc1_d = pydrcTB.TBcid(1, 23)
cid_dwc2_r = pydrcTB.TBcid(1, 25)
cid_dwc2_l = pydrcTB.TBcid(1, 27)
cid_dwc2_u = pydrcTB.TBcid(1, 29)
cid_dwc2_d = pydrcTB.TBcid(1, 31)

# Module 1 S ch. channel IDs 
cid_M1T1_S = pydrcTB.TBcid(1, 1)
cid_M1T2_S = pydrcTB.TBcid(1, 3)
cid_M1T3_S = pydrcTB.TBcid(1, 7)
cid_M1T4_S = pydrcTB.TBcid(1, 5)
# Module 1 C ch. channel IDs
cid_M1T1_C = pydrcTB.TBcid(1, 9)
cid_M1T2_C = pydrcTB.TBcid(1, 11)
cid_M1T3_C = pydrcTB.TBcid(1, 15)
cid_M1T4_C = pydrcTB.TBcid(1, 13)

# Module 2 S ch. channel IDs 
cid_M2T1_S = pydrcTB.TBcid(2, 1)
cid_M2T6_S = pydrcTB.TBcid(2, 2)
cid_M2T2_S = pydrcTB.TBcid(2, 3)
cid_M2T7_S = pydrcTB.TBcid(2, 4)
cid_M2T3_S = pydrcTB.TBcid(2, 5)
cid_M2T8_S = pydrcTB.TBcid(2, 6)
cid_M2T4_S = pydrcTB.TBcid(2, 7)
cid_M2T9_S = pydrcTB.TBcid(2, 8)
# Module 2 C ch. channel IDs 
cid_M2T1_C = pydrcTB.TBcid(2, 9)
cid_M2T6_C = pydrcTB.TBcid(2, 10)
cid_M2T2_C = pydrcTB.TBcid(2, 11)
cid_M2T7_C = pydrcTB.TBcid(2, 12)
cid_M2T3_C = pydrcTB.TBcid(2, 13)
cid_M2T8_C = pydrcTB.TBcid(2, 14)
cid_M2T4_C = pydrcTB.TBcid(2, 15)
cid_M2T9_C = pydrcTB.TBcid(2, 16)

# Module 1 channel ID list (total 19 channels) 
m1_cid_list  = [cid_mc,     cid_tc,     cid_ps,
                cid_dwc1_r, cid_dwc1_l, cid_dwc1_u, cid_dwc1_d,
                cid_dwc2_r, cid_dwc2_l, cid_dwc2_u, cid_dwc2_d,
                cid_M1T1_S, cid_M1T2_S, cid_M1T3_S, cid_M1T4_S,
                cid_M1T1_C, cid_M1T2_C, cid_M1T3_C, cid_M1T4_C]
# Module 2 channel ID list (total 27 channels) 
m2_cid_list  = [cid_mc,     cid_tc,     cid_ps,
                cid_dwc1_r, cid_dwc1_l, cid_dwc1_u, cid_dwc1_d,
                cid_dwc2_r, cid_dwc2_l, cid_dwc2_u, cid_dwc2_d,
                cid_M2T1_S, cid_M2T2_S, cid_M2T3_S, cid_M2T4_S, cid_M2T6_S, cid_M2T7_S, cid_M2T8_S, cid_M2T9_S,
                cid_M2T1_C, cid_M2T2_C, cid_M2T3_C, cid_M2T4_C, cid_M2T6_C, cid_M2T7_C, cid_M2T8_C, cid_M2T9_C]

# For validating whole single run, do this :
TButils = pydrcTB.TButility()
TButils.loading("/gatbawi/dream/mapping/mapping_Aug2022TB.root")
validator = pydrcTB.TBvalid()
validator.setDataList(DataFiles)
validator.setNtupleList(Ntuples)
if (args.doFast) :
    validator.checkFastTrigNum()
else :
    validator.checkWaveTrigNum()
for cid in m1_cid_list : # Use m2_cid_list for Module 2 runs
    MID = cid.mid()
    ch = cid.channel()
    print(f"Validating MID : {MID} Ch : {ch}")
    cid = pydrcTB.TBcid(MID, ch)
    Data_hist_name = f"h_Data_MID{MID}Ch{ch}"
    Ntuples_hist_name = f"h_Ntuple_MID{MID}Ch{ch}"
    Ratio_hist_name = f"h_Ratio_MID{MID}Ch{ch}"
    if args.doFast :
        valid = validator.simpleValidFast(cid)
        if not valid : 
            Data_hist = validator.drawFastHistFromData(cid, Data_hist_name)
            Ntup_hist = validator.drawFastHistFromNtuple(cid, Ntuples_hist_name)
            validator.drawRatio(Data_hist, Ntup_hist, Ratio_hist_name, outDir)
    else :
        valid = validator.simpleValidWave(cid)
        if not valid :
            Data_hist = validator.drawWaveHistFromData(cid, Data_hist_name)
            Ntup_hist = validator.drawWaveHistFromNtuple(cid, Ntuples_hist_name)                
            validator.drawRatio(Data_hist, Ntup_hist, Ratio_hist_name, outDir)
