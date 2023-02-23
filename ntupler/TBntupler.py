import pydrcTB
import os
import argparse
from glob import glob

# Command :
# For fastmode ntuple : python3 TBntuplize.py -rn (run number) -fe 100000 -f -o /path/to/output/fast/ntuple/
# For waveform ntuple : python3 TBntuplize.py -rn (run number) -fe 5000 -o /path/to/output/wave/ntuple/
# Remember to change output directory to your own SE_UserHome directory!

parser = argparse.ArgumentParser()

parser.add_argument('--runNumber', '-rn', action='store', type=str, required=True, help='Run number')
parser.add_argument('--maxEntry', '-me', action='store', type=int, default=-1, help='Maximum # of entries process, unlimited if -1')
parser.add_argument('--fileEntry', '-fe', action='store', type=int, default=-1, help='# of entries per single ntuple, unlimited if -1')
parser.add_argument('--doFast', '-f', action='store_true', default=False, help='If true, make fast mode ntuple')
parser.add_argument('--output', '-o', action='store', type=str, default="./", help='Directory to write output ntuples')
args = parser.parse_args()

rn=args.runNumber
mode = "Fast" if args.doFast else "Wave"

FileNumbers = []
for mid in range(1,16):
    DataPath=f"/gatbawi/dream/data/HDD_Run_{rn}_validated/Run_{rn}_{mode}/Run_{rn}_{mode}_MID_{mid}/"
    FileNumbers.append( len(glob(DataPath + "*.dat")) )
    maxFileNum = min(FileNumbers)

DataFiles = []
for fn in range(maxFileNum) :
    FilesPerFN = []
    for mid in range(1,16) :
        fileName=f"/gatbawi/dream/data/HDD_Run_{rn}_validated/Run_{rn}_{mode}/Run_{rn}_{mode}_MID_{mid}/Run_{rn}_{mode}_MID_{mid}_FILE_{fn}.dat"
        if(os.path.isfile(fileName)) :
            FilesPerFN.append(fileName)
        else:
            print(fileName,"not found")
    DataFiles.append(FilesPerFN)

if args.output.endswith("/") : 
    outputPath = args.output + f"Run_{rn}/{mode}/"
else :
    outputPath = args.output + f"/Run_{rn}/{mode}/"
if not os.path.exists(outputPath) :
    os.makedirs(outputPath)
output_name = outputPath + f"ntuple_Run_{rn}_{mode}" # ex) For fast mode, ntuple_Run_393_Fast_0.root, ntuple_Run_393_Fast_1.root ... will be stored in directory named "args.output + /Run_393/Fast/"

reader = pydrcTB.TBread()
if args.doFast : 
    reader.ntuplizeFastmode(DataFiles, output_name, args.maxEntry, args.fileEntry)
else :
    reader.ntuplizeWaveform(DataFiles, output_name, args.maxEntry, args.fileEntry)