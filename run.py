import subprocess
import sys
with open("condor_script/sub_condor_temp.sub","r") as f:
  lines=f.readlines()
rns=sys.argv[1].split(",")
for rn in rns:
  with open("condor_script/sub_condor_{}.sub".format(rn),"w") as f:
    for line in lines:
      f.write(line.replace("RUN_NUMBER",rn))
  subprocess.call(["condor_submit","condor_script/sub_condor_{}.sub".format(rn)])
 
