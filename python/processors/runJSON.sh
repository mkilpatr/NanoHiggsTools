#!/bin/bash
## declare an array variable
declare -a arr=("dyll" "diboson" "wjets" "ggHHto2b2tau" "ggHto2tau" "vbfHto2tau")

## now loop through the above array
for i in "${arr[@]}"
do

   echo "$i"
   ls /eos/uscms/store/user/mkilpatr/13TeV/nanoaod_2018_skim_diHiggs_033021/${i}_tree.root > testFile_json.txt
   cat testFile_json.txt
   python Stop0l_postproc_diHiggs.py -i testFile_json.txt -e 2018 -p json -s ${i}_100K -m 100000
done

