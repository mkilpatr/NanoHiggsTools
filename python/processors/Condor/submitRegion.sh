#!/bin/bash

location=$1

python SubmitLPC.py -f ../Stop0l_postproc_diHiggs.py -c ../sampleSets_PostProcessed_2018.cfg -o /store/user/mkilpatr/13TeV/${location}_tot/ -e 2018 -p json
python SubmitLPC.py -f ../Stop0l_postproc_diHiggs.py -c ../sampleSets_PostProcessed_2018.cfg -o /store/user/mkilpatr/13TeV/${location}_emu/ -e 2018 -p json -r emu
python SubmitLPC.py -f ../Stop0l_postproc_diHiggs.py -c ../sampleSets_PostProcessed_2018.cfg -o /store/user/mkilpatr/13TeV/${location}_ehad/ -e 2018 -p json -r ehad
python SubmitLPC.py -f ../Stop0l_postproc_diHiggs.py -c ../sampleSets_PostProcessed_2018.cfg -o /store/user/mkilpatr/13TeV/${location}_muhad/ -e 2018 -p json -r muhad
python SubmitLPC.py -f ../Stop0l_postproc_diHiggs.py -c ../sampleSets_PostProcessed_2018.cfg -o /store/user/mkilpatr/13TeV/${location}_hadhad/ -e 2018 -p json -r hadhad
