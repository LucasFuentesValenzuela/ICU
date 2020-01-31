#!/bin/bash

#AT='~/routing/Build/Devel/Launchers/./AssignTraffic'
AT='./AssignTraffic'
FLAGS='-v -so -a ch -ord input'
NUM_ITER=100
DUMMY='-dummy 1351' #to be removed for elastic demand!
PARAMS='-eco 0.8 -f custom_bpr'

GRAPH_PATH='~/ASL/Code/amod_fw/amod_fw/data/elastic_Pb1/exp/1_0/graph.gr.bin'
OD_PATH='~/ASL/Code/amod_fw/amod_fw/data/elastic_Pb1/exp/1_0/od.csv'
OUTPUT_PATH='~/ASL/Code/TestFiles/output'
FLOW_PATH='~/ASL/Code/TestFiles/flow'

$AT $FLAGS -n $NUM_ITER $DUMMY $PARAMS -i GRAPH_PATH -od OD_PATH -o OUTPUT_PATH -fp FLOW_PATH
