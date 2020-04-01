#!/bin/bash


#TODO: include tolerances as parameters too

#default parameters values
NI=100
LR=1000000
NO=20
SD='unstucking_balance'

PATH='10Nodes'
$CONDA_PYTHON_EXE run_FW_outer.py -p $PATH -L $LR -ni $NI -no $NO -sd $SD -ev 0 -sc 'rp' -fu 0.0
