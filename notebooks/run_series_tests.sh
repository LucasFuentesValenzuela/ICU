#!/bin/bash


#TODO: include tolerances as parameters too

#default parameters values
NI=1000
LR=10000
NO=20

# PATH='2Nodes'
SD='quals/ni_comparison'
# $CONDA_PYTHON_EXE run_FW_outer.py -p $PATH -L $LR -ni $NI -no $NO -sd $SD -ev 0 -sc 'rp' -fu 1 

# PATH='10Nodes'
# $CONDA_PYTHON_EXE run_FW_outer.py -p $PATH -L $LR -ni $NI -no $NO -sd $SD -ev 0 -sc 'rp' -fu 1

PATH='25Nodes'
for L in 10 100 1000 ; do $CONDA_PYTHON_EXE run_FW_outer.py -p $PATH -L $L -ni $NI -no $NO -sd $SD -ev 0 -sc 'rp' -fu 0; done


# PATH='100Nodes'
# $CONDA_PYTHON_EXE run_FW_outer.py -p $PATH -L $LR -ni $NI -no $NO -sd $SD -ev 0 -sc 'rp' -fu 1

# PATH='200Nodes'
# $CONDA_PYTHON_EXE run_FW_outer.py -p $PATH -L $LR -ni $NI -no $NO -sd $SD -ev 0 -sc 'rp' -fu 1