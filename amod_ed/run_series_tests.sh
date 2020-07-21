#!/bin/bash


#TODO: include tolerances as parameters too

#default parameters values
NI=5000
LR=10000
NO=50
SD='Li_comparison_ni_5000'

DATA_PATH='3Nodes'
for Li in 10 100 1000 10000 100000 1000000 10000000
do 
$CONDA_PYTHON_EXE run_FW_outer.py -p $DATA_PATH -L $Li -ni $NI -no $NO -sd $SD -ev 0 -sc 'rp' -fu 0.0
done
# $CONDA_PYTHON_EXE test.py


# PATH='25Nodes'
# for Li in 1000 100000 10000000
# do 
# $CONDA_PYTHON_EXE run_FW_outer.py -p $PATH -L $Li -ni $NI -no $NO -sd $SD -ev 0 -sc 'rp' -fu 0.0
# done

# PATH='100Nodes'
# for Li in 10 100 1000 10000 100000 1000000 10000000
# do 
# $CONDA_PYTHON_EXE run_FW_outer.py -p $PATH -L $Li -ni $NI -no $NO -sd $SD -ev 0 -sc 'rp' -fu 0.0
# done
