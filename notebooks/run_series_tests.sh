#!/bin/bash

#default parameters values
NI=5000
LR=10000
NO=50



PATH='25Nodes'
SD='rel_prog/ni_comparison'

for ni in 10 100 1000 10000; do $CONDA_PYTHON_EXE run_FW_outer.py -p $PATH -L $LR -ni $ni -no $NO -sd $SD; done
# $CONDA_PYTHON_EXE generate_figure.py -p $PATH -pc $SD

PATH='25Nodes'
SD='rel_prog/L_comparison'

for L in 10 100 1000 10000; do $CONDA_PYTHON_EXE run_FW_outer.py -p $PATH -L $L -ni $NI -no $NO -sd $SD; done
# $CONDA_PYTHON_EXE generate_figure.py -p $PATH -pc $SD

# #EXPERIMENT 2
# PATH='25Nodes'
# SD='ni_comparison'
# NO=20

# for ni in 10 100 1000 10000; do $CONDA_PYTHON_EXE run_FW_outer.py -p $PATH -L $L -ni $ni -no $NO -sd $SD; done
# $CONDA_PYTHON_EXE generate_figure.py -p $PATH -pc $SD

# #EXPERIMENT 3
# PATH='100Nodes'
# SD='L_comparison'
# NO=15 #YOU CERTAINLY NEED MORE OPERATIONS FOR 100 NODES

# for L in 1000 10000; do $CONDA_PYTHON_EXE run_FW_outer.py -p $PATH -L $L -ni $NI -no $NO -sd $SD; done
# $CONDA_PYTHON_EXE generate_figure.py -p $PATH -pc $SD

# #EXPERIMENT 4
# PATH='100Nodes'
# SD='ni_comparison'
# NO=15

# for ni in 10 100 1000; do $CONDA_PYTHON_EXE run_FW_outer.py -p $PATH -L $L -ni $ni -no $NO -sd $SD; done
# $CONDA_PYTHON_EXE generate_figure.py -p $PATH -pc $SD