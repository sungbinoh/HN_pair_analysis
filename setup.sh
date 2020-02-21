#!/bin/bash
export PLOTTER_WORKING_DIR=`pwd`
export path_2016='/data6/Users/suoh/SKFlatOutput/Run2Legacy_v4/HN_pair_all/2016/'
export path_2017='/data6/Users/suoh/SKFlatOutput/Run2Legacy_v4/HN_pair_all/2017/'
export path_2018='/data6/Users/suoh/SKFlatOutput/Run2Legacy_v4/HN_pair_all/2018/'
export FILE_PATH=$PLOTTER_WORKING_DIR/rootfiles/
export PLOT_PATH=$PLOTTER_WORKING_DIR/plots/
export SCRIPT_DIR=$PLOTTER_WORKING_DIR/script/
export ROOT_INCLUDE_PATH=$ROOT_INCLUDE_PATH:$PLOTTER_WORKING_DIR/include/:$PLOTTER_WORKING_DIR/src/

source $PLOTTER_WORKING_DIR/bin/BashColorSets.sh

alias lqout='cd '$PLOT_PATH

mkdir -p $FILE_PATH/$CATANVERSION
