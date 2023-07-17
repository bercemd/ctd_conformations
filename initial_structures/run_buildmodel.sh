#!/bin/bash

export MODINSTALL10v1="path for MODELLER"
export KEY_MODELLER10v1="XXXXXXX" # the key for modeller
export LIBS_LIB10v1="$MODINSTALL10v1/modlib/libs.lib"
alias  mod="mod10.1"
export MODELLEREXEC="$MODINSTALL10v1/bin/mod10.1_x86_64-intel8"
export LD_LIBRARY_PATH="$MODINSTALL10v1/lib/x86_64-intel8":$LD_LIBRARY_PATH
export PATH="~/miniconda3/bin:$PATH"
export PYTHONPATH="~/miniconda3/lib/python3.7"
export PYTHONHOME="~/miniconda3/lib"

export MMTSBDIR="path for MMTSB"

export CHARMMEXEC='path for CHARMM'

#Generation of 2CTD model

$CHARMMEXEC < build.inp > build.out

convpdb.pl -chainfromseg ctd.pdb > ctd.chain.pdb
loopModel.pl -models 5 -loop A8:YSPTSPS ctd.chain.pdb
