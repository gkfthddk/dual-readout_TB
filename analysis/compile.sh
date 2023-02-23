#! /bin/bash

# To compile TBxxx.cc, type the commnand below
# 
# bash compile.sh TBxxx.cc (or TBxxx for short)
#
# This will generate TBxxx.exe executiable

ext=${1##*.}
fname=`basename ${1} .${ext}`

echo "Compiling $fname.cc to $fname.exe"
g++ -I$INSTALL_DIR_PATH/include -L$INSTALL_DIR_PATH/lib64 $INSTALL_DIR_PATH/lib64/libdrcTB.so `root-config --cflags --libs` ${fname}.cc -o ${fname}.exe
echo "Done!"