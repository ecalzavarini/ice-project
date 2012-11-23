#!/bin/bash

if [ -z $1 ]
then
echo "Running with local"
else
echo "Running exemple" $1
mv param.in .param.in
mv define.h .define.h
namedir="../exemples/"$1
if [ -d $namedir ]
then
echo "The exemple" $1 "will be computed"
cp $namedir/param.in ./
cp $namedir/define.h ./
else
echo "The exemple" $1 "does not exist"
exit 
fi
fi

#STARTING
echo "Cleaning"
make clean
rm ./*.dat
rm -rf RUN
echo "Compiling"
make
echo "Running"
./ice2d


if [ -z $1 ]
then
echo "Local ended"
else
echo "Ended exemple" $1
mv .param.in param.in
mv .define.h define.h
fi

