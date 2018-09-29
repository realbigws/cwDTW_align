#!/bin/bash

executable=cwDTW_align

if [ -f "$executable" ] 
then
	echo " executable files '$executable' already compiled "
	exit 1
fi


mkdir -p Release
cd Release
cmake -DCMAKE_BUILD_TYPE=Release ..
make
mv bin/wletdtw ../$executable
cd ../
rm -rf Release


