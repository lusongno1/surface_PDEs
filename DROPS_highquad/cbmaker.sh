#! /bin/bash
rm ./bin -rf
mkdir bin
cd bin
cmake ../src -G "CodeBlocks - Unix Makefiles"
cd ..
cp ./svrun.sh ./bin/surfactant
