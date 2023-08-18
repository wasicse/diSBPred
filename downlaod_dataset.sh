#! /bin/bash
destinationFolder="diSBPred"

mkdir -p Database
if [ ! -d "./Database/nr" ] 
then
    wget -nc https://www.cs.uno.edu/~mkabir3/$destinationFolder/nr.tar.gz -P ./Database/
    cd Database
    tar -xvf nr.tar.gz
fi
