#! /bin/bash
destinationFolder="diSBPred"

mkdir -p Database
if [ ! -d "./Database/nr" ] 
then
    wget -nc https://huggingface.co/wasicse/featuresExtractfiles/resolve/main/nr.tar.gz -P ./Database/
    cd Database
    tar -xvf nr.tar.gz
fi
rm -rf nr.tar.gz    