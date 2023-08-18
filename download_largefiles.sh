#! /bin/bash
destinationFolder="diSBPred"
mkdir -p tempDownloads
cd tempDownloads
file="../Largefile.txt"
while read -r line; do
    echo -e "$line\n"
    line2="${line:2}"
    if [ ! -d "../$line2" ] 
    then
        filename=$(basename "$line")
        wget -nc https://www.cs.uno.edu/~mkabir3/$destinationFolder/$filename.tar.gz -P ./
        tar -xzvf $filename.tar.gz
        mv $line "../$line2" 
    fi
done <$file

rm -rf Scripts
rm  -rf Tools