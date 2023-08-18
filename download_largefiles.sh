#! /bin/bash
destinationFolder="diSBPred"
mkdir -p tempDownloads
cd tempDownloads
file="../Largefile.txt"
while read -r line; do
    echo -e "$line\n"
    filename=$(basename "$line")
    wget -nc https://www.cs.uno.edu/~mkabir3/$destinationFolder/$filename.tar.gz -P ./
    tar -xzvf $filename.tar.gz
    line2="${line:2}"
    mv $line "../$line2" 
done <$file

rm Scripts
rm Tools