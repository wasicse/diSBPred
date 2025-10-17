#! /bin/bash
destinationFolder="diSBPred"
mkdir -p tempDownloads
cd tempDownloads
file="../Largefile.txt"
while read -r line; do
    echo -e "$line\n"
    line2="${line:2}"
    if [ ! -f "../$line2" ] 
    then
        echo "File does not exists! ../$line2"
        filename=$(basename "$line")
        wget -nc https://huggingface.co/datasets/wasicse/diSBPred/resolve/main/$filename.tar.gz -P ./
        tar -xzvf $filename.tar.gz
        mv $filename "../$line2" 
    fi
done <$file

rm -rf Scripts
rm  -rf Tools
cd ..