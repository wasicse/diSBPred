#! /bin/bash
destinationFolder="diSBPred"

file="Largefile.txt"
while read -r line; do
    echo -e "$line\n"
    filename=$(basename "$line")
    wget https://www.cs.uno.edu/~mkabir3/$destinationFolder/$filename -P $line 
done <$file