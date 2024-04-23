#! /bin/bash

# Input arguments
input_fasta=$1
output_dir=$2

# if no input arguments
if [ -z "$input_fasta" ] || [ -z "$output_dir" ]
then
	echo "Please provide input fasta file and output directory"
	exit 1
fi

echo "Input fasta file: $input_fasta"
echo "Output directory: $output_dir"
mkdir -p $output_dir
chmod  777 $output_dir


# ./run_diSBPred_Docker.sh $(pwd)/Input/input.txt Output


ESMpath="/opt/diSBPred"
docker run --rm  -it \
	-v $input_fasta:$ESMpath/Input/FASTA \
	-v $(pwd)/Output:$ESMpath/$output_dir:rw \
	--entrypoint /bin/bash \
	wasicse/disbpred:latest 	

# Can't compile with either gfortran or ifort, aborting



