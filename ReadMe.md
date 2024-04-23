# diSBPred
diSBPred: A Machine Learning based Approach for Disulfide Bond Prediction
Authors: Avdesh Mishra, Md Wasi Ul Kabir and Md Tamjidul Hoque

*******************************************************************
diSBPred is a software used to predict disulfide bonds from protein sequence.

 ## Dataset
Under "Dataset" directory you will be able to find a list file "list_1859_Uniprot_Cutoff25_33_Intra_Cutoff25_Combined_Cutoff25_TotalProts1866.txt", which contains
the UniProt id's of the proteins used to benchmark diSBPred.

# Getting Started

## Retrieve the code

```
git clone https://github.com/wasicse/diSBPred.git

```

## Install Dependencies

We have tested diSBPred on Ubuntu 20.04. You would need to install the following software before replicating this framework in your local or server machine. 

1. pyenv latest version

2. Python version 3.7.4

3. Poetry version 1.1.13

You can install them using the following command:

```
./install_dependencies.sh
```
## Run diSBPred

To run the program, executes the following command to run diSBPred from the script directory. The tool takes input from Input/input.txt file.

```
./run_diSBPred.sh
```
## Run with Docker
- To run the diSBPred tool with docker, you can either build the docker image using dockerfile or pull the docker image from the registry.

#### Build Docker image 

```
docker build -t wasicse/disbpred https://github.com/wasicse/disbpred.git#master    
```
 #### (Alternatively) Pull image from Docker registry.

- Pull the image from the registry.
 ```
docker pull wasicse/disbpred:latest
```
#### Run run_diSBPred using Docker image
- Create the run_diSBPred container. The script will mount input fasta, and output directories from the current (run_diSBPred) directory (downloaded from GitHub) into the docker container.

```
./run_diSBPred_Docker.sh $(pwd)/Input/input.txt Output
```

## Authors

Md Wasi Ul Kabir, Md Tamjidul Hoque. For any issue please contact: Md Tamjidul Hoque, thoque@uno.edu 

## References
1. Mishra, Avdesh, et al. “diSBPred: A Machine Learning Based Approach for Disulfide Bond Prediction.” Computational Biology and Chemistry, vol. 91, Apr. 2021, p. 107436. ScienceDirect, https://doi.org/10.1016/j.compbiolchem.2021.107436.











