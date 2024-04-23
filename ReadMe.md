# diSBPred
diSBPred: A Machine Learning based Approach for Disulfide Bond Prediction
Authors: Avdesh Mishra, Md Wasi Ul Kabir and Md Tamjidul Hoque

*******************************************************************
diSBPred is a software used to predict disulfide bonds from protein sequence

 ## Dataset
Under "Data" directory you will be able to find a list file "list_1859_Uniprot_Cutoff25_33_Intra_Cutoff25_Combined_Cutoff25_TotalProts1866.txt", which contains
the UniProt id's of the proteins used to benchmark diSBPred.

- Retrieve the code

```
git clone https://github.com/wasicse/diSBPred.git

```
## INPUT
Make sure the input has a new line at the end.

### Dependencies

We have tested Dispredict3.0 on Ubuntu 20.04. You would need to install the following software before replicating this framework in your local or server machine. 

1. pyenv latest version
    ```
    curl https://pyenv.run | bash
    exec $SHELL
    ```
    For more details, visit: https://github.com/pyenv/pyenv-installer

1. Python version 3.7.4

    ```
    pyenv install miniconda3-4.7.12
    pyenv local miniconda3-4.7.12 
    ```

    For more details, visit: https://github.com/pyenv/pyenv

2. Poetry version 1.1.13

    ```
    curl -sSL https://install.python-poetry.org | python3 - --version 1.1.13
    ```
    For more details, visit: https://python-poetry.org/docs/

You can install them using the following command:

```
Install_pyenv_poetry.sh
```
## Run diSBPred

To run the program, first install all required libraries by running the following command:

```
cd diSBPred
poetry install
poetry shell
```

Then execute the following command to run Dispredict3.0 from the script directory.

```
./run_diSBPred.sh
```

## Correspondence

Please address your questions to:
	Dr. Md Tamjidul Hoque
	email: thoque@uno.edu
	Department of Computer Science
	University of New Orleans
	2000 Lakeshore Dr., New Orleans, LA, 70148












