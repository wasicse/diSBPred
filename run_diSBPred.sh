#! /bin/bash

# parse input
mkdir -p ../Input/FASTA
cd Scripts
../.venv/bin/python parseFasta.py 
cd -

# # Download large files
# ./download_largefiles.sh
# ./download_dataset.sh
# #Install dependencies
# source ./.venv/bin/activate
# ./.venv/bin/poetry install

# Home=$PWD

# # ask for confirmation before removing the old results
# read -p "Removing old files? " -n 1 -r
# echo    # (optional) move to a new line
# if [[ $REPLY =~ ^[Yy]$ ]]
# then
    
#     echo "Removing old Results"
#     rm -f -r Scripts/CombFeatures
#     rm -f -r SingleCysWindowedFile
#     rm -f -r Scripts/IndividualCysPredictionProb
#     rm -f -r Scripts/CysPairWindowedFiles_IncCysProb
#     rm -f -r Scripts/Predictions_IncSingCysProb_5folddata


#     rm -f -r Tools/DisPredict_v2.0/Software/PSSM/*
#     rm -f -r Tools/DisPredict_v2.0/Software/Features/*
#     rm -f -r Tools/DisPredict_v2.0/Software/Output/log/*
#     rm -f -r Tools/DisPredict_v2.0/Software/Output/prediction/*

#     rm -f -r Tools/DisPredict_v2.0/Software/Input
#     rm -f -r Tools/BalancedSSP/Software/PSSM/*
#     rm -f -r Tools/BalancedSSP/Software/Features/*
#     rm -f -r Tools/BalancedSSP/Software/Output/log/*
#     rm -f -r Tools/BalancedSSP/Software/Output/prediction/*

#     rm -f -r Tools/BalancedSSP/Software/Input
#     echo "Old results are not removed"
#     # exit 1
# fi


# mkdir -p Scripts/CombFeatures
# mkdir -p Scripts/SingleCysWindowedFile
# mkdir -p Scripts/IndividualCysPredictionProb
# mkdir -p Scripts/CysPairWindowedFiles_IncCysProb
# mkdir -p Scripts/Predictions_IncSingCysProb_5folddata


# echo "Running DisPredict_v2.0"
# cp -r ./Input Tools/DisPredict_v2.0/Software/
# cd Tools/DisPredict_v2.0/Software/Scripts
# ./run_DisPredict_v2.0 2>&1 | tee output_DisPredict.txt

# cd $Home
# echo "Running BalancedSSP"
# cp -r ./Input Tools/BalancedSSP/Software/
# cd Tools/BalancedSSP/Software/Scripts
# ./run_BalancedSSP 2>&1 | tee output_BalancedSSP.txt

# cd $Home
# cd Scripts

# echo "Collecting Features and Start Predictions"
# javac *.java 
# java CreateWindowsForSingleCysteineForCombinedData
# python LightGBM_CollectIndividualCysProb_ws31.py
# java CreateWindowsForCysteinePairForCombinedData
# python Stacking_CysPairLoadModel_Test_ws13.py