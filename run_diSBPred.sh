#! /bin/bash

./download_largefiles.sh

#Install dependencies
source ./.venv/bin/activate
./.venv/bin/poetry install

Home=$PWD

echo "Removing old Results"
rm -r Scripts/CombFeatures
rm -r SingleCysWindowedFile
rm -r Scripts/IndividualCysPredictionProb
rm -r Scripts/CysPairWindowedFiles_IncCysProb
rm -r Scripts/Predictions_IncSingCysProb_5folddata


rm -r Tools/DisPredict_v2.0/Software/PSSM/*
rm -r Tools/DisPredict_v2.0/Software/Features/*
rm -r Tools/DisPredict_v2.0/Software/Output/log/*
rm -r Tools/DisPredict_v2.0/Software/Output/prediction/*

rm -r Tools/DisPredict_v2.0/Software/Input
rm -r Tools/BalancedSSP/Software/PSSM/*
rm -r Tools/BalancedSSP/Software/Features/*
rm -r Tools/BalancedSSP/Software/Output/log/*
rm -r Tools/BalancedSSP/Software/Output/prediction/*

rm -r Tools/BalancedSSP/Software/Input

mkdir Scripts/CombFeatures
mkdir Scripts/SingleCysWindowedFile
mkdir Scripts/IndividualCysPredictionProb
mkdir Scripts/CysPairWindowedFiles_IncCysProb
mkdir Scripts/Predictions_IncSingCysProb_5folddata


echo "Running DisPredict_v2.0"
cp -r ./Input Tools/DisPredict_v2.0/Software/
cd Tools/DisPredict_v2.0/Software/Scripts
./run_DisPredict_v2.0 2>&1 | tee output_DisPredict.txt

cd $Home
echo "Running BalancedSSP"
cp -r ./Input Tools/BalancedSSP/Software/
cd Tools/BalancedSSP/Software/Scripts
./run_BalancedSSP 2>&1 | tee output_BalancedSSP.txt

cd $Home
cd Scripts

echo "Collecting Features and Start Predictions"
javac *.java 
java CreateWindowsForSingleCysteineForCombinedData
python LightGBM_CollectIndividualCysProb_ws31.py
java CreateWindowsForCysteinePairForCombinedData
python Stacking_CysPairLoadModel_Test_ws13.py