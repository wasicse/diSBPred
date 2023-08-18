import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;


public class CreateWindowsForSingleCysteineForCombinedData {
	
	public static void main(String[] args) throws IOException {	
		
		double[][] aaFactors = new double[][]{
				
				{-0.591, -1.302, -0.733, 1.570, -0.146}, // A
				{1.538,	-0.055,	1.502, 0.440, 2.897},	// R
				{0.945, 0.828, 1.299, -0.169, 0.933},	// N
				{1.050, 0.302, -3.656, -0.259, -3.242}, // D
				{-1.343, 0.465, -0.862, -1.020, -0.255}, // C
				{0.931, -0.179, -3.005, -0.503, -1.853}, // Q
				{1.357, -1.453, 1.477, 0.113, -0.837},	// E
				{-0.384, 1.652, 1.330, 1.045, 2.064},	// G
				{0.336, -0.417, -1.673, -1.474, -0.078}, // H
				{-1.239, -0.547, 2.131, 0.393, 0.816},	// I
				{-1.019, -0.987, -1.505, 1.266, -0.912}, // L
				{1.831, -0.561, 0.533, -0.277, 1.648},	// K
				{-0.663, -1.524, 2.219, -1.005, 1.212},	// M
				{-1.006, -0.590, 1.891, -0.397, 0.412},	// F
				{0.189, 2.081, -1.628, 0.421, -1.392},	// P
				{-0.228, 1.399, -4.760, 0.670, -2.647}, // S
				{-0.032, 0.326, 2.213, 0.908, 1.313}, // T
				{-0.595, 0.009, 0.672, -2.128, -0.184}, // W
				{0.260, 0.830, 3.097, -0.838, 1.512}, // Y
				{-1.337, -0.279, -0.544, 1.242, -1.262}, // V				
				
		};
		
		
		int num_of_features = 59; // you should provide exact number of features here
		
		File currentDir = new File(new File(".").getAbsolutePath());
		String curr_dir_path = currentDir.getCanonicalPath();
		
		
		int ws= 31;

		String idListPathTest = "../Input/id_list.txt";
		String combinedFeatureDirTest = "./CombFeatures/";
		String singleCysWindowedFilesDirTest = "./SingleCysWindowedFile/";
		
		String balSSPOutputPath = "../Tools/BalancedSSP/Software/Output/prediction/";
		String dispredictOutputPath = "../Tools/DisPredict_v2.0/Software/Output/prediction/";
		String dispredictFeaturesPath = "../Tools/DisPredict_v2.0/Software/Features/";

		// turn on the commented lines below to obtain combined feature file for test dataset
		createCombinedFeaturesFiles(idListPathTest, balSSPOutputPath, dispredictOutputPath, dispredictFeaturesPath, combinedFeatureDirTest);
		
		createWindowForTest(idListPathTest, combinedFeatureDirTest, singleCysWindowedFilesDirTest, ws, num_of_features, aaFactors);
		

	}
		
	
	public static void createWindowForTest(String idListPathTest, String combinedFeatureDirTest, String singleCysWindowedFilesDirTest, int ws, int num_of_features, double[][] aaFactors) throws IOException{

		
		// parse the list file which contains protein ids
		
		BufferedReader listFileRd = BufferReaderAndWriter.getReader(new File(idListPathTest));
		String proteinID = "";
		
	
		while ((proteinID = listFileRd.readLine())!=null) {
	
			// loading file in memory
			List<String> featuresArray = new ArrayList<String>();
			
			System.out.println(proteinID);
				
			String featureFilePath = combinedFeatureDirTest+"/"+proteinID+".61features";
			
			BufferedReader featureFileReader = BufferReaderAndWriter.getReader(new File(featureFilePath));
			
			String featureLine = "";
			
			while((featureLine = featureFileReader.readLine())!=null){
			
				String[] element = featureLine.split(",");
				if(element[0].equalsIgnoreCase("disProb(1)")){
					
					continue;					
					
				}
				
				featuresArray.add(featureLine);
			
			}
			
			featureFileReader.close();
			
			String singleCysWindowedFileName = singleCysWindowedFilesDirTest+proteinID+"_singlecys_W_"+ws+".txt";
			BufferedWriter windowedFileWriter = BufferReaderAndWriter
					.getWriter(new File(singleCysWindowedFileName));
		
						
			// start creating features based on window size for each residue
			
			int beforeAndAfterRowsToBeIncluded = ws/2;
			
			
			for (int j = 0; j < featuresArray.size(); j++) {
				
				String featuresString = "";
			//	int featureCounter = 0;
				
				String[] thisIndexFeature = featuresArray.get(j).trim().split(",");
				int aaIndex = Integer.parseInt(thisIndexFeature[1]);
				
				// System.out.println(aaIndex);
				if(aaIndex != 5){
					
					continue;
					
				}
				// System.out.println("aaIndex");
	//			/*=================== Ignore the cystein pairs, either of which can not form consecutive residues equal to the window size =============*/
	//			int beforeWindowCheck = j - beforeAndAfterRowsToBeIncluded;
	//			int afterWindowCheck = j + beforeAndAfterRowsToBeIncluded;
	//			
	//			if((beforeWindowCheck < 0) || (afterWindowCheck >= featuresArray.size())){
	//				System.out.println(proteinID+"\t"+featuresArray.size()+"\t"+j+"\t"+j);
	//				continue;					
	//				
	//			}
				
				if(beforeAndAfterRowsToBeIncluded > 0){		// if window size is > 1
					
					int beforeFeatureCounter = beforeAndAfterRowsToBeIncluded;
					
					for(int k = 0; k < beforeAndAfterRowsToBeIncluded; k++){		// append the features to the left
																		
						int beforeIndex = j-beforeFeatureCounter;
						
						if(beforeIndex < 0){										// if for terminal residues features does not exist
							
							for(int l = 0; l < num_of_features; l++){
							//	featureCounter++;
//								featuresString = featuresString+"9999.9999,";
								featuresString = featuresString+"0,";
								
							}
							
						}else {
							
							String[] beforeFeatureExist = featuresArray.get(beforeIndex).trim().split(",");
							
							for(int m = 0; m < beforeFeatureExist.length; m++){
							//	featureCounter++;
							// 	featuresString = featuresString+featureCounter+":"+beforeFeatureExist[m]+" ";
							//	featuresString = featuresString+beforeFeatureExist[m]+" ";
	//							String[] feature = beforeFeatureExist[m].split(":");
								
								if((m >= 2 && m <= 8)){	 // ignore the physical parameter features
									
									continue;
									
								}
								
								if(m == 1){	// aaFactors features
									
									featuresString = featuresString+beforeFeatureExist[m]+",";
									int aaInd =  Integer.parseInt(beforeFeatureExist[m]);
																		
									double[] aaFact = {0, 0, 0, 0, 0};
									
									if(!(aaInd == 0)){
										
										aaFact =  aaFactors[aaInd-1];
										
									}
									
									for(int q = 0; q < 5; q++){
																				
										featuresString = featuresString+aaFact[q]+",";
										
										
									}
									
									continue;
									
									
								}
								
								featuresString = featuresString+beforeFeatureExist[m]+",";
								
								
							}
							
							
						}
						
						beforeFeatureCounter--;
												
					}
					
					String[] currentIndexFeature = featuresArray.get(j).trim().split(",");
					
					for(int n = 0; n < currentIndexFeature.length; n++){
					//	featureCounter++;
					//	featuresString = featuresString+featureCounter+":"+currentIndexFeature[n]+" ";
					//	featuresString = featuresString+currentIndexFeature[n]+" ";
						if(n == 1){	// skip aa index for cysteine residues
							
							continue;
						}
						
						if((n >= 2 && n <= 8)){	 // ignore the physical parameter features
							
							continue;
							
						}
						
	//					String[] feature = currentIndexFeature[n].split(":");
						featuresString = featuresString+currentIndexFeature[n]+",";
						
					}
					
	//				featuresString = currentIndexFeature[0]+","+featuresString;
					
										
					int afterFeatureCounter = 0;
					
					for(int k = 0; k < beforeAndAfterRowsToBeIncluded; k++){		// append the features to the right
						
						afterFeatureCounter++;
						
						int afterIndex = j+afterFeatureCounter;
						
						if(afterIndex > featuresArray.size()-1){
							
							for(int l = 0; l < num_of_features; l++){
							//	featureCounter++;
							//	featuresString = featuresString+featureCounter+":0 ";
								
//								featuresString = featuresString+"9999.9999,";
								featuresString = featuresString+"0,";
								
							}
							
						}else {
							
							String[] afterFeatureExist = featuresArray.get(afterIndex).trim().split(",");
							
							for(int m = 0; m < afterFeatureExist.length; m++){
							//	featureCounter++;
							//	featuresString = featuresString+featureCounter+":"+afterFeatureExist[m]+" ";
							//	featuresString = featuresString+afterFeatureExist[m]+" ";
	//							String[] feature = afterFeatureExist[m].split(":");
								
								if((m >= 2 && m <= 8)){	 // ignore the physical parameter features
									
									continue;
									
								}
								
								if(m == 1){	// aaFactors features
									
									featuresString = featuresString+afterFeatureExist[m]+",";
									int aaInd =  Integer.parseInt(afterFeatureExist[m]);
																		
									double[] aaFact = {0, 0, 0, 0, 0};
									
									if(!(aaInd == 0)){
										
										aaFact =  aaFactors[aaInd-1];
										
									}
									
									for(int q = 0; q < 5; q++){
																				
										featuresString = featuresString+aaFact[q]+",";										
										
									}
									
									continue;
									
									
								}
								
								featuresString = featuresString+afterFeatureExist[m]+",";
								
							}
							
																					
						}
												
					}
					
					
				}else { // if window size is 1
					
					String[] thisFeature = featuresArray.get(j).trim().split(",");
					
					for(int p = 0; p < thisFeature.length; p++){
						
					//	featureCounter++;
					//	featuresString = featuresString+featureCounter+":"+thisFeature[p]+" ";
					//	featuresString = featuresString+thisFeature[p]+" ";
	//					String[] feature = thisFeature[p].split(":");
						
						if(p == 1){
							
							continue;
						}
						
						if((p >= 2 && p <= 8)){	 // ignore the physical parameter features
							
							continue;
							
						}
						
						featuresString = featuresString+thisFeature[p]+",";							
						
					}
					
	//				featuresString = thisFeature[0]+","+featuresString.trim();
					
										
				}
				
									
				featuresString = "0"+","+featuresString.trim();
					
				
				
				if(featuresString.contains("inf") || featuresString.contains("nan") ||featuresString.contains("null")){
					
					System.out.println("protein contains inf "+proteinID);
					System.exit(0);
				}
				windowedFileWriter.write(featuresString.substring(0, featuresString.length()-1));
//				windowedFileWriter.write(featuresString);
				windowedFileWriter.write("\n");
				windowedFileWriter.flush();
				
			}
			
			windowedFileWriter.close();
			
		}
		
		listFileRd.close();
		
		
		
	}
	
	
	public static void createCombinedFeaturesFiles(String idListPath, String balSSPOutputPath, String dispredictOutputPath, String dispredictFeaturesPath, String outputFeatureDir) throws IOException{

		// combine the ".57features" file with the disorder probability (Dispredict_V2.0) and the secondary structure probability from BalancedSSP
		BufferedReader idListRd = BufferReaderAndWriter.getReader(new File(idListPath));
		String id = "";
						
					
		while((id = idListRd.readLine())!=null){
			
//			System.out.println(id);
			
			String disPredIpFeaPath = dispredictFeaturesPath+id+"/"+id+".57pfeatures";
			BufferedReader disPredFeatureRd = BufferReaderAndWriter.getReader(new File(disPredIpFeaPath));
			
			
			String disPredProbPath = dispredictOutputPath+id+"/DisPredict2/"+id+".drp";
			BufferedReader disPredOutputRd = BufferReaderAndWriter.getReader(new File(disPredProbPath));
			String disPredOutputLine = "";
			
			String balSSPProbPath = balSSPOutputPath+id+"/MetaSSPred/"+id+".MetaSSpred";
			BufferedReader balSSPOutputRd = BufferReaderAndWriter.getReader(new File(balSSPProbPath));
			
			
			String outputFile = outputFeatureDir+id+".61features";
			BufferedWriter opWr = BufferReaderAndWriter.getWriter(new File(outputFile));
			opWr.write("disProb(1),");
			
			
			if((disPredFeatureRd!=null) && (disPredOutputRd !=null) && (balSSPOutputRd != null)){
				
				for(int i = 0; i < 6; i++){ // skip the header information of disorder probability file
					
					disPredOutputRd.readLine();
					
				}
				
				String featuresHeader = disPredFeatureRd.readLine(); // skip header lines from .57pfeatures file
				opWr.write(featuresHeader.substring(7, featuresHeader.length())); 
				
				for(int i = 0; i < 4; i++){ // skip header lines from .MetaSSPred file
					
					balSSPOutputRd.readLine();
					
				}
				
				
			}
			
			opWr.write(",pBeta,pCoil,pHelix\n");
			
			while((disPredOutputLine = disPredOutputRd.readLine())!=null){
				
				if(disPredOutputLine.equalsIgnoreCase("END")){
					
					continue;
				}
				
				String[] disPredElements = disPredOutputLine.split("\\s+");
				String disOrdProb = disPredElements[2];
				
				String[] feature57 = disPredFeatureRd.readLine().split("\\s+");
				String subFeature = "";
				for(int i = 1; i < feature57.length; i++){
					
					subFeature +=feature57[i]+",";
					
				}
				
				String[] metaSSFeatures = balSSPOutputRd.readLine().split("\\s+");
				String ssFeature = metaSSFeatures[metaSSFeatures.length-3]+","+metaSSFeatures[metaSSFeatures.length-2]+","+metaSSFeatures[metaSSFeatures.length-1]+"\n";
				
				String combinedFeatures = disOrdProb+","+subFeature+ssFeature;
				
				if(combinedFeatures.contains("nan") || combinedFeatures.contains("null") || combinedFeatures.contains("inf")){
					
					System.out.println("contains nan or null "+id);
					break;
					
				}
				
				opWr.write(combinedFeatures);				
				
				
				
			}
			
			opWr.flush();
			opWr.close();
			disPredOutputRd.close();
			disPredFeatureRd.close();
			balSSPOutputRd.close();	
			
			
			
		}
		
		idListRd.close();
		
		
		
	}
	

}
