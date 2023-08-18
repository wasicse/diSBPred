import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;


public class CreateWindowsForCysteinePairForCombinedData {
	
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
		
				
		int num_of_features = 61; // you should provide exact number of features here
		
		File currentDir = new File(new File(".").getAbsolutePath());
		String curr_dir_path = currentDir.getCanonicalPath();
		
		int ws = 13;
		String singleCysPredProbPathTest = "./IndividualCysPredictionProb/";
		String idListPathTest = "../Input/id_list.txt";
		String combinedFeatureDirTest = "./CombFeatures/";
		String cysPairWindowedFilesDirTest = "./CysPairWindowedFiles_IncCysProb/";
		String fastaFileDirTest = "../Input/FASTA/";
		createWindowForCysPairForTestData(idListPathTest, fastaFileDirTest, cysPairWindowedFilesDirTest, combinedFeatureDirTest, singleCysPredProbPathTest, ws, num_of_features, aaFactors);
		






	}
	
	
	public static void createWindowForCysPair(String idListPath, String cysPairWindowedFileName, String combinedFeatureDir, String singleCysPredProbPath, String cysConnectionPath, int window_size, int num_of_features, double[][] aaFactors) throws IOException{
		
		BufferedWriter winWr = BufferReaderAndWriter.getWriter(new File(cysPairWindowedFileName));		
//		winWr.write("class,");
//		
//		for(int i = 0; i < 1051; i++){
//			
//			if(i == 1050){
//				
//				winWr.write("V"+(i+1));
//				
//			}else{
//				
//				winWr.write("V"+(i+1)+",");
//				
//			}			
//			
//		}
//		
//		winWr.write("\n");
		
		List<String> formatted_features = new ArrayList<String>();
		List<String> positive_features = new ArrayList<String>();
		List<String> negative_features = new ArrayList<String>();
		
		List<Integer> positiveCysPairDistance = new ArrayList<Integer>();
		List<Integer> negativeCysPairDistance = new ArrayList<Integer>();
		
		double positiveClassCount = 0;
		double negativeClassCount = 0;
		double positiveCysPairsCount = 0;
		double negativeCysPairsCount = 0;
		
		
		BufferedReader br_id = BufferReaderAndWriter.getReader(new File(idListPath));
		String pdbId = "";
		
		while((pdbId = br_id.readLine())!=null){
			
			System.out.println("pdbId "+pdbId);
			// parse the connectivity file to obtain cys pairs and the total number of cysteins present in the structure
			BufferedReader brConn = BufferReaderAndWriter.getReader(new File(cysConnectionPath+pdbId+".ssbond"));
			String connLine = "";
			
			List<String> cysPairs = new ArrayList<String>(); // contains the disulfide cysteine pairs
			int totalCys = -1;  // total number of cys present in the structure
			
			while((connLine = brConn.readLine())!=null){
				
				if(connLine.startsWith("RltvPos")){
					
					String[] elements = connLine.split("\\s+");
					cysPairs.add(elements[1]);
									
				}else if(connLine.startsWith("TotNumCys")){
					
					String[] elements = connLine.split("\\s+");
					totalCys = Integer.parseInt(elements[1]);
					
				}
				
			}
			
			brConn.close();
			
			/* ================== create total possible combination of cysteine pairs for the structure ======================== */

			List<String> cysCombs = new ArrayList<String>();	// holds all possible combination of the CYS in a protein
			for(int i = 1; i < totalCys; i++){
				
				for(int j = i+1; j <= totalCys; j++){
					
					cysCombs.add("c"+i+"-"+"c"+j);
					
					
				}
				
				
			}
			
			
/* ================== parse the single cysteine predict file and collect cysteine bonding prob to use it as a feature ======================== */
			
//			Map<String, Double> cysBondProb = new HashMap<String, Double>();
//			
//			String cysBondPredFilePath = singleCysPredProbPath+pdbId+".predict";
//			BufferedReader cysBondPredFileRd = BufferReaderAndWriter.getReader(new File(cysBondPredFilePath));
//			String cysBondLine = "";
//			int cysId = 0;
////			int classOneCount = 0;
////			System.out.println(pdbId);
//			while((cysBondLine = cysBondPredFileRd.readLine())!=null){
//				
//				if(cysBondLine.startsWith("predClass")){
//					
//					continue;
//					
//				}
//				cysId++;
//				String[] values = cysBondLine.split("\\s+");
//				cysBondProb.put("c"+cysId, Double.parseDouble(values[2]));
//				
////				if(Integer.parseInt(values[0]) == 1){
////					
////					classOneCount++;
////					
////				}
//				
//				
//			}
//						
//			cysBondPredFileRd.close();
			
			
/* ================== parse the features file ======================== */
			
			// loading file in memory
			List<String> featuresArray = new ArrayList<String>();	// holds features for each residues
			Map<String, Integer> cysInd = new HashMap<String, Integer>();	 	// holds the cysteine relative position and sequence index e.g. key = c1 and Val = 10
						
			String featureFilePath = combinedFeatureDir+pdbId+".61features";
			
			BufferedReader featureFileReader = BufferReaderAndWriter.getReader(new File(featureFilePath));
			
			String featureLine = "";
			int cysRelPos = 0;
			
			while((featureLine = featureFileReader.readLine())!=null){
			
				
				if(featureLine.startsWith("disProb(1)")){
					
					continue;
				}
				
				featuresArray.add(featureLine);
				
				String[] featureElements = featureLine.split(",");
				String aa_type = featureElements[1];
				if(Integer.parseInt(aa_type) == 5){
					
					cysRelPos++;					
					cysInd.put("c"+cysRelPos, featuresArray.size()-1); // sequence position starts from 0
										
				}			
			
			}
			
			featureFileReader.close(); 
			
/* ================== Collect Windowed Feature for all possible combination of cys pair and assign the class label ======================== */
			
			for(int i = 0; i < cysCombs.size(); i++){
				
				String cysComb = cysCombs.get(i);
				int classType = 0;
				
				if(cysPairs.contains(cysComb)){
					
					classType = +1;
					positiveClassCount += 1;
					
				}else{
					
					classType = 0;
					negativeClassCount += 1;
					
				}			
				
				
				String[] cysCombArray = cysComb.split("-");
				String cysCombFirst = cysCombArray[0];
				String cysCombSec = cysCombArray[1];
				
//				double firstCysBondProb = cysBondProb.get(cysCombFirst);
//				double secCysBondProb = cysBondProb.get(cysCombSec);
				
				int cysFirstSeqInd = cysInd.get(cysCombFirst);
				int cysSecSeqInd = cysInd.get(cysCombSec);
				
				int cysDistance = Math.abs(cysSecSeqInd - cysFirstSeqInd);
				
				if(classType == 1 && cysDistance < 100){
					
//					positiveCysPairDistance.add(cysDistance);
					positiveCysPairsCount += 1; 
					
				}else if(classType == 0 && cysDistance < 100){
					
//					negativeCysPairDistance.add(cysDistance);
					negativeCysPairsCount += 1;
					
				}
				
				
				
				if(classType == 0 && cysDistance >= 100){
					
					continue;
					
				}
				
				// for 5-fold negative sample and positive samples turn this condition on
//				if(classType == 0 && negativeCysPairsCount > 43935){
//					
//					continue;
//					
//				}
				// for new dataset
//				if(classType == 0 && negativeCysPairsCount > 40280){
//					
//					continue;
//					
//				}
				
				// for balanced negative and positive samples turn this condition on
//				if(classType == 0 && negativeCysPairsCount > 8787){
//					
//					continue;
//					
//				}
				// for new dataset
				if(classType == 0 && negativeCysPairsCount > 8056){
				
					continue;
				
				}
				
				
				StringBuilder sb = new StringBuilder();
				StringBuilder sb_before = new StringBuilder();
				StringBuilder sb_after = new StringBuilder();
				StringBuilder sb_cys_site = new StringBuilder();
				
				sb.append(Integer.toString(classType)+",");
				
				int beforeAndAfterRowsToBeIncluded = window_size/2;
				
//				/*=================== Ignore the cystein pairs, either of which can not form consecutive residues equal to the window size =============*/
//				int beforeWindowCheck = cysFirstSeqInd - beforeAndAfterRowsToBeIncluded;
//				int afterWindowCheck = cysSecSeqInd + beforeAndAfterRowsToBeIncluded;
//				
//				if((beforeWindowCheck < 0) || (afterWindowCheck >= featuresArray.size())){
//					System.out.println(pdbId+"\t"+featuresArray.size()+"\t"+cysFirstSeqInd+"\t"+cysSecSeqInd);
//					continue;					
//					
//				}
				
				if(beforeAndAfterRowsToBeIncluded > 0){		// if window size is > 1
					
					int beforeFeatureCounter = beforeAndAfterRowsToBeIncluded;
					int afterFeatureCounter = 0;
					
					List<Integer> aaIndexVectorFirst = new ArrayList<Integer>();
					List<Integer> aaIndexVectorSec = new ArrayList<Integer>();
					
					for(int k = 0; k < beforeAndAfterRowsToBeIncluded; k++){		// append the features to the left
						afterFeatureCounter++;												
						int firstBeforeIndex = cysFirstSeqInd-beforeFeatureCounter;
						int secBeforeIndex = cysSecSeqInd-beforeFeatureCounter;
						int firstAfterIndex = cysFirstSeqInd+afterFeatureCounter;
						int secAfterIndex = cysSecSeqInd+afterFeatureCounter;
						List<Double> firstBeforeFeatures = new ArrayList<Double>();
						List<Double> secBeforeFeatures = new ArrayList<Double>();
						List<Double> firstAfterFeatures = new ArrayList<Double>();
						List<Double> secAfterFeatures = new ArrayList<Double>();
												
												
						if(firstBeforeIndex < 0){										// if for terminal residues features does not exist
							
							for(int l = 0; l < num_of_features; l++){
							//	featureCounter++;
								firstBeforeFeatures.add(0.0);
								
							}
							aaIndexVectorFirst.add(0);
							
						}else {
							
							String[] beforeFeatureExist = featuresArray.get(firstBeforeIndex).trim().split(",");
							
							for(int m = 0; m < beforeFeatureExist.length; m++){
								
								if(m == 1){ // ignore AA based feature
																		
									aaIndexVectorFirst.add(Integer.parseInt(beforeFeatureExist[m]));
								//	continue;
								}
							//	featureCounter++;
							// 	featuresString = featuresString+featureCounter+":"+beforeFeatureExist[m]+" ";
							//	featuresString = featuresString+beforeFeatureExist[m]+" ";
								
								firstBeforeFeatures.add(Double.parseDouble(beforeFeatureExist[m]));
								
							}							
							
						}
						
						if(secBeforeIndex < 0){
							
							for(int l = 0; l < num_of_features; l++){
								//	featureCounter++;
									secBeforeFeatures.add(0.0);
									
							}
							
							aaIndexVectorSec.add(0);
							
							
						}else{
							
							String[] beforeFeatureExist = featuresArray.get(secBeforeIndex).trim().split(",");
							
							for(int m = 0; m < beforeFeatureExist.length; m++){
								if(m == 1){	// ignore aa based feature
																		
									aaIndexVectorSec.add(Integer.parseInt(beforeFeatureExist[m]));
//									continue;
								}
							//	featureCounter++;
							// 	featuresString = featuresString+featureCounter+":"+beforeFeatureExist[m]+" ";
							//	featuresString = featuresString+beforeFeatureExist[m]+" ";
								
								secBeforeFeatures.add(Double.parseDouble(beforeFeatureExist[m]));
								
							}
							
						}
						
						if(firstAfterIndex > featuresArray.size()-1){										// if for terminal residues features does not exist
							
							for(int l = 0; l < num_of_features; l++){
							//	featureCounter++;
								firstAfterFeatures.add(0.0);
								
							}
							
							aaIndexVectorFirst.add(0);
							
						}else {
							
							String[] afterFeatureExist = featuresArray.get(firstAfterIndex).trim().split(",");
							
							for(int m = 0; m < afterFeatureExist.length; m++){
								if(m == 1){	// ignore aa based feature
																		
									aaIndexVectorFirst.add(Integer.parseInt(afterFeatureExist[m]));
//									continue;
								}
							//	featureCounter++;
							// 	featuresString = featuresString+featureCounter+":"+beforeFeatureExist[m]+" ";
							//	featuresString = featuresString+beforeFeatureExist[m]+" ";
								
								firstAfterFeatures.add(Double.parseDouble(afterFeatureExist[m]));
								
							}
							
							
						}
						
						if(secAfterIndex > featuresArray.size()-1){
							
							for(int l = 0; l < num_of_features; l++){
								//	featureCounter++;
									secAfterFeatures.add(0.0);
									
							}
							
							aaIndexVectorSec.add(0);
							
							
						}else{
							
							String[] afterFeatureExist = featuresArray.get(secAfterIndex).trim().split(",");
							
							for(int m = 0; m < afterFeatureExist.length; m++){
								if(m == 1){	// ignore aa based feature
																		
									aaIndexVectorSec.add(Integer.parseInt(afterFeatureExist[m]));
//									continue;
								}
							//	featureCounter++;
							// 	featuresString = featuresString+featureCounter+":"+beforeFeatureExist[m]+" ";
							//	featuresString = featuresString+beforeFeatureExist[m]+" ";
								
								secAfterFeatures.add(Double.parseDouble(afterFeatureExist[m]));
								
							}
							
						}
						
						for(int p = 0; p < firstBeforeFeatures.size(); p++){ // firstBeforeFeatures contains all the features except the class label information
							
//							if((p >= 2 && p <= 8) || p > 28){
//								
//								continue;
//								
//							}
							
							if((p >= 2 && p <= 8)){	 // ignore the physical parameter features
								
//								double diffFeature = Math.abs((firstBeforeFeatures.get(p) - secBeforeFeatures.get(p)));
//								double sumFeature = Math.abs((firstBeforeFeatures.get(p) + secBeforeFeatures.get(p)));
//								sb_before.append(diffFeature+","+sumFeature+",");
								continue;
								
							}
							
							if(p == 0){	// disorder feature
								
								double diffFeature = Math.abs((firstBeforeFeatures.get(p) - secBeforeFeatures.get(p)));
								double sumFeature = Math.abs((firstBeforeFeatures.get(p) + secBeforeFeatures.get(p)));
								sb_before.append(diffFeature+","+sumFeature+",");
								
//								double diffFeature = Math.abs((firstBeforeFeatures.get(p) - secBeforeFeatures.get(p))*(firstBeforeFeatures.get(p) - secBeforeFeatures.get(p)));
//								double sumFeature = Math.abs((firstBeforeFeatures.get(p) + secBeforeFeatures.get(p))*(firstBeforeFeatures.get(p) + secBeforeFeatures.get(p)));
//								sb_before.append(diffFeature+","+sumFeature+",");
								
							}
							
							if(p == 1){	// aaFactors features
								
								int firstAAIndex =  firstBeforeFeatures.get(p).intValue();
								int secAAIndex = secBeforeFeatures.get(p).intValue();
								
								double[] firstAAFactors = {0, 0, 0, 0, 0};
								double[] secAAFactors = {0, 0, 0, 0, 0};
								
								if(!(firstAAIndex == 0)){
									
									firstAAFactors =  aaFactors[firstAAIndex-1];
									
								}
								
								if(!(secAAIndex == 0)){
									
									secAAFactors = aaFactors[secAAIndex-1];
									
								}
								
								for(int q = 0; q < 5; q++){
									
									double diffFeature = Math.abs(firstAAFactors[q] - secAAFactors[q]);
									double sumFeature = Math.abs(firstAAFactors[q] + secAAFactors[q]);
									sb_before.append(diffFeature+","+sumFeature+",");
									
//									double diffFeature = Math.abs((firstAAFactors[q] - secAAFactors[q])*(firstAAFactors[q] - secAAFactors[q]));
//									double sumFeature = Math.abs((firstAAFactors[q] + secAAFactors[q])*(firstAAFactors[q] + secAAFactors[q]));
//									sb_before.append(diffFeature+","+sumFeature+",");
									
									
								}				
								
								
							}						
														
							
							if(p >= 9 && p <= 28){	// PSSM features
								
								double diffFeature = Math.abs((firstBeforeFeatures.get(p) - secBeforeFeatures.get(p)));
								double sumFeature = Math.abs((firstBeforeFeatures.get(p) + secBeforeFeatures.get(p)));
								sb_before.append(diffFeature+","+sumFeature+",");
								
//								double diffFeature = Math.abs((firstBeforeFeatures.get(p) - secBeforeFeatures.get(p))*(firstBeforeFeatures.get(p) - secBeforeFeatures.get(p)));
//								double sumFeature = Math.abs((firstBeforeFeatures.get(p) + secBeforeFeatures.get(p))*(firstBeforeFeatures.get(p) + secBeforeFeatures.get(p)));
//								sb_before.append(diffFeature+","+sumFeature+",");
															
								
							}
							
							if(p > 28){ // all other features
								
								double diffFeature = Math.abs((firstBeforeFeatures.get(p) - secBeforeFeatures.get(p)));
								double sumFeature = Math.abs((firstBeforeFeatures.get(p) + secBeforeFeatures.get(p)));
								sb_before.append(diffFeature+","+sumFeature+",");
								
//								double diffFeature = Math.abs((firstBeforeFeatures.get(p) - secBeforeFeatures.get(p))*(firstBeforeFeatures.get(p) - secBeforeFeatures.get(p)));
//								double sumFeature = Math.abs((firstBeforeFeatures.get(p) + secBeforeFeatures.get(p))*(firstBeforeFeatures.get(p) + secBeforeFeatures.get(p)));
//								sb_before.append(diffFeature+","+sumFeature+",");
								
							}
							
							
						}
						
						for(int p = 0; p < firstAfterFeatures.size(); p++){
							
//							if((p >= 2 && p <= 8) || p > 28){
//								
//								continue;
//								
//							}
							
							if((p >= 2 && p <= 8)){	// ignore the physical parameter feature
							
//								double diffFeature = Math.abs((firstAfterFeatures.get(p) - secAfterFeatures.get(p)));
//								double sumFeature = Math.abs((firstAfterFeatures.get(p) + secAfterFeatures.get(p)));
//								sb_after.append(diffFeature+","+sumFeature+",");
								continue;
								
							}
							
							if(p == 0){	// disorder feature
								
								double diffFeature = Math.abs((firstAfterFeatures.get(p) - secAfterFeatures.get(p)));
								double sumFeature = Math.abs((firstAfterFeatures.get(p) + secAfterFeatures.get(p)));
								sb_after.append(diffFeature+","+sumFeature+",");
								
//								double diffFeature = Math.abs((firstAfterFeatures.get(p) - secAfterFeatures.get(p))*(firstAfterFeatures.get(p) - secAfterFeatures.get(p)));
//								double sumFeature = Math.abs((firstAfterFeatures.get(p) + secAfterFeatures.get(p))*(firstAfterFeatures.get(p) + secAfterFeatures.get(p)));
//								sb_after.append(diffFeature+","+sumFeature+",");
								
							}
							
							if(p == 1){	// aaFactors features
								
								int firstAAIndex =  firstAfterFeatures.get(p).intValue();
								int secAAIndex = secAfterFeatures.get(p).intValue();
								
								double[] firstAAFactors = {0, 0, 0, 0, 0};
								double[] secAAFactors = {0, 0, 0, 0, 0};
								if(!(firstAAIndex == 0)){
									
									firstAAFactors =  aaFactors[firstAAIndex-1];
									
								}
								
								if(!(secAAIndex == 0)){
									
									secAAFactors = aaFactors[secAAIndex-1];
									
								}								
								
								
								for(int q = 0; q < 5; q++){
									
									double diffFeature = Math.abs(firstAAFactors[q] - secAAFactors[q]);
									double sumFeature = Math.abs(firstAAFactors[q] + secAAFactors[q]);
									sb_after.append(diffFeature+","+sumFeature+",");
									
//									double diffFeature = Math.abs((firstAAFactors[q] - secAAFactors[q])*(firstAAFactors[q] - secAAFactors[q]));
//									double sumFeature = Math.abs((firstAAFactors[q] + secAAFactors[q])*(firstAAFactors[q] + secAAFactors[q]));
//									sb_after.append(diffFeature+","+sumFeature+",");
									
									
								}				
								
								
							}						
														
							
							if(p >= 9 && p <= 28){	// PSSM features
								
								double diffFeature = Math.abs((firstAfterFeatures.get(p) - secAfterFeatures.get(p)));
								double sumFeature = Math.abs((firstAfterFeatures.get(p) + secAfterFeatures.get(p)));
								sb_after.append(diffFeature+","+sumFeature+",");
								
//								double diffFeature = Math.abs((firstAfterFeatures.get(p) - secAfterFeatures.get(p))*(firstAfterFeatures.get(p) - secAfterFeatures.get(p)));
//								double sumFeature = Math.abs((firstAfterFeatures.get(p) + secAfterFeatures.get(p))*(firstAfterFeatures.get(p) + secAfterFeatures.get(p)));
//								sb_after.append(diffFeature+","+sumFeature+",");
																
								
							}
							
							if(p > 28){ // all other features
								
								double diffFeature = Math.abs((firstAfterFeatures.get(p) - secAfterFeatures.get(p)));
								double sumFeature = Math.abs((firstAfterFeatures.get(p) + secAfterFeatures.get(p)));
								sb_after.append(diffFeature+","+sumFeature+",");
								
//								double diffFeature = Math.abs((firstAfterFeatures.get(p) - secAfterFeatures.get(p))*(firstAfterFeatures.get(p) - secAfterFeatures.get(p)));
//								double sumFeature = Math.abs((firstAfterFeatures.get(p) + secAfterFeatures.get(p))*(firstAfterFeatures.get(p) + secAfterFeatures.get(p)));
//								sb_after.append(diffFeature+","+sumFeature+",");
								
							}
							
							
						}
						
												
						beforeFeatureCounter--;
												
					}
					
					String[] firstFeature = featuresArray.get(cysFirstSeqInd).trim().split(",");
					String[] secFeature = featuresArray.get(cysSecSeqInd).trim().split(",");
					List<Double> firstCysFeatures = new ArrayList<Double>();
					List<Double> secCysFeatures = new ArrayList<Double>();
					
					for(int p = 0; p < firstFeature.length; p++){
											
						
//						if(p == 35 || p == 36){	// ignore dpsi and MG features
//							
//							continue;
//						}
						
						// sb.append(feature[1]+",");
						firstCysFeatures.add(Double.parseDouble(firstFeature[p]));
						
					}
					
					for(int p = 0; p < secFeature.length; p++){
									
						
//						if(p == 35 || p == 36){	// ignore dpsi and MG features
//							
//							continue;
//						}
									
						// sb.append(feature[1]+",");
						secCysFeatures.add(Double.parseDouble(secFeature[p]));
						
					}
					
					for(int p = 0; p < firstCysFeatures.size(); p++){
						
						
//						if(p > 20){
//							
//							continue;
//							
//						}
						
						if(p >= 1 && p <= 8){ // ignore the residue index and physical properties of the cysteines as feature because they will be same for both the cysteines
							
							continue;
						}
						
						if(p == 0){	// disorder feature
							
							double diffFeature = Math.abs((firstCysFeatures.get(p) - secCysFeatures.get(p)));
							double sumFeature = Math.abs((firstCysFeatures.get(p) + secCysFeatures.get(p)));
							sb_cys_site.append(diffFeature+","+sumFeature+",");
							
//							double diffFeature = Math.abs((firstCysFeatures.get(p) - secCysFeatures.get(p))*(firstCysFeatures.get(p) - secCysFeatures.get(p)));
//							double sumFeature = Math.abs((firstCysFeatures.get(p) + secCysFeatures.get(p))*(firstCysFeatures.get(p) + secCysFeatures.get(p)));
//							sb_cys_site.append(diffFeature+","+sumFeature+",");
							
						}
						
						if(p >= 9 && p <= 28){	// PSSM features
							
							double diffFeature = Math.abs((firstCysFeatures.get(p) - secCysFeatures.get(p)));
							double sumFeature = Math.abs((firstCysFeatures.get(p) + secCysFeatures.get(p)));
							sb_cys_site.append(diffFeature+","+sumFeature+",");
							
//							double diffFeature = Math.abs((firstCysFeatures.get(p) - secCysFeatures.get(p))*(firstCysFeatures.get(p) - secCysFeatures.get(p)));
//							double sumFeature = Math.abs((firstCysFeatures.get(p) + secCysFeatures.get(p))*(firstCysFeatures.get(p) + secCysFeatures.get(p)));
//							sb_cys_site.append(diffFeature+","+sumFeature+",");
												
							
						}
						
						if(p > 28){	// consider all other features too
							
							double diffFeature = Math.abs((firstCysFeatures.get(p) - secCysFeatures.get(p)));
							double sumFeature = Math.abs((firstCysFeatures.get(p) + secCysFeatures.get(p)));
							sb_cys_site.append(diffFeature+","+sumFeature+",");
							
//							double diffFeature = Math.abs((firstCysFeatures.get(p) - secCysFeatures.get(p))*(firstCysFeatures.get(p) - secCysFeatures.get(p)));
//							double sumFeature = Math.abs((firstCysFeatures.get(p) + secCysFeatures.get(p))*(firstCysFeatures.get(p) + secCysFeatures.get(p)));
//							sb_cys_site.append(diffFeature+","+sumFeature+",");
							
						}
						
						
					}
					
//					double cysBondProbDiff = Math.abs(firstCysBondProb - secCysBondProb);
//					double cysBondProbSum = Math.abs(firstCysBondProb + secCysBondProb);
//					sb_cys_site.append(cysBondProbDiff+","+cysBondProbSum+",");
//					cysBondProbDiff = 2* cysBondProbDiff;
//					cysBondProbSum = 2* cysBondProbSum;
					
					
					// first append before features then middle cysteine features and finally the after features
					sb.append(sb_before.toString()+sb_cys_site.toString()+sb_after.toString());
					
//					if(aaIndexVectorFirst.size() != 2){
//						
//						System.out.println(aaIndexVectorFirst.size());
//					}
					
//					if(aaIndexVectorSec.size() != 2){
//						
//						System.out.println(pdbId);
//						System.out.println(aaIndexVectorSec.size());
//						
//					}
					
					
					for(int p = 0; p < aaIndexVectorFirst.size(); p++){
						
						sb.append(aaIndexVectorFirst.get(p)+",");
						
					}
					
					for(int p = 0; p < aaIndexVectorSec.size(); p++){
						
						sb.append(aaIndexVectorSec.get(p)+",");
						
					}
					
				//	double cod_nor_length = cysDistance/(double)featuresArray.size();
				//	double cod_nor_log = cysDistance/Math.log((double)featuresArray.size());
				//	double cod_nor_log10 = cysDistance/Math.log10((double)featuresArray.size());
//					cysDistance = 2*cysDistance;
					sb.append(cysDistance);
				//	sb.append(cod_nor_length+",");
				//	sb.append(cod_nor_log+",");
				//	sb.append(cod_nor_log10);
					sb.append("\n");
					
					
				}else{ // if window size is 1
					
					String[] firstFeature = featuresArray.get(cysFirstSeqInd).trim().split(",");
					String[] secFeature = featuresArray.get(cysSecSeqInd).trim().split(",");
					List<Double> firstCysFeatures = new ArrayList<Double>();
					List<Double> secCysFeatures = new ArrayList<Double>();
					
					for(int p = 0; p < firstFeature.length; p++){
										
												
//						if(p == 35 || p == 36){
//							
//							continue;
//						}
						
						// sb.append(feature[1]+",");
						firstCysFeatures.add(Double.parseDouble(firstFeature[p]));
						
					}
					
					for(int p = 0; p < secFeature.length; p++){
								
						
//						if(p == 35 || p == 36){
//							
//							continue;
//						}
						
						// sb.append(feature[1]+",");
						secCysFeatures.add(Double.parseDouble(secFeature[p]));
						
					}
					
					for(int p = 0; p < firstCysFeatures.size(); p++){
						
						
						if(p >= 1 && p <= 8){ // ignore the physical properties of the cysteines as feature because they will be same for both the cysteins
							
							continue;
						}
						
						if(p == 0){	// disorder feature
							
							double diffFeature = Math.abs((firstCysFeatures.get(p) - secCysFeatures.get(p)));
							double sumFeature = Math.abs((firstCysFeatures.get(p) + secCysFeatures.get(p)));
							sb_cys_site.append(diffFeature+","+sumFeature+",");
							
//							double diffFeature = Math.abs((firstCysFeatures.get(p) - secCysFeatures.get(p))*(firstCysFeatures.get(p) - secCysFeatures.get(p)));
//							double sumFeature = Math.abs((firstCysFeatures.get(p) + secCysFeatures.get(p))*(firstCysFeatures.get(p) + secCysFeatures.get(p)));
//							sb_cys_site.append(diffFeature+","+sumFeature+",");
							
						}
						
						if(p >= 9 && p <= 28){	// PSSM features
							
							double diffFeature = Math.abs((firstCysFeatures.get(p) - secCysFeatures.get(p)));
							double sumFeature = Math.abs((firstCysFeatures.get(p) + secCysFeatures.get(p)));
							sb_cys_site.append(diffFeature+","+sumFeature+",");
							
//							double diffFeature = Math.abs((firstCysFeatures.get(p) - secCysFeatures.get(p))*(firstCysFeatures.get(p) - secCysFeatures.get(p)));
//							double sumFeature = Math.abs((firstCysFeatures.get(p) + secCysFeatures.get(p))*(firstCysFeatures.get(p) + secCysFeatures.get(p)));
//							sb_cys_site.append(diffFeature+","+sumFeature+",");
												
							
						}
						
						if(p >28){	// all other features
							
							double diffFeature = Math.abs((firstCysFeatures.get(p) - secCysFeatures.get(p)));
							double sumFeature = Math.abs((firstCysFeatures.get(p) + secCysFeatures.get(p)));
							sb_cys_site.append(diffFeature+","+sumFeature+",");
							
//							double diffFeature = Math.abs((firstCysFeatures.get(p) - secCysFeatures.get(p))*(firstCysFeatures.get(p) - secCysFeatures.get(p)));
//							double sumFeature = Math.abs((firstCysFeatures.get(p) + secCysFeatures.get(p))*(firstCysFeatures.get(p) + secCysFeatures.get(p)));
//							sb_cys_site.append(diffFeature+","+sumFeature+",");
												
							
						}
						
					}
					
				
//					double cysBondProbDiff = Math.abs(firstCysBondProb - secCysBondProb);
//					double cysBondProbSum = Math.abs(firstCysBondProb + secCysBondProb);
//					sb_cys_site.append(cysBondProbDiff+","+cysBondProbSum+",");
					
					sb.append(sb_cys_site.toString());
					
				//	double cod_nor_length = cysDistance/(double)featuresArray.size();
				//	double cod_nor_log = cysDistance/Math.log((double)featuresArray.size());
				//	double cod_nor_log10 = cysDistance/Math.log10((double)featuresArray.size());
					
				//	sb.append(","+cysDistance);
				//	sb.append(","+cod_nor_length);
				//	sb.append(","+cod_nor_log);
				//	sb.append(","+(cod_nor_log10));
					sb.append(cysDistance);
					sb.append("\n");
					
					
				}
				
//				if(classType == 1){
//					
//					positive_features.add(sb.toString());
//					
//				}else if(classType == 0){
//					
//					negative_features.add(sb.toString());
//					
//				}
				
				winWr.write(sb.toString());
				winWr.flush();
			
			}
				
		}
		
		br_id.close();
		
//		for(int p = 0; p < positive_features.size(); p++){
//			
//			formatted_features.add(positive_features.get(p));
//			
//		}
//		
//		if(negative_features.size() >= 5*positive_features.size()){
//			
//			for(int p = 0; p < positive_features.size()*5; p++){
//				
//				formatted_features.add(negative_features.get(p));
//				
//			}
//			
//			
//		}else{
//			
//			for(int p = 0; p < negative_features.size(); p++){
//				
//				formatted_features.add(negative_features.get(p));
//				
//			}
//						
//		}
		
		
		
		
//		System.out.println("number of positive samples "+positive_features.size());
//		System.out.println("number of negative samples "+positive_features.size()*5);
//		System.out.println("actual number of negative samples "+negative_features.size());
//		System.out.println("number of total samples "+formatted_features.size());
		System.out.println("positiveCysPairsCount "+positiveCysPairsCount);
		System.out.println("positiveCysClassCount "+positiveClassCount);
		System.out.println("negativeCysPairsCount "+negativeCysPairsCount);
		System.out.println("negativeCysClassCount "+negativeClassCount);
		System.out.println("percentage of positive samples with distance < 100  "+((positiveCysPairsCount/positiveClassCount)*100));
		System.out.println("percentage of negative samples with distance < 100  "+((negativeCysPairsCount/negativeClassCount)*100));
		
		
		
//		Collections.shuffle(formatted_features);
//		for(int s = 0; s < formatted_features.size(); s++){
//			
//			winWr.write(formatted_features.get(s));
//			winWr.flush();
//			
//		}
		
		winWr.close();
		
		
	}
	
	public static void createWindowForCysPairForTestData(String idListPathTest, String fastaFileDirTest, String cysPairWindowedFilesDirTest, String combinedFeatureDirTest, String singleCysPredProbPathTest, int ws, int num_of_features, double[][] aaFactors) throws IOException{
		
				
//		winWr.write("class,");
//		
//		for(int i = 0; i < 1051; i++){
//			
//			if(i == 1050){
//				
//				winWr.write("V"+(i+1));
//				
//			}else{
//				
//				winWr.write("V"+(i+1)+",");
//				
//			}			
//			
//		}
//		
//		winWr.write("\n");
		
		List<String> formatted_features = new ArrayList<String>();
		List<String> positive_features = new ArrayList<String>();
		List<String> negative_features = new ArrayList<String>();
		
		List<Integer> positiveCysPairDistance = new ArrayList<Integer>();
		List<Integer> negativeCysPairDistance = new ArrayList<Integer>();
		
		double positiveClassCount = 0;
		double negativeClassCount = 0;
		double positiveCysPairsCount = 0;
		double negativeCysPairsCount = 0;
		
		
		BufferedReader br_id = BufferReaderAndWriter.getReader(new File(idListPathTest));
		String pdbId = "";
		
		while((pdbId = br_id.readLine())!=null){
			
			System.out.println("pdbId "+pdbId);
			// parse the connectivity file to obtain cys pairs and the total number of cysteins present in the structure
//			BufferedReader brConn = BufferReaderAndWriter.getReader(new File(cysConnectionPath+pdbId+".ssbond"));
//			String connLine = "";
//			
//			List<String> cysPairs = new ArrayList<String>(); // contains the disulfide cysteine pairs
//			int totalCys = -1;  // total number of cys present in the structure
//			
//			while((connLine = brConn.readLine())!=null){
//				
//				if(connLine.startsWith("RltvPos")){
//					
//					String[] elements = connLine.split("\\s+");
//					cysPairs.add(elements[1]);
//									
//				}else if(connLine.startsWith("TotNumCys")){
//					
//					String[] elements = connLine.split("\\s+");
//					totalCys = Integer.parseInt(elements[1]);
//					
//				}
//				
//			}
//			
//			brConn.close();
			
			// for test data cysteine connection files are not available so ... totalCys is extracted from fasta sequence
			int totalCys = 0;
			String fastaFile = fastaFileDirTest+pdbId+".fasta";
			// System.out.println(fastaFile);
			BufferedReader fastaRd = BufferReaderAndWriter.getReader(new File(fastaFile));
			String fastaLine = "";
			while((fastaLine = fastaRd.readLine())!=null){
				
				if(fastaLine.startsWith(">")){
					continue;
				}
				
				char[] aa = fastaLine.toCharArray();
				for(int i = 0; i < aa.length; i++){
					
					if(aa[i] == 'C'){
						
						totalCys += 1;
						
					}
					
				}
				
			}
			
			fastaRd.close();
			
			
			/* ================== create total possible combination of cysteine pairs for the structure ======================== */

			List<String> cysCombs = new ArrayList<String>();	// holds all possible combination of the CYS in a protein
			for(int i = 1; i < totalCys; i++){
				
				for(int j = i+1; j <= totalCys; j++){
					
//					int cysDistance = Math.abs(j - i);
//					
//					if(cysDistance >= 100){
//						
//						continue;
//					}
					
					cysCombs.add("c"+i+"-"+"c"+j);
//					System.out.println("c"+i+"-"+"c"+j);
					
					
				}
				
				
			}
			
			
/* ================== parse the single cysteine predict file and collect cysteine bonding prob to use it as a feature ======================== */
			
			Map<String, Double> cysBondProb = new HashMap<String, Double>();
			
			String cysBondPredFilePath = singleCysPredProbPathTest+pdbId+".predict";
			BufferedReader cysBondPredFileRd = BufferReaderAndWriter.getReader(new File(cysBondPredFilePath));
			String cysBondLine = "";
			int cysId = 0;
//			int classOneCount = 0;
//			System.out.println(pdbId);
			while((cysBondLine = cysBondPredFileRd.readLine())!=null){
				
				if(cysBondLine.startsWith("predClass")){
					
					continue;
					
				}
				cysId++;
				String[] values = cysBondLine.split("\\s+");
				cysBondProb.put("c"+cysId, Double.parseDouble(values[2]));
				
//				if(Integer.parseInt(values[0]) == 1){
//					
//					classOneCount++;
//					
//				}
				
				
			}
						
			cysBondPredFileRd.close();
			
			
/* ================== parse the features file ======================== */
			
			// loading file in memory
			List<String> featuresArray = new ArrayList<String>();	// holds features for each residues
			Map<String, Integer> cysInd = new HashMap<String, Integer>();	 	// holds the cysteine relative position and sequence index e.g. key = c1 and Val = 10
						
			String featureFilePath = combinedFeatureDirTest+pdbId+".61features";
			
			BufferedReader featureFileReader = BufferReaderAndWriter.getReader(new File(featureFilePath));
			
			String featureLine = "";
			int cysRelPos = 0;
			
			while((featureLine = featureFileReader.readLine())!=null){
			
				
				if(featureLine.startsWith("disProb(1)")){
					
					continue;
				}
				
				featuresArray.add(featureLine);
				
				String[] featureElements = featureLine.split(",");
				String aa_type = featureElements[1];
				if(Integer.parseInt(aa_type) == 5){
					
					cysRelPos++;					
					cysInd.put("c"+cysRelPos, featuresArray.size()-1);
										
				}			
			
			}
			
			featureFileReader.close();
			
			// create windowed output file
			String cysPairWindowedFileName = cysPairWindowedFilesDirTest+pdbId+"_cyspair_incsingcysprob_ws_"+ws+".txt";
			BufferedWriter winWr = BufferReaderAndWriter.getWriter(new File(cysPairWindowedFileName));
			
/* ================== Collect Windowed Feature for all possible combination of cys pair and assign the class label ======================== */
			
			for(int i = 0; i < cysCombs.size(); i++){
				
				String cysComb = cysCombs.get(i);
				int classType = 0;
				
//				if(cysPairs.contains(cysComb)){
//					
//					classType = +1;
//					positiveClassCount += 1;
//					
//				}else{
//					
//					classType = 0;
//					negativeClassCount += 1;
//					
//				}			
				
				
				String[] cysCombArray = cysComb.split("-");
				String cysCombFirst = cysCombArray[0];
				String cysCombSec = cysCombArray[1];
				
				double firstCysBondProb = cysBondProb.get(cysCombFirst);
				double secCysBondProb = cysBondProb.get(cysCombSec);
				
				int cysFirstSeqInd = cysInd.get(cysCombFirst);
				int cysSecSeqInd = cysInd.get(cysCombSec);
				
				int cysDistance = Math.abs(cysSecSeqInd - cysFirstSeqInd);
				
				if(cysDistance >= 100){
				
					continue;
				
				}
				
				System.out.println(cysComb);
				
//				if(classType == 1 && cysDistance < 100){
//					
////					positiveCysPairDistance.add(cysDistance);
//					positiveCysPairsCount += 1; 
//					
//				}else if(classType == 0 && cysDistance < 100){
//					
////					negativeCysPairDistance.add(cysDistance);
//					negativeCysPairsCount += 1;
//					
//				}
//				
//				
//				
//				if(classType == 0 && cysDistance >= 100){
//					
//					continue;
//					
//				}
//				
//				// for 5-fold negative sample and positive samples turn this condition on
////				if(classType == 0 && negativeCysPairsCount > 43935){
////					
////					continue;
////					
////				}
//				
//				// for balanced negative and positive samples turn this condition on
//				if(classType == 0 && negativeCysPairsCount > 8787){
//					
//					continue;
//					
//				}
				
				
				StringBuilder sb = new StringBuilder();
				StringBuilder sb_before = new StringBuilder();
				StringBuilder sb_after = new StringBuilder();
				StringBuilder sb_cys_site = new StringBuilder();
				
				sb.append(Integer.toString(classType)+",");
				
				int beforeAndAfterRowsToBeIncluded = ws/2;
				
//				/*=================== Ignore the cystein pairs, either of which can not form consecutive residues equal to the window size =============*/
//				int beforeWindowCheck = cysFirstSeqInd - beforeAndAfterRowsToBeIncluded;
//				int afterWindowCheck = cysSecSeqInd + beforeAndAfterRowsToBeIncluded;
//				
//				if((beforeWindowCheck < 0) || (afterWindowCheck >= featuresArray.size())){
//					System.out.println(pdbId+"\t"+featuresArray.size()+"\t"+cysFirstSeqInd+"\t"+cysSecSeqInd);
//					continue;					
//					
//				}
				
				if(beforeAndAfterRowsToBeIncluded > 0){		// if window size is > 1
					
					int beforeFeatureCounter = beforeAndAfterRowsToBeIncluded;
					int afterFeatureCounter = 0;
					
					List<Integer> aaIndexVectorFirst = new ArrayList<Integer>();
					List<Integer> aaIndexVectorSec = new ArrayList<Integer>();
					
					for(int k = 0; k < beforeAndAfterRowsToBeIncluded; k++){		// append the features to the left
						afterFeatureCounter++;												
						int firstBeforeIndex = cysFirstSeqInd-beforeFeatureCounter;
						int secBeforeIndex = cysSecSeqInd-beforeFeatureCounter;
						int firstAfterIndex = cysFirstSeqInd+afterFeatureCounter;
						int secAfterIndex = cysSecSeqInd+afterFeatureCounter;
						List<Double> firstBeforeFeatures = new ArrayList<Double>();
						List<Double> secBeforeFeatures = new ArrayList<Double>();
						List<Double> firstAfterFeatures = new ArrayList<Double>();
						List<Double> secAfterFeatures = new ArrayList<Double>();
												
												
						if(firstBeforeIndex < 0){										// if for terminal residues features does not exist
							
							for(int l = 0; l < num_of_features; l++){
							//	featureCounter++;
								firstBeforeFeatures.add(0.0);
								
							}
							aaIndexVectorFirst.add(0);
							
						}else {
							
							String[] beforeFeatureExist = featuresArray.get(firstBeforeIndex).trim().split(",");
							
							for(int m = 0; m < beforeFeatureExist.length; m++){
								
								if(m == 1){ // ignore AA based feature
																		
									aaIndexVectorFirst.add(Integer.parseInt(beforeFeatureExist[m]));
								//	continue;
								}
							//	featureCounter++;
							// 	featuresString = featuresString+featureCounter+":"+beforeFeatureExist[m]+" ";
							//	featuresString = featuresString+beforeFeatureExist[m]+" ";
								
								firstBeforeFeatures.add(Double.parseDouble(beforeFeatureExist[m]));
								
							}
							
							
						}
						
						if(secBeforeIndex < 0){
							
							for(int l = 0; l < num_of_features; l++){
								//	featureCounter++;
									secBeforeFeatures.add(0.0);
									
							}
							
							aaIndexVectorSec.add(0);
							
							
						}else{
							
							String[] beforeFeatureExist = featuresArray.get(secBeforeIndex).trim().split(",");
							
							for(int m = 0; m < beforeFeatureExist.length; m++){
								if(m == 1){	// ignore aa based feature
																		
									aaIndexVectorSec.add(Integer.parseInt(beforeFeatureExist[m]));
//									continue;
								}
							//	featureCounter++;
							// 	featuresString = featuresString+featureCounter+":"+beforeFeatureExist[m]+" ";
							//	featuresString = featuresString+beforeFeatureExist[m]+" ";
								
								secBeforeFeatures.add(Double.parseDouble(beforeFeatureExist[m]));
								
							}
							
						}
						
						if(firstAfterIndex > featuresArray.size()-1){										// if for terminal residues features does not exist
							
							for(int l = 0; l < num_of_features; l++){
							//	featureCounter++;
								firstAfterFeatures.add(0.0);
								
							}
							
							aaIndexVectorFirst.add(0);
							
						}else {
							
							String[] afterFeatureExist = featuresArray.get(firstAfterIndex).trim().split(",");
							
							for(int m = 0; m < afterFeatureExist.length; m++){
								if(m == 1){	// ignore aa based feature
																		
									aaIndexVectorFirst.add(Integer.parseInt(afterFeatureExist[m]));
//									continue;
								}
							//	featureCounter++;
							// 	featuresString = featuresString+featureCounter+":"+beforeFeatureExist[m]+" ";
							//	featuresString = featuresString+beforeFeatureExist[m]+" ";
								
								firstAfterFeatures.add(Double.parseDouble(afterFeatureExist[m]));
								
							}
							
							
						}
						
						if(secAfterIndex > featuresArray.size()-1){
							
							for(int l = 0; l < num_of_features; l++){
								//	featureCounter++;
									secAfterFeatures.add(0.0);
									
							}
							
							aaIndexVectorSec.add(0);
							
							
						}else{
							
							String[] afterFeatureExist = featuresArray.get(secAfterIndex).trim().split(",");
							
							for(int m = 0; m < afterFeatureExist.length; m++){
								if(m == 1){	// ignore aa based feature
																		
									aaIndexVectorSec.add(Integer.parseInt(afterFeatureExist[m]));
//									continue;
								}
							//	featureCounter++;
							// 	featuresString = featuresString+featureCounter+":"+beforeFeatureExist[m]+" ";
							//	featuresString = featuresString+beforeFeatureExist[m]+" ";
								
								secAfterFeatures.add(Double.parseDouble(afterFeatureExist[m]));
								
							}
							
						}
						
						for(int p = 0; p < firstBeforeFeatures.size(); p++){ // firstBeforeFeatures contains all the features except the class label information
							
//							if((p >= 2 && p <= 8) || p > 28){
//								
//								continue;
//								
//							}
							
							if((p >= 2 && p <= 8)){	 // ignore the physical parameter features
								
//								double diffFeature = Math.abs((firstBeforeFeatures.get(p) - secBeforeFeatures.get(p)));
//								double sumFeature = Math.abs((firstBeforeFeatures.get(p) + secBeforeFeatures.get(p)));
//								sb_before.append(diffFeature+","+sumFeature+",");
								continue;
								
							}
							
							if(p == 0){	// disorder feature
								
								double diffFeature = Math.abs((firstBeforeFeatures.get(p) - secBeforeFeatures.get(p)));
								double sumFeature = Math.abs((firstBeforeFeatures.get(p) + secBeforeFeatures.get(p)));
								sb_before.append(diffFeature+","+sumFeature+",");
								
//								double diffFeature = Math.abs((firstBeforeFeatures.get(p) - secBeforeFeatures.get(p))*(firstBeforeFeatures.get(p) - secBeforeFeatures.get(p)));
//								double sumFeature = Math.abs((firstBeforeFeatures.get(p) + secBeforeFeatures.get(p))*(firstBeforeFeatures.get(p) + secBeforeFeatures.get(p)));
//								sb_before.append(diffFeature+","+sumFeature+",");
								
							}
							
							if(p == 1){	// aaFactors features
								
								int firstAAIndex =  firstBeforeFeatures.get(p).intValue();
								int secAAIndex = secBeforeFeatures.get(p).intValue();
								
								double[] firstAAFactors = {0, 0, 0, 0, 0};
								double[] secAAFactors = {0, 0, 0, 0, 0};
								if(!(firstAAIndex == 0)){
									
									firstAAFactors =  aaFactors[firstAAIndex-1];
									
								}
								
								if(!(secAAIndex == 0)){
									
									secAAFactors = aaFactors[secAAIndex-1];
									
								}
								
								for(int q = 0; q < 5; q++){
									
									double diffFeature = Math.abs(firstAAFactors[q] - secAAFactors[q]);
									double sumFeature = Math.abs(firstAAFactors[q] + secAAFactors[q]);
									sb_before.append(diffFeature+","+sumFeature+",");
									
//									double diffFeature = Math.abs((firstAAFactors[q] - secAAFactors[q])*(firstAAFactors[q] - secAAFactors[q]));
//									double sumFeature = Math.abs((firstAAFactors[q] + secAAFactors[q])*(firstAAFactors[q] + secAAFactors[q]));
//									sb_before.append(diffFeature+","+sumFeature+",");
									
									
								}				
								
								
							}						
														
							
							if(p >= 9 && p <= 28){	// PSSM features
								
								double diffFeature = Math.abs((firstBeforeFeatures.get(p) - secBeforeFeatures.get(p)));
								double sumFeature = Math.abs((firstBeforeFeatures.get(p) + secBeforeFeatures.get(p)));
								sb_before.append(diffFeature+","+sumFeature+",");
								
//								double diffFeature = Math.abs((firstBeforeFeatures.get(p) - secBeforeFeatures.get(p))*(firstBeforeFeatures.get(p) - secBeforeFeatures.get(p)));
//								double sumFeature = Math.abs((firstBeforeFeatures.get(p) + secBeforeFeatures.get(p))*(firstBeforeFeatures.get(p) + secBeforeFeatures.get(p)));
//								sb_before.append(diffFeature+","+sumFeature+",");
															
								
							}
							
							if(p > 28){ // all other features
								
								double diffFeature = Math.abs((firstBeforeFeatures.get(p) - secBeforeFeatures.get(p)));
								double sumFeature = Math.abs((firstBeforeFeatures.get(p) + secBeforeFeatures.get(p)));
								sb_before.append(diffFeature+","+sumFeature+",");
								
//								double diffFeature = Math.abs((firstBeforeFeatures.get(p) - secBeforeFeatures.get(p))*(firstBeforeFeatures.get(p) - secBeforeFeatures.get(p)));
//								double sumFeature = Math.abs((firstBeforeFeatures.get(p) + secBeforeFeatures.get(p))*(firstBeforeFeatures.get(p) + secBeforeFeatures.get(p)));
//								sb_before.append(diffFeature+","+sumFeature+",");
								
							}
							
							
						}
						
						for(int p = 0; p < firstAfterFeatures.size(); p++){
							
//							if((p >= 2 && p <= 8) || p > 28){
//								
//								continue;
//								
//							}
							
							if((p >= 2 && p <= 8)){	// ignore the physical parameter feature
							
//								double diffFeature = Math.abs((firstAfterFeatures.get(p) - secAfterFeatures.get(p)));
//								double sumFeature = Math.abs((firstAfterFeatures.get(p) + secAfterFeatures.get(p)));
//								sb_after.append(diffFeature+","+sumFeature+",");
								continue;
								
							}
							
							if(p == 0){	// disorder feature
								
								double diffFeature = Math.abs((firstAfterFeatures.get(p) - secAfterFeatures.get(p)));
								double sumFeature = Math.abs((firstAfterFeatures.get(p) + secAfterFeatures.get(p)));
								sb_after.append(diffFeature+","+sumFeature+",");
								
//								double diffFeature = Math.abs((firstAfterFeatures.get(p) - secAfterFeatures.get(p))*(firstAfterFeatures.get(p) - secAfterFeatures.get(p)));
//								double sumFeature = Math.abs((firstAfterFeatures.get(p) + secAfterFeatures.get(p))*(firstAfterFeatures.get(p) + secAfterFeatures.get(p)));
//								sb_after.append(diffFeature+","+sumFeature+",");
								
							}
							
							if(p == 1){	// aaFactors features
								
								int firstAAIndex =  firstAfterFeatures.get(p).intValue();
								int secAAIndex = secAfterFeatures.get(p).intValue();
								
								double[] firstAAFactors = {0, 0, 0, 0, 0};
								double[] secAAFactors = {0, 0, 0, 0, 0};
								if(!(firstAAIndex == 0)){
									
									firstAAFactors =  aaFactors[firstAAIndex-1];
									
								}
								
								if(!(secAAIndex == 0)){
									
									secAAFactors = aaFactors[secAAIndex-1];
									
								}								
								
								
								for(int q = 0; q < 5; q++){
									
									double diffFeature = Math.abs(firstAAFactors[q] - secAAFactors[q]);
									double sumFeature = Math.abs(firstAAFactors[q] + secAAFactors[q]);
									sb_after.append(diffFeature+","+sumFeature+",");
									
//									double diffFeature = Math.abs((firstAAFactors[q] - secAAFactors[q])*(firstAAFactors[q] - secAAFactors[q]));
//									double sumFeature = Math.abs((firstAAFactors[q] + secAAFactors[q])*(firstAAFactors[q] + secAAFactors[q]));
//									sb_after.append(diffFeature+","+sumFeature+",");
									
									
								}				
								
								
							}						
														
							
							if(p >= 9 && p <= 28){	// PSSM features
								
								double diffFeature = Math.abs((firstAfterFeatures.get(p) - secAfterFeatures.get(p)));
								double sumFeature = Math.abs((firstAfterFeatures.get(p) + secAfterFeatures.get(p)));
								sb_after.append(diffFeature+","+sumFeature+",");
								
//								double diffFeature = Math.abs((firstAfterFeatures.get(p) - secAfterFeatures.get(p))*(firstAfterFeatures.get(p) - secAfterFeatures.get(p)));
//								double sumFeature = Math.abs((firstAfterFeatures.get(p) + secAfterFeatures.get(p))*(firstAfterFeatures.get(p) + secAfterFeatures.get(p)));
//								sb_after.append(diffFeature+","+sumFeature+",");
																
								
							}
							
							if(p > 28){ // all other features
								
								double diffFeature = Math.abs((firstAfterFeatures.get(p) - secAfterFeatures.get(p)));
								double sumFeature = Math.abs((firstAfterFeatures.get(p) + secAfterFeatures.get(p)));
								sb_after.append(diffFeature+","+sumFeature+",");
								
//								double diffFeature = Math.abs((firstAfterFeatures.get(p) - secAfterFeatures.get(p))*(firstAfterFeatures.get(p) - secAfterFeatures.get(p)));
//								double sumFeature = Math.abs((firstAfterFeatures.get(p) + secAfterFeatures.get(p))*(firstAfterFeatures.get(p) + secAfterFeatures.get(p)));
//								sb_after.append(diffFeature+","+sumFeature+",");
								
							}
							
							
						}
						
												
						beforeFeatureCounter--;
												
					}
					
					String[] firstFeature = featuresArray.get(cysFirstSeqInd).trim().split(",");
					String[] secFeature = featuresArray.get(cysSecSeqInd).trim().split(",");
					List<Double> firstCysFeatures = new ArrayList<Double>();
					List<Double> secCysFeatures = new ArrayList<Double>();
					
					for(int p = 0; p < firstFeature.length; p++){
											
						
//						if(p == 35 || p == 36){	// ignore dpsi and MG features
//							
//							continue;
//						}
						
						// sb.append(feature[1]+",");
						firstCysFeatures.add(Double.parseDouble(firstFeature[p]));
						
					}
					
					for(int p = 0; p < secFeature.length; p++){
									
						
//						if(p == 35 || p == 36){	// ignore dpsi and MG features
//							
//							continue;
//						}
									
						// sb.append(feature[1]+",");
						secCysFeatures.add(Double.parseDouble(secFeature[p]));
						
					}
					
					for(int p = 0; p < firstCysFeatures.size(); p++){
						
						
//						if(p > 20){
//							
//							continue;
//							
//						}
						
						if(p >= 1 && p <= 8){ // ignore the residue index and physical properties of the cysteines as feature because they will be same for both the cysteines
							
							continue;
						}
						
						if(p == 0){	// disorder feature
							
							double diffFeature = Math.abs((firstCysFeatures.get(p) - secCysFeatures.get(p)));
							double sumFeature = Math.abs((firstCysFeatures.get(p) + secCysFeatures.get(p)));
							sb_cys_site.append(diffFeature+","+sumFeature+",");
							
//							double diffFeature = Math.abs((firstCysFeatures.get(p) - secCysFeatures.get(p))*(firstCysFeatures.get(p) - secCysFeatures.get(p)));
//							double sumFeature = Math.abs((firstCysFeatures.get(p) + secCysFeatures.get(p))*(firstCysFeatures.get(p) + secCysFeatures.get(p)));
//							sb_cys_site.append(diffFeature+","+sumFeature+",");
							
						}
						
						if(p >= 9 && p <= 28){	// PSSM features
							
							double diffFeature = Math.abs((firstCysFeatures.get(p) - secCysFeatures.get(p)));
							double sumFeature = Math.abs((firstCysFeatures.get(p) + secCysFeatures.get(p)));
							sb_cys_site.append(diffFeature+","+sumFeature+",");
							
//							double diffFeature = Math.abs((firstCysFeatures.get(p) - secCysFeatures.get(p))*(firstCysFeatures.get(p) - secCysFeatures.get(p)));
//							double sumFeature = Math.abs((firstCysFeatures.get(p) + secCysFeatures.get(p))*(firstCysFeatures.get(p) + secCysFeatures.get(p)));
//							sb_cys_site.append(diffFeature+","+sumFeature+",");
												
							
						}
						
						if(p > 28){	// consider all other features too
							
							double diffFeature = Math.abs((firstCysFeatures.get(p) - secCysFeatures.get(p)));
							double sumFeature = Math.abs((firstCysFeatures.get(p) + secCysFeatures.get(p)));
							sb_cys_site.append(diffFeature+","+sumFeature+",");
							
//							double diffFeature = Math.abs((firstCysFeatures.get(p) - secCysFeatures.get(p))*(firstCysFeatures.get(p) - secCysFeatures.get(p)));
//							double sumFeature = Math.abs((firstCysFeatures.get(p) + secCysFeatures.get(p))*(firstCysFeatures.get(p) + secCysFeatures.get(p)));
//							sb_cys_site.append(diffFeature+","+sumFeature+",");
							
						}
						
						
					}
					
					double cysBondProbDiff = Math.abs(firstCysBondProb - secCysBondProb);
					double cysBondProbSum = Math.abs(firstCysBondProb + secCysBondProb);
					sb_cys_site.append(cysBondProbDiff+","+cysBondProbSum+",");
//					cysBondProbDiff = 2* cysBondProbDiff;
//					cysBondProbSum = 2* cysBondProbSum;
					
					
					// first append before features then middle cysteine features and finally the after features
					sb.append(sb_before.toString()+sb_cys_site.toString()+sb_after.toString());
					
//					if(aaIndexVectorFirst.size() != 2){
//						
//						System.out.println(aaIndexVectorFirst.size());
//					}
					
//					if(aaIndexVectorSec.size() != 2){
//						
//						System.out.println(pdbId);
//						System.out.println(aaIndexVectorSec.size());
//						
//					}
					
					
					for(int p = 0; p < aaIndexVectorFirst.size(); p++){
						
						sb.append(aaIndexVectorFirst.get(p)+",");
						
					}
					
					for(int p = 0; p < aaIndexVectorSec.size(); p++){
						
						sb.append(aaIndexVectorSec.get(p)+",");
						
					}
					
				//	double cod_nor_length = cysDistance/(double)featuresArray.size();
				//	double cod_nor_log = cysDistance/Math.log((double)featuresArray.size());
				//	double cod_nor_log10 = cysDistance/Math.log10((double)featuresArray.size());
//					cysDistance = 2*cysDistance;
					sb.append(cysDistance);
				//	sb.append(cod_nor_length+",");
				//	sb.append(cod_nor_log+",");
				//	sb.append(cod_nor_log10);
					sb.append("\n");
					
					
				}else{ // if window size is 1
					
					String[] firstFeature = featuresArray.get(cysFirstSeqInd).trim().split(",");
					String[] secFeature = featuresArray.get(cysSecSeqInd).trim().split(",");
					List<Double> firstCysFeatures = new ArrayList<Double>();
					List<Double> secCysFeatures = new ArrayList<Double>();
					
					for(int p = 0; p < firstFeature.length; p++){
										
												
//						if(p == 35 || p == 36){
//							
//							continue;
//						}
						
						// sb.append(feature[1]+",");
						firstCysFeatures.add(Double.parseDouble(firstFeature[p]));
						
					}
					
					for(int p = 0; p < secFeature.length; p++){
								
						
//						if(p == 35 || p == 36){
//							
//							continue;
//						}
						
						// sb.append(feature[1]+",");
						secCysFeatures.add(Double.parseDouble(secFeature[p]));
						
					}
					
					for(int p = 0; p < firstCysFeatures.size(); p++){
						
						
						if(p >= 1 && p <= 8){ // ignore the physical properties of the cysteines as feature because they will be same for both the cysteins
							
							continue;
						}
						
						if(p == 0){	// disorder feature
							
							double diffFeature = Math.abs((firstCysFeatures.get(p) - secCysFeatures.get(p)));
							double sumFeature = Math.abs((firstCysFeatures.get(p) + secCysFeatures.get(p)));
							sb_cys_site.append(diffFeature+","+sumFeature+",");
							
//							double diffFeature = Math.abs((firstCysFeatures.get(p) - secCysFeatures.get(p))*(firstCysFeatures.get(p) - secCysFeatures.get(p)));
//							double sumFeature = Math.abs((firstCysFeatures.get(p) + secCysFeatures.get(p))*(firstCysFeatures.get(p) + secCysFeatures.get(p)));
//							sb_cys_site.append(diffFeature+","+sumFeature+",");
							
						}
						
						if(p >= 9 && p <= 28){	// PSSM features
							
							double diffFeature = Math.abs((firstCysFeatures.get(p) - secCysFeatures.get(p)));
							double sumFeature = Math.abs((firstCysFeatures.get(p) + secCysFeatures.get(p)));
							sb_cys_site.append(diffFeature+","+sumFeature+",");
							
//							double diffFeature = Math.abs((firstCysFeatures.get(p) - secCysFeatures.get(p))*(firstCysFeatures.get(p) - secCysFeatures.get(p)));
//							double sumFeature = Math.abs((firstCysFeatures.get(p) + secCysFeatures.get(p))*(firstCysFeatures.get(p) + secCysFeatures.get(p)));
//							sb_cys_site.append(diffFeature+","+sumFeature+",");
												
							
						}
						
						if(p >28){	// all other features
							
							double diffFeature = Math.abs((firstCysFeatures.get(p) - secCysFeatures.get(p)));
							double sumFeature = Math.abs((firstCysFeatures.get(p) + secCysFeatures.get(p)));
							sb_cys_site.append(diffFeature+","+sumFeature+",");
							
//							double diffFeature = Math.abs((firstCysFeatures.get(p) - secCysFeatures.get(p))*(firstCysFeatures.get(p) - secCysFeatures.get(p)));
//							double sumFeature = Math.abs((firstCysFeatures.get(p) + secCysFeatures.get(p))*(firstCysFeatures.get(p) + secCysFeatures.get(p)));
//							sb_cys_site.append(diffFeature+","+sumFeature+",");
												
							
						}
						
					}
					
				
					double cysBondProbDiff = Math.abs(firstCysBondProb - secCysBondProb);
					double cysBondProbSum = Math.abs(firstCysBondProb + secCysBondProb);
					sb_cys_site.append(cysBondProbDiff+","+cysBondProbSum+",");
					
					sb.append(sb_cys_site.toString());
					
				//	double cod_nor_length = cysDistance/(double)featuresArray.size();
				//	double cod_nor_log = cysDistance/Math.log((double)featuresArray.size());
				//	double cod_nor_log10 = cysDistance/Math.log10((double)featuresArray.size());
					
				//	sb.append(","+cysDistance);
				//	sb.append(","+cod_nor_length);
				//	sb.append(","+cod_nor_log);
				//	sb.append(","+(cod_nor_log10));
					sb.append(cysDistance);
					sb.append("\n");
					
					
				}
				
//				if(classType == 1){
//					
//					positive_features.add(sb.toString());
//					
//				}else if(classType == 0){
//					
//					negative_features.add(sb.toString());
//					
//				}
				
				winWr.write(sb.toString());
				winWr.flush();
			
			}
			
			winWr.close();
				
		}
		
		br_id.close();
		
//		for(int p = 0; p < positive_features.size(); p++){
//			
//			formatted_features.add(positive_features.get(p));
//			
//		}
//		
//		if(negative_features.size() >= 5*positive_features.size()){
//			
//			for(int p = 0; p < positive_features.size()*5; p++){
//				
//				formatted_features.add(negative_features.get(p));
//				
//			}
//			
//			
//		}else{
//			
//			for(int p = 0; p < negative_features.size(); p++){
//				
//				formatted_features.add(negative_features.get(p));
//				
//			}
//						
//		}
		
		
		
		
//		System.out.println("number of positive samples "+positive_features.size());
//		System.out.println("number of negative samples "+positive_features.size()*5);
//		System.out.println("actual number of negative samples "+negative_features.size());
//		System.out.println("number of total samples "+formatted_features.size());
		
//		System.out.println("positiveCysPairsCount "+positiveCysPairsCount);
//		System.out.println("positiveCysClassCount "+positiveClassCount);
//		System.out.println("negativeCysPairsCount "+negativeCysPairsCount);
//		System.out.println("negativeCysClassCount "+negativeClassCount);
//		System.out.println("percentage of positive samples with distance < 100  "+((positiveCysPairsCount/positiveClassCount)*100));
//		System.out.println("percentage of negative samples with distance < 100  "+((negativeCysPairsCount/negativeClassCount)*100));
//		
		
		
//		Collections.shuffle(formatted_features);
//		for(int s = 0; s < formatted_features.size(); s++){
//			
//			winWr.write(formatted_features.get(s));
//			winWr.flush();
//			
//		}	
		
		
	}
	
	
	

}
