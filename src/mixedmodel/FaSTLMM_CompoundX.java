package mixedmodel;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import javax.tools.DocumentationTool.Location;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularValueDecomposition;

public class FaSTLMM_CompoundX extends FaSTLMM_Compound{

	public FaSTLMM_CompoundX(Phenotype ori_phenotype, String genotype_file_hdf5, String global_kinship_file) {
		super(ori_phenotype, genotype_file_hdf5, global_kinship_file);
		// TODO Auto-generated constructor stub
	}
	
	boolean load_one_res_by_location(int chr, int location, String emmax_res){
		Boolean selected_SNP = false;
		try {
			BufferedReader bReader = new BufferedReader(new FileReader(emmax_res));
			String lineString = bReader.readLine();
			while(lineString !=null) {
				if(lineString.startsWith(Integer.toString(chr)+","+Integer.toString(location)+",")) {
					String[] lineStrings = lineString.split(",");
					double p_value = Double.parseDouble(lineStrings[2]);
					if(Double.compare(p_value, 0.05/genotype.num_sites_total)<=0) {
						selected_SNP=true;
						break;
					}
				}
				lineString = bReader.readLine();
			}
			bReader.close();
		}catch (Exception e) {
			e.printStackTrace();
			// TODO: handle exception
		}
		return selected_SNP;
	}
	
	public static String getKeyFromValue(Map<String,Double> map, Double value) {
		Set<String> keys = map.keySet();
		for(String key:keys) {
			if(map.get(key).equals(value)) {
				return key;
			}
		}
		return null;
	}
	
	public double[][] select_SNPs_base_emmax_res(double[][] snp_in_this_region, int m, int chr, int start,int end, String emmax_res) {
		//select top m snps from p.values
		System.out.println(chr+";"+start+";"+end+";"+emmax_res);
		if(snp_in_this_region.length >=m && m>=0) {
			HashMap<String,String> p_location = new HashMap<>();
			int snps_have_p_values=0;
			try {
				BufferedReader bReader = new BufferedReader(new FileReader(emmax_res));
				String lineString = bReader.readLine();
				while(lineString !=null) {
					if(!lineString.startsWith("#")) {
						String[] lineStrings = lineString.split(",");
						if(Integer.parseInt(lineStrings[0])==chr) {
							String loc=lineStrings[1];
							if(Integer.parseInt(loc)>=start && Integer.parseInt(loc)<=end) {
								String p=lineStrings[2];
								//no two identical locations in one chr
								if(p_location.containsKey(p)) {
									String old_loc = p_location.get(p);
									if(!old_loc.equals(loc)) {
										String new_loc = old_loc+"_"+loc;
										//System.out.println(new_loc);
										p_location.put(p, new_loc);
										snps_have_p_values++;
									}
								}else {
									p_location.put(p,loc);
									snps_have_p_values++;
								}
							}
						}else if(Integer.parseInt(lineStrings[0])>chr) {
							break;
						}
					}
					lineString = bReader.readLine();
				}
				bReader.close();
			}catch (Exception e) {
				e.printStackTrace();
				// TODO: handle exception
			}
			System.out.println("p_location size: "+p_location.size());
			Set<String> keySet = p_location.keySet();
			double[] pvalues = new double[keySet.size()];
			int pp=0;
			for(String key:keySet) {
				pvalues[pp]=Double.parseDouble(key);
				pp++;
			}
			Arrays.sort(pvalues);
			int used_snps=0;
			if(snps_have_p_values>=m) {
				double[][] selected_snp_in_this_region = new double[m][];
				int index=0;
				while(used_snps<m && index<pvalues.length) {
					String tag_locations = p_location.get(Double.toString(pvalues[index]));
					System.out.println("#"+tag_locations);
					if(tag_locations.contains("_")) {
						String[] tag_locationStrings = tag_locations.split("_");	
						if(used_snps+tag_locationStrings.length<=m) {
							//1.All snps in this tag locations need output
							for(int k=0;k<tag_locationStrings.length;k++) {
								//System.out.print(tag_locationStrings[k]);
								selected_snp_in_this_region[used_snps] = genotype.load_one_variant_by_location(chr-1, Integer.parseInt(tag_locationStrings[k]));
								used_snps++;}
						}else {
							//2.Some snps in this tag locations need output
							int cur_used_snps=used_snps;
							for(int k=0;k<m-cur_used_snps;k++) {
								//System.out.print(tag_locationStrings[k]);
								selected_snp_in_this_region[used_snps] = genotype.load_one_variant_by_location(chr-1, Integer.parseInt(tag_locationStrings[k]));
								used_snps++;}
						}
					}else {
						selected_snp_in_this_region[used_snps] = genotype.load_one_variant_by_location(chr-1, Integer.parseInt(tag_locations));
						used_snps++;
					}
					index++;
				}
				return selected_snp_in_this_region;		
			}else { //the number of snps with p.value is less than m
				double[][] selected_snp_in_this_region = new double[snps_have_p_values][];
				for(int p_index=0;p_index<pvalues.length;p_index++) {
					String tag_locations = p_location.get(Double.toString(pvalues[p_index]));
					System.out.println("#"+tag_locations);
					if(tag_locations.contains("_")) {
						String[] tag_locationStrings = tag_locations.split("_");	
						if(used_snps+tag_locationStrings.length<=snps_have_p_values) {
							//1.All snps in this tag locations need output
							for(int k=0;k<tag_locationStrings.length;k++) {
								//System.out.print(tag_locationStrings[k]);
								selected_snp_in_this_region[used_snps] = genotype.load_one_variant_by_location(chr-1, Integer.parseInt(tag_locationStrings[k]));
								used_snps++;}
						}else {
							//2.Some snps in this tag locations need output
							int cur_used_snps=used_snps;
							for(int k=0;k<snps_have_p_values-cur_used_snps;k++) {
								//System.out.print(tag_locationStrings[k]);
								selected_snp_in_this_region[used_snps] = genotype.load_one_variant_by_location(chr-1, Integer.parseInt(tag_locationStrings[k]));
								used_snps++;}
						}
					}else {
						selected_snp_in_this_region[used_snps] = genotype.load_one_variant_by_location(chr-1, Integer.parseInt(tag_locations));
						used_snps++;
					}
				}
				return selected_snp_in_this_region;
			}				
		}else if(snp_in_this_region.length <m && snp_in_this_region.length >0) {
			return snp_in_this_region;
		}else{
			return null;
		}
	}
	
	public void analysis_specified_compounds(String mlr_output_file, int[][][] compounds, boolean plot, String emmax_res, int m){
		try{
			System.out.println(phenotype.phe_id);		
			BufferedWriter bw=new BufferedWriter(new FileWriter(mlr_output_file));
			bw.write("#SampleSize="+this.sample_size+ "; REML="+this.Likelihood_null+"; " +
					"h0-heritability="+this.sigma_g_null/(this.sigma_e_null+this.sigma_g_null) +"; delta="+this.delta_null+"\n");
			bw.write("compound(chr:start:end_chr:start_end), best_ml-NULL, p-value, variance_local, variance_global, variance_e\n");	
			final double[] all_one=new double[this.sample_size];
			Arrays.fill(all_one, 1);
			final double[] all_one_gt=this.transformation_matrix.operate(all_one); //Returns the result of multiplying this by the vector v.
			int total=0;
			int[] low_counts={0};			
			for(int compound_index=0;compound_index<compounds.length;compound_index++){
				String compound_name="";
				total++;
				//int start=1, end=start+tiling_win_size-1;
				ArrayList<double[][]> data_buffer=new ArrayList<double[][]>();
				ArrayList<double[][]> data_selected_buffer=new ArrayList<double[][]>();
				int num_snps_in_compound=0;
				for(int region_index=0;region_index<compounds[compound_index].length;region_index++){					
					int chr=compounds[compound_index][region_index][0];
					int start=compounds[compound_index][region_index][1];
					int end=compounds[compound_index][region_index][2];
					double[][] snp_in_this_region=genotype.load_variants_in_region(chr-1, start, end);
					data_buffer.add(snp_in_this_region);
					if(region_index==compounds[compound_index].length-1)
						compound_name=compound_name+chr+":"+start+":"+end;
					else compound_name=compound_name+chr+":"+start+":"+end+"_";
					num_snps_in_compound+=snp_in_this_region.length;
					//Select SNPs based on the p.value
					//For location from start to end, search one by one, requires results in both genotype and res file
					double[][]	snp_selected_in_this_region = select_SNPs_base_emmax_res(snp_in_this_region, m, chr,start,end,emmax_res);
					if(snp_selected_in_this_region!=null) {data_selected_buffer.add(snp_selected_in_this_region);}
				}//At the end of this for loop, we finish loading SNPs in all regions
				//If all regions have one SNPs rest, we construct the new genotype matrix; Else, we move to the next compound
				System.out.println("data buffer: "+data_buffer.size());
				System.out.println("data selected buffer: "+data_selected_buffer.size());
				double[][] combined_data_p1=new double[num_snps_in_compound][];
				int combined_index=0;
				for(int region_index=0;region_index<compounds[compound_index].length;region_index++){	
					double[][] snp_in_this_region=data_buffer.get(region_index);
					for(int k=0;k<snp_in_this_region.length;k++){
						combined_data_p1[combined_index++]=snp_in_this_region[k];
					}
				}
				ArrayList<double[]> component_productArrayList = new ArrayList<double[]>();
				if(data_selected_buffer.size()==2) {
					double[][] snp_in_1st_selected_region=data_selected_buffer.get(0);
					double[][] snp_in_2nd_selected_region=data_selected_buffer.get(1);	
					for(int i=0; i<snp_in_1st_selected_region.length;i++) {
						for(int j=0; j<snp_in_2nd_selected_region.length;j++) {
							double[] component_prodcut =  new double[snp_in_1st_selected_region[0].length];
							boolean valid_snp = false;
							for(int k=0;k<snp_in_1st_selected_region[0].length;k++) {
								try {
									component_prodcut[k]=snp_in_1st_selected_region[i][k]*snp_in_2nd_selected_region[j][k];
								}catch(NullPointerException npe) {
									System.out.println(compound_name+","+emmax_res+", 1: "+snp_in_1st_selected_region.length+", 2: "+snp_in_2nd_selected_region.length+", k:"+k);
								}
								if(component_prodcut[k]!=0) {
									valid_snp = true;
								}
							}
							if(valid_snp) {
								component_productArrayList.add(component_prodcut);
							}
						}
					}
				}
				double[][] combined_data_final;
				combined_data_final=new double[num_snps_in_compound+component_productArrayList.size()][];
				for(int i=0;i<num_snps_in_compound;i++) {
					combined_data_final[i]=combined_data_p1[i];
				}
				for(int j=0;j<component_productArrayList.size();j++) {
					combined_data_final[num_snps_in_compound+j]=component_productArrayList.get(j);
				}					
				double[][] data_matching_phenotype=this.select_subset_of_data(combined_data_final);
				if(data_matching_phenotype.length<2){
					//System.out.println("NO SNP: "+chr+"/"+start);
					continue;
				}
//				for(int i=0;i<data_matching_phenotype.length;i++) {
//					for(int j=0;j<data_matching_phenotype[0].length;j++) {
//						bw.write(data_matching_phenotype[i][j]+",");
//					}bw.write("\n");
//				}
				run_a_compound(data_matching_phenotype, all_one_gt, low_counts, bw, compound_name);
			}bw.close();
			System.out.println("Total #compounds="+total+"; #Low-rank regions="+low_counts[0]);
			//if(plot)LocalKinshipAnalyzer.make_plot_one_phen(mlr_output_file);
		}catch(Exception e){e.printStackTrace();}
	}
	
	public void run_a_compound(double[][] data_matching_phenotype, double[] all_one_gt, int[] low_full_counts, BufferedWriter bw, 
			String compound_name){
		//System.out.println(chr+"/"+start);
		try{
			FaSTLMM.normalize_terms_withSqrtNumVarDivided(data_matching_phenotype);
			final double[] all_one=new double[this.sample_size];
			Arrays.fill(all_one, 1);
			double[][] transformed_data=FaSTLMM.transform_X_by_Ut(transformation_matrix, data_matching_phenotype);
			double[][] X_for_null=new double[1][];
			X_for_null[0]=all_one_gt.clone();
			
			VarComp_Result result_null=FaSTLMM_FullRank.fullRankSolver_FixedEffects_ML(X_for_null, this.Y_global_transformed.clone(), 
					all_one, 0);
			double ml_null=result_null.ml;
			
			
			if(data_matching_phenotype.length<sample_size){ //low rank
				//System.out.println("Low-Rank chr"+(chr+1)+" win "+start);
				low_full_counts[0]++;
				double[] Y_new_original=this.Y_global_transformed;
				//FaSTLMM.normalize_terms_withSqrtNumVarDivided(transformed_data); NO NEED TO NORMALIZE AGAIN!!!
				SingularValueDecomposition local_svd=new SingularValueDecomposition(
						(new Array2DRowRealMatrix(transformed_data)).transpose());	
				
				double[] S=local_svd.getSingularValues();
				for(int k=0;k<S.length;k++)S[k]=(S[k]*S[k]);
				RealMatrix local_U1t= local_svd.getUT();
				double[] transformed_Y=local_U1t.operate(Y_new_original);
				RealMatrix local_IminusU1U1t=local_svd.getU().multiply(local_U1t).scalarMultiply(-1);
				for(int i=0;i<this.sample_size;i++)local_IminusU1U1t.setEntry(i, i, 1+local_IminusU1U1t.getEntry(i, i));
				double[] Y_transformed2=local_IminusU1U1t.operate(Y_new_original);
				double[][] X_new_original=new double[1][];
				X_new_original[0]=all_one_gt;
				//double[][] transformed_X=new double[1][];
				//transformed_X[0]=local_U1t.operate(all_one);						
				VarComp_Result result=FaSTLMM_LowRank.lowRankSolver_FixedEffects(X_new_original, Y_new_original, 
						transformed_Y,Y_transformed2, S, ml_null, local_U1t, local_IminusU1U1t);								
				if(result!=null){
					double total_variance=result.sigma_g+result.sigma_e;//*(1+this.delta_null);
					bw.write(compound_name + ", "+(result.ml-ml_null)+", "+result.pvalue+", "+
							result.sigma_g/total_variance+", "+(result.sigma_e/total_variance)*(1/(1+this.delta_null))+", "+ 
							(result.sigma_e/total_variance)*(this.delta_null/(1+this.delta_null))+"\n");
					//System.out.println("Succ finsihed chr"+(chr+1)+" win "+start);
				}
			}else if(data_matching_phenotype.length>=sample_size){	// full-rank						
				RealMatrix local_kinship= VariantsDouble.calculate_RRM_local_kinship_RealMatrix_noNomalization(transformed_data);
				//RealMatrix local_kinship=KinshipMatrix.re_scale_kinship_matrix(local_kinship0, true);
				double[] Y_new_original=this.Y_global_transformed;
				EigenDecomposition local_eigen= null;
				try{local_eigen= new EigenDecomposition((Array2DRowRealMatrix)local_kinship,0);}
				catch(Exception e){}
				if(local_eigen!=null){
					double[] S=local_eigen.getRealEigenvalues();
					RealMatrix local_kinship_Ut= local_eigen.getVT();
					double[] transformed_Y=local_kinship_Ut.operate(Y_new_original);
					double[][] transformed_X=new double[1][];
					transformed_X[0]=local_kinship_Ut.operate(all_one_gt);
					//double ml_null=myMathLib.StatFuncs.null_logL(transformed_Y,transformed_X[0]);
					//System.out.println(ml_null+","+this.Likelihood_null);
					VarComp_Result result=FaSTLMM_FullRank.fullRankSolver_FixedEffects_ML(transformed_X, transformed_Y, 
							S, ml_null);
					
					if(result!=null){
						double total_variance=result.sigma_g+result.sigma_e;//*(1+this.delta_null);
						bw.write(compound_name+", "+(result.ml-ml_null)+", "+result.pvalue+", "+
								result.sigma_g/total_variance+", "+(result.sigma_e/total_variance)*(1/(1+this.delta_null))+", "+ 
								(result.sigma_e/total_variance)*(this.delta_null/(1+this.delta_null))+"\n");
				}				
				//	System.out.println("Succ finsihed chr"+(chr+1)+" win "+start);
				}
			}	
		}catch(Exception e){e.printStackTrace();}
	}
}
