package real_data_analysis;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.Random;

//import org.apache.commons.io.FileSystemUtils;
//import org.apache.commons.lang.WordUtils;

import mixedmodel.AssociationResults;
import mixedmodel.EMMAX;
import mixedmodel.FaSTLMM_LocalKinship;
import mixedmodel.KinshipMatrix;
import mixedmodel.LocalKinshipAnalyzer;
import mixedmodel.MultiPhenotype;
import mixedmodel.Phenotype;
import mixedmodel.RareAnalyzerAggregate;
import mixedmodel.Regions;
import mixedmodel.VarianceComponent;
import mixedmodel.VariantsDouble;
import myMathLib.Test;
import myPlotLab.MyLineChart;

import phenotype_simulation.Simulator;


public class EntranceSim {
	
	public static double emmax_arabidopsis_threshold=0.019679907/1416064;
	public static double aggregate_arabidopsis_threshold=0.01212773/12929;
	public static double local_gene_arabidopsis_threshold=0.052253921/12929;
	public static double local_win50k_arabidopsis_threshold=0.058879085/4800;

	//public static double[] heritability={0.9,0.7,0.5,0.3,0.1};	
	//public static double[] concentration_rare={0.5,0.4,0.3,0.2,0.1};
	public static double[] heritability={0.1,0.3,0.5,0.7,0.9};
	public static double[] concentration_rare={0.1,0.2,0.3,0.4,0.5};
	public static double[] rare_maf={0.01,0.02,0.03,0.04,0.05};	//for old plots
//	public static double[] rare_maf={0.01};	
	public static double[] common_maf={0.05,0.10};
//	public static double[] num_of_causal_common={1,2,3,5,8,11,15,20,25,30};
	public static double[] num_of_causal_common={1,2,3,4,5};
	
	public static void have_a_look(){
		String data_folder="/Volumes/Projects/DATA/arabidopsis_genomes/Genotype250K/";
		String genotype=data_folder+"call_method_54.tair9.csv.num.hdf5";
		String k_file1=data_folder+"call_method_54.tair9.csv.num.kinship.RRMkinship.RRM";
		String k_file2=data_folder+"call_method_54.tair9.csv.num.kinship.IBSkinship.rescaled.ibs";
		
		String simulation_folder="/Volumes/Projects/Local-kinship/simulations/";
		String phenotype_file=simulation_folder+"additive.tsv";
		String vo_file=simulation_folder+"h.vo.e100.80.csv";
		
		double min_maf=0.4;
		int chr=1;
		int start=15657492-50000;
		int end=15660261+50000;
		
		
		String model1=Simulator.additiveInfinitesimal;
		String model2=Simulator.additive_exp_pow;
		int num_of_causal=100;
		
		int num_replicates=10;
		double exp_ratio=0.6;
		
		Simulator sim=new Simulator(genotype, 0.80, min_maf, 0.5);
//		sim.additive(model2, num_of_causal, exp_ratio, num_replicates, phenotype_file);
		sim.additive(model1, num_of_causal, -1, num_replicates, false, phenotype_file);
//		sim.heterogeneity(Simulator.heterogeneityRecessive, num_of_causal, num_replicates, 1, phenotype_file);
//		sim.canalization(num_replicates, 20, phenotype_file);
		
		KinshipMatrix kinship1=new KinshipMatrix(k_file1);
		KinshipMatrix kinship2=new KinshipMatrix(k_file2);
		MultiPhenotype phenotypes=new MultiPhenotype(phenotype_file);
		
		VarianceComponent.VO_sinle_kinship_fastlmm(kinship1, kinship2, phenotypes, vo_file);
		VarianceComponent.VO_sinle_kinship_emma(kinship1, kinship2, phenotypes, new VariantsDouble(genotype), vo_file+".emma");
	}
	
	public static void sim_local_trans(String input_genotype, String true_signal_file, String temp_pheno_file, 
			int num_replicates, double heritability, double min_maf, double max_maf, int total_num_of_causal,
			int num_allele_needed_to_have_efffet, int win_size, String tmep_output_folder, String global_kinship_file,
			String regions_file,
			int chr, int start_loc, int end_loc){
		try{
			Simulator sim=new Simulator(input_genotype, heritability, min_maf, max_maf, chr, start_loc, end_loc);
			sim.heterogeneity(Simulator.heterogeneityEpistasis, total_num_of_causal, num_replicates, 
					num_allele_needed_to_have_efffet, false, temp_pheno_file);			
			MultiPhenotype phenotypeS=new MultiPhenotype(temp_pheno_file);
			for(int phe_index=0;phe_index<phenotypeS.num_of_pheno;phe_index++){
				String out_phe_file=tmep_output_folder+"Local_VO."+phe_index+"."+phenotypeS.phenotypes[phe_index].phe_id+".w"+
						win_size+".csv";
				String out_phe_file_gene=tmep_output_folder+"Local_VO."+phe_index+"."+phenotypeS.phenotypes[phe_index].phe_id+".gene"+".csv";
				FaSTLMM_LocalKinship calculator=new FaSTLMM_LocalKinship(phenotypeS.phenotypes[phe_index], input_genotype, 
						global_kinship_file);
				int[][][] regions=(new Regions(regions_file, calculator.genotype)).region_coordinates;
				calculator.analysis_specified_region(out_phe_file_gene, regions, true);
				calculator.analysis_tiling_win(out_phe_file, win_size, true);
				
			}
		}catch(Exception e){e.printStackTrace();}
	}
	
	public static void sim_three(String input_genotype, String true_signal_file, String temp_pheno_file, 
			int num_replicates, double heritability, double min_maf, double max_maf, int total_num_of_causal,
			int num_allele_needed_to_have_efffet, int win_size, String tmep_output_folder, String global_kinship_file,
			String regions_file,
			int chr, int start_loc, int end_loc){
		try{
			Simulator sim=new Simulator(input_genotype, heritability, min_maf, max_maf, chr, start_loc, end_loc);
			sim.heterogeneity(Simulator.heterogeneityEpistasis, total_num_of_causal, num_replicates, 
					num_allele_needed_to_have_efffet, false, temp_pheno_file);			
			MultiPhenotype phenotypeS=new MultiPhenotype(temp_pheno_file);
			for(int phe_index=0;phe_index<phenotypeS.num_of_pheno;phe_index++){
				String out_phe_file=tmep_output_folder+"Local_VO."+phe_index+"."+phenotypeS.phenotypes[phe_index].phe_id+".w"+
						win_size+".csv";  
				String out_phe_file_gene=tmep_output_folder+"Local_VO."+phe_index+"."+phenotypeS.phenotypes[phe_index].phe_id+".gene"+".csv";
//				FaSTLMM_LocalKinship calculator=new FaSTLMM_LocalKinship(phenotypeS.phenotypes[phe_index], input_genotype, 
//						global_kinship_file);
//				int[][][] regions=(new Regions(regions_file, calculator.genotype)).region_coordinates;
//				calculator.analysis_specified_region(out_phe_file_gene, regions, true);
//				calculator.analysis_tiling_win(out_phe_file, win_size, true);
				RareAnalyzerAggregate rare=new RareAnalyzerAggregate(input_genotype, regions_file, phenotypeS.phenotypes[phe_index],
						global_kinship_file, max_maf);
				rare.rare_association(0.01, tmep_output_folder, true);
				
			}
		}catch(Exception e){e.printStackTrace();}
	}
	
//	public static void sim_local_grid(String input_genotype, String true_signal_file, String temp_pheno_file, 
//			int num_replicates, double heritability, double min_maf, double max_maf, int total_num_of_causal,
//			int num_allele_needed_to_have_efffet, int win_size, String tmep_output_folder, String global_kinship_file, String local_kinship_files_folder,
//			int chr, int start_loc, int end_loc){
//		try{
//			Simulator sim=new Simulator(input_genotype, heritability, min_maf, max_maf, chr, start_loc, end_loc);
//			sim.heterogeneity(Simulator.heterogeneityEpistasis, total_num_of_causal, num_replicates, 
//					num_allele_needed_to_have_efffet, temp_pheno_file);			
//			MultiPhenotype phenotypeS=new MultiPhenotype(temp_pheno_file);
//			LocalKinshipAnalyzer local_k=new LocalKinshipAnalyzer(input_genotype, win_size, null);
//			for(int phe_index=0;phe_index<phenotypeS.num_of_pheno;phe_index++){
//				String out_phe_file=tmep_output_folder+"Local_VO."+phe_index+"."+phenotypeS.phenotypes[phe_index].phe_id+".w"+
//						win_size+".csv";
//				
//				local_k.local_VO(phenotypeS.phenotypes[phe_index], input_genotype, global_kinship_file, local_kinship_files_folder,
//						out_phe_file, 0.1, true, "grid");
//			}
//		}catch(Exception e){e.printStackTrace();}
//	}
	
	/*
	 * MAIN Rare function for power
	 */
	public static void RARE_power_coparion(String working_folder, String input_genotype, String regions_file, String global_kinship_file,
			int sim_win, int analysis_tiling_win_size, int round_index){		
		VariantsDouble genotype=new VariantsDouble(input_genotype);
		Regions regions=new Regions(regions_file, genotype);
		Random rg = new Random();
		//
		try{			
			for(int h_index=0;h_index<heritability.length;h_index++){
				for(int maf_index=0;maf_index<rare_maf.length;maf_index++){
					double max_maf=rare_maf[maf_index];
					double min_maf=0;
					if(maf_index!=0) min_maf=rare_maf[maf_index-1];					
					for(int concentration_index=0;concentration_index<concentration_rare.length;concentration_index++){
						int chr=(int)(rg.nextDouble()*genotype.num_chrs);
						if(chr==genotype.num_chrs)chr--;
						int region_index=(int)(rg.nextDouble()*regions.region_coordinates[chr].length);
						if(region_index==regions.region_coordinates[chr].length)region_index--;
						int start_loc=regions.region_coordinates[chr][region_index][0]-sim_win;
						int end_loc=regions.region_coordinates[chr][region_index][1]+sim_win;
						
						String simulated_phenotype_file_folder=working_folder+"/SIM.h"+heritability[h_index]+".min"+min_maf+".max"+max_maf+
								".con"+concentration_rare[concentration_index]+".w"+sim_win+
								".chr"+(chr+1)+".start"+start_loc+".end"+end_loc+".round"+round_index+"/";
						File the_folder=new File(simulated_phenotype_file_folder);
						if(!the_folder.exists())the_folder.mkdir();
						Simulator sim=new Simulator(input_genotype, heritability[h_index], min_maf, max_maf, chr, start_loc, end_loc);
						String phenotype_file_simulated=simulated_phenotype_file_folder+"simulated_phenotype.tsv";
						sim.heterogeneity2(Simulator.heterogeneityEpistasis, concentration_rare[concentration_index], 
								1, 4, false, phenotype_file_simulated);
						sim.additive2(Simulator.additiveBidirection, concentration_rare[concentration_index], -1, 1, false, phenotype_file_simulated);
						sim.liabilityThreshold2(Simulator.liabilityThreshold, concentration_rare[concentration_index], 1, false, phenotype_file_simulated);
						
						BufferedWriter bw=new BufferedWriter(new FileWriter(simulated_phenotype_file_folder+"/summary.txt"));
						
						if(phenotype_check_passed(phenotype_file_simulated)){
							running_single_round_rare_QUICK(regions_file, analysis_tiling_win_size, input_genotype,  
									phenotype_file_simulated, simulated_phenotype_file_folder, global_kinship_file, 
									bw, chr, start_loc+sim_win,end_loc-sim_win);
							
							running_single_round_rare_QUICK(regions_file, analysis_tiling_win_size, input_genotype,  
									phenotype_file_simulated, simulated_phenotype_file_folder, global_kinship_file, 
									bw, chr, start_loc,end_loc);
							
							running_single_round_rare_QUICK(regions_file, analysis_tiling_win_size, input_genotype,  
									phenotype_file_simulated, simulated_phenotype_file_folder, global_kinship_file, 
									bw, chr, start_loc+sim_win-analysis_tiling_win_size/2,end_loc-sim_win+analysis_tiling_win_size/2);
						}
						bw.close();
						
					}
				}
			}			
		}catch(Exception e){e.printStackTrace();}
	}
	
	/*
	 * MAIN Common function for power
	 */
	public static void COMMON_power_coparion(String working_folder, String input_genotype, String regions_file, String global_kinship_file,
			int sim_win, int analysis_tiling_win_size, int round_index, String analysis_type){		
		VariantsDouble genotype=new VariantsDouble(input_genotype);
		Regions regions=new Regions(regions_file, genotype);
		Random rg = new Random();
		//
		try{			
			for(int h_index=0;h_index<heritability.length;h_index++){				
					double max_maf=common_maf[1];
					double min_maf=common_maf[0];				
					for(int num_index=0;num_index<num_of_causal_common.length;num_index++){
						int chr=(int)(rg.nextDouble()*genotype.num_chrs);
						if(chr==genotype.num_chrs)chr--;
						int region_index=(int)(rg.nextDouble()*regions.region_coordinates[chr].length);
						if(region_index==regions.region_coordinates[chr].length)region_index--;
						int start_loc=regions.region_coordinates[chr][region_index][0]-sim_win;
						int end_loc=regions.region_coordinates[chr][region_index][1]+sim_win;
						
						String simulated_phenotype_file_folder=working_folder+"/SIM.h"+heritability[h_index]+".min"+min_maf+".max"+max_maf+
								".num"+num_of_causal_common[num_index]+".w"+sim_win+
								".chr"+(chr+1)+".start"+start_loc+".end"+end_loc+".round"+round_index+"/";
						File the_folder=new File(simulated_phenotype_file_folder);
						if(!the_folder.exists())the_folder.mkdir();
						Simulator sim=new Simulator(input_genotype, heritability[h_index], min_maf, max_maf, 
								chr, start_loc, end_loc);
						String phenotype_file_simulated=simulated_phenotype_file_folder+"simulated_phenotype.tsv";
						if(analysis_type.equals("dorminant"))
							sim.heterogeneity(Simulator.heterogeneityDorminant, (int)num_of_causal_common[num_index], 
								1, 1, false, phenotype_file_simulated);
						else if(analysis_type.equals("epistasis")) 
							sim.heterogeneity(Simulator.heterogeneityEpistasis, (int)num_of_causal_common[num_index], 
								1, 4, false, phenotype_file_simulated);
						else if(analysis_type.equals("additive")) 
							sim.additive(Simulator.additiveBidirection, (int)num_of_causal_common[num_index], -1, 1, false, phenotype_file_simulated);
						else if(analysis_type.equals("liability")) 
							sim.liabilityThreshold(Simulator.liabilityThreshold, (int)num_of_causal_common[num_index], 1, false, phenotype_file_simulated);
						
						//BufferedWriter bw=new BufferedWriter(new FileWriter(simulated_phenotype_file_folder+"/summary.txt"));
						
						if(phenotype_check_passed(phenotype_file_simulated)){
							running_single_round_common_QUICK(regions_file, analysis_tiling_win_size, input_genotype,  
									phenotype_file_simulated, simulated_phenotype_file_folder, global_kinship_file, 
									simulated_phenotype_file_folder+"/summary.txt", chr, start_loc,end_loc);
						}
						//bw.close();
						
					}
				
			}			
		}catch(Exception e){e.printStackTrace();}
	}
	
	/*
	 * MAIN Rare function for power, adding background genetic effects and multiple regressions.
	 */
	public static void RARE_power_coparion_bk(String working_folder, String input_genotype, String regions_file, String global_kinship_file,
			int sim_win, int analysis_tiling_win_size, int round_index){		
		VariantsDouble genotype=new VariantsDouble(input_genotype);
		Regions regions=new Regions(regions_file, genotype);
		Random rg = new Random();
		//
		try{			
			for(int h_index=0;h_index<heritability.length;h_index++){
				for(int maf_index=0;maf_index<rare_maf.length;maf_index++){
					double max_maf=rare_maf[maf_index];
					double min_maf=0;
					if(maf_index!=0) min_maf=rare_maf[maf_index-1];					
					for(int concentration_index=0;concentration_index<concentration_rare.length;concentration_index++){
						int chr=(int)(rg.nextDouble()*genotype.num_chrs);
						if(chr==genotype.num_chrs)chr--;
						int region_index=(int)(rg.nextDouble()*regions.region_coordinates[chr].length);
						if(region_index==regions.region_coordinates[chr].length)region_index--;
						int start_loc=regions.region_coordinates[chr][region_index][0]-sim_win;
						int end_loc=regions.region_coordinates[chr][region_index][1]+sim_win;
						
						String simulated_phenotype_file_folder=working_folder+"/SIM.h"+heritability[h_index]+".min"+min_maf+".max"+max_maf+
								".con"+concentration_rare[concentration_index]+".w"+sim_win+
								".chr"+(chr+1)+".start"+start_loc+".end"+end_loc+".round"+round_index+"/";
						File the_folder=new File(simulated_phenotype_file_folder);
						if(!the_folder.exists())the_folder.mkdir();
						Simulator sim=new Simulator(input_genotype, heritability[h_index], min_maf, max_maf, chr, start_loc, end_loc);
						String phenotype_file_simulated=simulated_phenotype_file_folder+"simulated_phenotype.tsv";
						sim.heterogeneity2(Simulator.heterogeneityEpistasis, concentration_rare[concentration_index], 
								1, 4, false, phenotype_file_simulated);
						sim.additive2(Simulator.additiveBidirection, concentration_rare[concentration_index], -1, 1, false, phenotype_file_simulated);
						sim.liabilityThreshold2(Simulator.liabilityThreshold, concentration_rare[concentration_index], 1, false, phenotype_file_simulated);
						
						BufferedWriter bw=new BufferedWriter(new FileWriter(simulated_phenotype_file_folder+"/summary.txt"));
						
						if(phenotype_check_passed(phenotype_file_simulated)){
							running_single_round_rare_QUICK(regions_file, analysis_tiling_win_size, input_genotype,  
									phenotype_file_simulated, simulated_phenotype_file_folder, global_kinship_file, 
									bw, chr, start_loc+sim_win,end_loc-sim_win);
							
							running_single_round_rare_QUICK(regions_file, analysis_tiling_win_size, input_genotype,  
									phenotype_file_simulated, simulated_phenotype_file_folder, global_kinship_file, 
									bw, chr, start_loc,end_loc);
							
							running_single_round_rare_QUICK(regions_file, analysis_tiling_win_size, input_genotype,  
									phenotype_file_simulated, simulated_phenotype_file_folder, global_kinship_file, 
									bw, chr, start_loc+sim_win-analysis_tiling_win_size/2,end_loc-sim_win+analysis_tiling_win_size/2);
						}
						bw.close();
						
					}
				}
			}			
		}catch(Exception e){e.printStackTrace();}
	}
	
	/*
	 * MAIN Common function for power, adding background genetic effects and multiple regressions
	 */
	public static void COMMON_power_coparion_bk(String working_folder, String input_genotype, String regions_file, String global_kinship_file,
			int sim_win, int analysis_tiling_win_size, int round_index){		
		VariantsDouble genotype=new VariantsDouble(input_genotype);
		Regions regions=new Regions(regions_file, genotype);
		Random rg = new Random();
		//
		try{			
			for(int h_index=0;h_index<heritability.length;h_index++){				
				double max_maf=common_maf[1];
				double min_maf=common_maf[0];				
				for(int num_index=0;num_index<num_of_causal_common.length;num_index++){
					int chr=(int)(rg.nextDouble()*genotype.num_chrs);
					if(chr==genotype.num_chrs)chr--;
					int region_index=(int)(rg.nextDouble()*regions.region_coordinates[chr].length);
					if(region_index==regions.region_coordinates[chr].length)region_index--;
					int start_loc=regions.region_coordinates[chr][region_index][0]-sim_win;
					int end_loc=regions.region_coordinates[chr][region_index][1]+sim_win;

					String simulated_phenotype_file_folder=working_folder+"/SIM.h"+heritability[h_index]+".min"+min_maf+".max"+max_maf+
							".num"+num_of_causal_common[num_index]+".w"+sim_win+
							".chr"+(chr+1)+".start"+start_loc+".end"+end_loc+".round"+round_index+"/";
					File the_folder=new File(simulated_phenotype_file_folder);
					if(!the_folder.exists())the_folder.mkdir();
					Simulator sim=new Simulator(input_genotype, heritability[h_index], min_maf, max_maf, chr, start_loc, end_loc);
					String phenotype_file_simulated=simulated_phenotype_file_folder+"simulated_phenotype.tsv";
					//sim.heterogeneity(Simulator.heterogeneityEpistasis, (int)num_of_causal_common[num_index], 
					//								1, 1, phenotype_file_simulated);
					//sim.additive(Simulator.additiveBidirection, (int)num_of_causal_common[num_index], -1, 1, false, phenotype_file_simulated);
					sim.liabilityThreshold(Simulator.liabilityThreshold, (int)num_of_causal_common[num_index], 1, true, phenotype_file_simulated);

					BufferedWriter bw=new BufferedWriter(new FileWriter(simulated_phenotype_file_folder+"/summary.txt"));

					if(phenotype_check_passed(phenotype_file_simulated)){
						running_single_round_common2_QUICK(regions_file, analysis_tiling_win_size, input_genotype,  
								phenotype_file_simulated, simulated_phenotype_file_folder, global_kinship_file, 
								bw, chr, start_loc,end_loc);
					}
					bw.close();

				}
				
			}			
		}catch(Exception e){e.printStackTrace();}
	}
	
	/*
	 * type I error for both regions and windows:
	 * sort the p-values and take the top 0.05. then divide the cutoff with number of tests.  
	 * 
	 * type I error for four method:
	 * 
	 * Rare, EMMAX, Local-window, local-region
	 */
	public static void calculate_typeIerror(String regions_file, int win_size, String genotype_file, double typeIerror, 
			String random_phenotype_file, int round, String results_output_folder, String global_kinship_file){
		try{
			VariantsDouble genotype=new VariantsDouble(genotype_file);
			for(int round_index=0;round_index<round;round_index++){
				Simulator.generate_random_phenotype(random_phenotype_file, genotype.sample_ids);
				Phenotype phenotype=new MultiPhenotype(random_phenotype_file).phenotypes[0];
				FaSTLMM_LocalKinship calculator=new FaSTLMM_LocalKinship(phenotype, genotype_file, 
						global_kinship_file);
				int[][][] regions=(new Regions(regions_file, calculator.genotype)).region_coordinates;
				// emmax, region
				EMMAX.emmax_analysis_regions(genotype_file, random_phenotype_file, global_kinship_file, results_output_folder, 
						1000000, 10, 0, 0.05, true, regions);
				//local
				String out_phe_file=results_output_folder+"Local_VO."+round_index+"."+phenotype.phe_id+".w"+
						win_size+".csv";
				String out_phe_file_gene=results_output_folder+"Local_VO."+round_index+"."+phenotype.phe_id+".gene"+".csv";								
				calculator.analysis_specified_region(out_phe_file_gene, regions, true);
				calculator.analysis_tiling_win(out_phe_file, win_size, true);
				//rare
				RareAnalyzerAggregate rare=new RareAnalyzerAggregate(genotype_file, regions_file, phenotype,
						global_kinship_file, 0.01);
				rare.rare_association(0.01, results_output_folder, true);				
			}			
		}catch(Exception e){e.printStackTrace();}
	}
	
	/*
	 * type I error for both regions and windows:
	 * sort the p-values and take the top 0.05. then divide the cutoff with number of tests.  
	 * 
	 * type I error for four method:
	 * 
	 * Rare, EMMAX, Local-window, local-region
	 */
	public static void running_single_round(String regions_file, int win_size, String genotype_file, 
			String simulated_phenotype_file, String results_output_folder, String global_kinship_file,
			int chr, int start, int end){
		try{
			//VariantsDouble genotype=new VariantsDouble(genotype_file);
			Phenotype phenotype=new MultiPhenotype(simulated_phenotype_file).phenotypes[0];
			FaSTLMM_LocalKinship calculator=new FaSTLMM_LocalKinship(phenotype, genotype_file, 
					global_kinship_file);
			int[][][] regions=(new Regions(regions_file, calculator.genotype)).region_coordinates;
			
			//aggregate
			String aggregate_outfile=results_output_folder+"/Aggregate."+phenotype.phe_id+".csv";
			RareAnalyzerAggregate rare=new RareAnalyzerAggregate(genotype_file, regions_file, phenotype,
					global_kinship_file, 0.01);
			rare.rare_association(0.01, results_output_folder, true);	
			// emmax, region
			String emmax_outputfile =results_output_folder+"/EMMAX.0_"+phenotype.phe_id+".top";
			EMMAX.emmax_analysis_regions(genotype_file, simulated_phenotype_file, global_kinship_file, results_output_folder,  
					1000000, 10, 0, 0.05, true, regions);
			//local
			String local_win_outfile=results_output_folder+"/Local_VO."+phenotype.phe_id+".w"+
					win_size+".csv";
			String local_gene_outfile=results_output_folder+"/Local_VO."+phenotype.phe_id+".gene"+".csv";								
			calculator.analysis_specified_region(local_gene_outfile, regions, true);
			calculator.analysis_tiling_win(local_win_outfile, win_size, true);
						
			BufferedWriter bw=new BufferedWriter(new FileWriter(results_output_folder+"/summary.txt"));
			bw.write("EMMAX\t"+check_power(emmax_outputfile, chr, start, end, EntranceSim.emmax_arabidopsis_threshold)+"\t"+(-Math.log10(EntranceSim.emmax_arabidopsis_threshold))+"\n");
			bw.write("Local_win\t"+check_power(local_win_outfile, chr, start-win_size, end, EntranceSim.local_win50k_arabidopsis_threshold)+"\t"+(-Math.log10(EntranceSim.local_win50k_arabidopsis_threshold))+"\n");
			bw.write("Local_gene\t"+check_power(local_gene_outfile, chr, start, end, EntranceSim.local_gene_arabidopsis_threshold)+"\t"+(-Math.log10(EntranceSim.local_gene_arabidopsis_threshold))+"\n");
			bw.write("Aggregate\t"+check_power(aggregate_outfile, chr, start, end, EntranceSim.aggregate_arabidopsis_threshold)+"\t"+(-Math.log10(EntranceSim.aggregate_arabidopsis_threshold))+"\n");
			bw.close();
		}catch(Exception e){e.printStackTrace();}
	}
	
	/*
	 * type I error for both regions and windows:
	 * sort the p-values and take the top 0.05. then divide the cutoff with number of tests.  
	 * 
	 * type I error for four method:
	 * 
	 * Rare, EMMAX, Local-window, local-region for specific regions ONLY. QUICK!
	 */
	public static void running_single_round_rare_QUICK(String regions_file, int analysis_win_size, String genotype_file, 
			String simulated_phenotype_file, String results_output_folder, String global_kinship_file,
			BufferedWriter bw, 
			int chr, int start, int end){
		try{
			//VariantsDouble genotype=new VariantsDouble(genotype_file);
			Phenotype phenotype=new MultiPhenotype(simulated_phenotype_file).phenotypes[0];
			FaSTLMM_LocalKinship calculator=new FaSTLMM_LocalKinship(phenotype, genotype_file, 
					global_kinship_file);
			int[][][] regions=(new Regions(regions_file, calculator.genotype)).region_coordinates;
			int[][][] region_only1=new int[regions.length][][];
			String[][] region_name=new String[regions.length][];
			for(int chr_index=0;chr_index<regions.length;chr_index++){
				if(chr_index==chr){
					region_only1[chr]=new int[1][2];
					region_only1[chr][0][0]=start;
					region_only1[chr][0][1]=end;
					region_name[chr]=new String[1];
					region_name[chr][0]="The_Gene";
				}else{
					region_only1[chr_index]=new int[0][];
					region_name[chr_index]=new String[0];
				}
			}
			//aggregate
			String aggregate_outfile=results_output_folder+"/Aggregate."+phenotype.phe_id+".csv";
			RareAnalyzerAggregate rare=new RareAnalyzerAggregate(genotype_file, region_only1, region_name, phenotype,
					global_kinship_file, 0.01);
			rare.rare_association(0.01, results_output_folder, false);	
//			// emmax, region
//			String emmax_outputfile =results_output_folder+"/EMMAX.0_"+phenotype.phe_id+".top";
//			EMMAX.emmax_analysis_regions(genotype_file, simulated_phenotype_file, global_kinship_file, results_output_folder,  
//					1000000, 10, 0, 0.05, false, region_only1);
			//local
			String local_win_outfile=results_output_folder+"/Local_VO."+phenotype.phe_id+".w"+
					analysis_win_size+".csv";
			String local_gene_outfile=results_output_folder+"/Local_VO."+phenotype.phe_id+".gene"+".csv";								
			calculator.analysis_specified_region(local_gene_outfile, region_only1, false);
//			calculator.analysis_tiling_win(local_win_outfile, analysis_win_size, true);
				
//			bw.write("EMMAX."+start+"."+end+"\t"+check_power(emmax_outputfile, chr, start, end, EntranceSim.emmax_arabidopsis_threshold)+"\t"+(-Math.log10(EntranceSim.emmax_arabidopsis_threshold))+"\n");
			bw.write("Local_gene."+start+"."+end+"\t"+check_power(local_gene_outfile, chr, start, end, EntranceSim.local_gene_arabidopsis_threshold)+"\t"+(-Math.log10(EntranceSim.local_gene_arabidopsis_threshold))+"\n");
			bw.write("Aggregate."+start+"."+end+"\t"+check_power(aggregate_outfile, chr, start, end, EntranceSim.aggregate_arabidopsis_threshold)+"\t"+(-Math.log10(EntranceSim.aggregate_arabidopsis_threshold))+"\n");
//			bw.close(); DO NOT CLOSE BW!!
			
		}catch(Exception e){e.printStackTrace();}
	}
	
	public static void running_single_round_common_QUICK(String regions_file, int analysis_win_size, String genotype_file, 
			String simulated_phenotype_file, String results_output_folder, String global_kinship_file,
			String output_summary_file, int chr, int start, int end){
		try{
			//VariantsDouble genotype=new VariantsDouble(genotype_file);
			Phenotype phenotype=new MultiPhenotype(simulated_phenotype_file).phenotypes[0];
			FaSTLMM_LocalKinship calculator=new FaSTLMM_LocalKinship(phenotype, genotype_file, 
					global_kinship_file);
			int[][][] regions=(new Regions(regions_file, calculator.genotype)).region_coordinates;
			int[][][] region_only1=new int[regions.length][][];
			String[][] region_name=new String[regions.length][];
			for(int chr_index=0;chr_index<regions.length;chr_index++){
				if(chr_index==chr){
					region_only1[chr]=new int[1][2];
					region_only1[chr][0][0]=start;
					region_only1[chr][0][1]=end;
					region_name[chr]=new String[1];
					region_name[chr][0]="The_Gene";
				}else{
					region_only1[chr_index]=new int[0][];
					region_name[chr_index]=new String[0];
				}
			}
			// emmax, region
			String emmax_outputfile =results_output_folder+"/EMMAX.0_"+phenotype.phe_id+".top";
			EMMAX.emmax_analysis_regions(genotype_file, simulated_phenotype_file, global_kinship_file, results_output_folder,  
					0.05, 10, 0, 0.05, false, region_only1);
			//local
			String local_gene_outfile=results_output_folder+"/Local_VO."+phenotype.phe_id+".gene"+".csv";								
			calculator.analysis_specified_region(local_gene_outfile, region_only1, false);
//			calculator.analysis_tiling_win(local_win_outfile, analysis_win_size, true);
			BufferedWriter bw=new BufferedWriter(new FileWriter(output_summary_file));	
			bw.write("EMMAX."+start+"."+end+"\t"+check_power_specific_region_only(emmax_outputfile, EntranceSim.emmax_arabidopsis_threshold)+"\t"+(-Math.log10(EntranceSim.emmax_arabidopsis_threshold))+"\n");
			bw.write("Local_gene."+start+"."+end+"\t"+check_power(local_gene_outfile, chr, start, end, EntranceSim.local_win50k_arabidopsis_threshold)+"\t"+(-Math.log10(EntranceSim.local_win50k_arabidopsis_threshold))+"\n");
			bw.close();
			
		}catch(Exception e){e.printStackTrace();}
	}
	
	public static void running_single_round_common2_QUICK(String regions_file, int analysis_win_size, String genotype_file, 
			String simulated_phenotype_file, String results_output_folder, String global_kinship_file,
			BufferedWriter bw, 
			int chr, int start, int end){
		try{
			//VariantsDouble genotype=new VariantsDouble(genotype_file);
			Phenotype phenotype=new MultiPhenotype(simulated_phenotype_file).phenotypes[0];
			FaSTLMM_LocalKinship calculator=new FaSTLMM_LocalKinship(phenotype, genotype_file, 
					global_kinship_file);
			int[][][] regions=(new Regions(regions_file, calculator.genotype)).region_coordinates;
			int[][][] region_only1=new int[regions.length][][];
			String[][] region_name=new String[regions.length][];
			for(int chr_index=0;chr_index<regions.length;chr_index++){
				if(chr_index==chr){
					region_only1[chr]=new int[1][2];
					region_only1[chr][0][0]=start;
					region_only1[chr][0][1]=end;
					region_name[chr]=new String[1];
					region_name[chr][0]="The_Gene";
				}else{
					region_only1[chr_index]=new int[0][];
					region_name[chr_index]=new String[0];
				}
			}
			// emmax, region
			String emmax_outputfile =results_output_folder+"/EMMAX.0_"+phenotype.phe_id+".top";
			EMMAX.emmax_analysis_regions(genotype_file, simulated_phenotype_file, global_kinship_file, results_output_folder,  
					1000000, 10, 0, 0.05, false, region_only1);
			EMMAX.run_stepwise_for_one_region(pvalue_threshold, emma, maf_plot_threshold, genotype_hdf5_file, round, chr, start_loc, end_loc, decomposed_array, emma_estimated)
			//local
			String local_gene_outfile=results_output_folder+"/Local_VO."+phenotype.phe_id+".gene"+".csv";								
			calculator.analysis_specified_region(local_gene_outfile, region_only1, false);
//			calculator.analysis_tiling_win(local_win_outfile, analysis_win_size, true);
				
			bw.write("EMMAX."+start+"."+end+"\t"+check_power_specific_region_only(emmax_outputfile, EntranceSim.emmax_arabidopsis_threshold)+"\t"+(-Math.log10(EntranceSim.emmax_arabidopsis_threshold))+"\n");
			bw.write("Local_gene."+start+"."+end+"\t"+check_power(local_gene_outfile, chr, start, end, EntranceSim.local_gene_arabidopsis_threshold)+"\t"+(-Math.log10(EntranceSim.local_gene_arabidopsis_threshold))+"\n");
			bw.write("Stepwise."+start+"."+end+"\t"+check_power_stepwise(local_gene_outfile, chr, start, end, EntranceSim.local_gene_arabidopsis_threshold)+"\t"+(-Math.log10(EntranceSim.local_gene_arabidopsis_threshold))+"\n");

//			bw.close(); DO NOT CLOSE BW!!
			
		}catch(Exception e){e.printStackTrace();}
	}
	
	public static double check_power(String pvalue_results, int chr, int start, int end, double pvalue_threshold){
		AssociationResults results=new AssociationResults(pvalue_results, pvalue_threshold);
		double most_significant=0;
		for(int i=0;i<results.num_of_var;i++){
			if(results.chr[i]==(chr+1) && results.location[i]>=start && results.location[i]<=end){
				double the_p=1000;
				if(results.pvalue[i]!=0)the_p=-Math.log10(results.pvalue[i]);
				if(the_p>most_significant)
					most_significant=the_p;
			}
		}
		return most_significant;
	}
	
	/*
	 * on the result file that only specific regions are tested.
	 */
	public static String check_power_specific_region_only(String pvalue_results, double pvalue_threshold){
		AssociationResults results=new AssociationResults(pvalue_results, pvalue_threshold);
		double most_significant=0;
		int most_sig_index=-1;
		for(int i=0;i<results.num_of_var;i++){			
			double the_p=1000;
			if(results.pvalue[i]!=0)the_p=-Math.log10(results.pvalue[i]);
			if(the_p>most_significant){
				most_significant=the_p;
				most_sig_index=i;
			}		
		}
		return most_significant+"/"+results.chr[most_sig_index]+"/"+results.location[most_sig_index];
	}
	
	public static boolean phenotype_check_passed(String phenotype_simulated){
		return true; //TODO: check whether this is really right!!!
		
//		MultiPhenotype pheno=new MultiPhenotype(phenotype_simulated);
//		double[] value=pheno.phenotypes[0].values;
//		double big=0;
//		for(int k=0;k<value.length;k++){
//			if(value[k]>0.5)big++;
//		}
//		System.out.println("==="+value.length+"/"+big+"====");
//		double ratio=big/value.length;
//		return (ratio>=0.1 && ratio<=0.9);
	}
	
	/*
	 * SIM.h0.9.min0.05.max0.1.con0.5.w10000.chr4.start2208379.end2231114.round42
	 * 
	 * public static double[] heritability={0.9,0.7,0.5,0.3,0.1};	
	public static double[] concentration_rare={0.5,0.4,0.3,0.2,0.1};
	public static double[] rare_maf={0.01,0.05,0.1};
	 */
	public static void make_plots_rare(String data_folder){
		double[][][] powers_small_local=new double[heritability.length][rare_maf.length][concentration_rare.length];
		double[][][] powers_exact_local=new double[heritability.length][rare_maf.length][concentration_rare.length];
		double[][][] powers_big_local=new double[heritability.length][rare_maf.length][concentration_rare.length];
		
		double[][][] powers_small_agg=new double[heritability.length][rare_maf.length][concentration_rare.length];
		double[][][] powers_exact_agg=new double[heritability.length][rare_maf.length][concentration_rare.length];
		double[][][] powers_big_agg=new double[heritability.length][rare_maf.length][concentration_rare.length];
		
		double[][][] counts=new double[heritability.length][rare_maf.length][concentration_rare.length];
		
		try{
			String plot_folder=data_folder+"_plot";
			File plot_folder_file=new File(plot_folder);
			if(!plot_folder_file.exists())plot_folder_file.mkdir();
			File root=new File(data_folder);
			File[] subfolders=root.listFiles();
			for(int k=0;k<subfolders.length;k++){
				if(k%1000==0)System.out.println(k/1000+"K");
				String name=subfolders[k].toString().split("/")[subfolders[k].toString().split("/").length-1];
				//System.out.println(k+":"+name);
				double h=Double.parseDouble(name.substring(5, 8));				
				int max_index_name=name.indexOf("max");
				int con_index_name=name.indexOf("con");
				double maf=Double.parseDouble(name.substring(max_index_name+3,con_index_name-1));				
				double con=Double.parseDouble(name.substring(con_index_name+3,con_index_name+6));
				//System.out.println(h+"/"+maf+"/"+con);
				int h_index=-1, maf_index=-1, con_index=-1;
				for(int p=0;p<heritability.length;p++)if(heritability[p]==h)h_index=p;
				for(int p=0;p<rare_maf.length;p++)if(rare_maf[p]==maf)maf_index=p;
				for(int p=0;p<concentration_rare.length;p++)if(concentration_rare[p]==con)con_index=p;
				if(h_index==-1||maf_index==-1||con_index==-1){System.out.println("=====return====="); return;}
				//System.out.println(h_index+"/"+maf_index+"/"+con_index);
				BufferedReader br=new BufferedReader(new FileReader(subfolders[k]+"/summary.txt"));
				String line=br.readLine();
				if(line!=null){
					counts[h_index][maf_index][con_index]++;
					if(!(Double.parseDouble(line.split("\t")[1])==0.0))
						powers_small_local[h_index][maf_index][con_index]++;
					if(!(Double.parseDouble(br.readLine().split("\t")[1])==0.0))
						powers_small_agg[h_index][maf_index][con_index]++;
					if(!(Double.parseDouble(br.readLine().split("\t")[1])==0.0))
						powers_exact_local[h_index][maf_index][con_index]++;
					if(!(Double.parseDouble(br.readLine().split("\t")[1])==0.0))
						powers_exact_agg[h_index][maf_index][con_index]++;
					if(!(Double.parseDouble(br.readLine().split("\t")[1])==0.0))
						powers_big_local[h_index][maf_index][con_index]++;
					if(!(Double.parseDouble(br.readLine().split("\t")[1])==0.0))
						powers_big_agg[h_index][maf_index][con_index]++;
				}				
			}
			for(int i=0;i<heritability.length;i++){
				for(int j=0;j<rare_maf.length;j++){
					for(int m=0;m<concentration_rare.length;m++){
						powers_small_local[i][j][m]/=counts[i][j][m];
						powers_small_agg[i][j][m]/=counts[i][j][m];
						powers_exact_local[i][j][m]/=counts[i][j][m];
						powers_exact_agg[i][j][m]/=counts[i][j][m];
						powers_big_local[i][j][m]/=counts[i][j][m];
						powers_big_agg[i][j][m]/=counts[i][j][m];
					}
				}		
			}
			String[] legend={"Local-kinship: Smaller Window","Aggregate: Smaller Window",
					"Local-kinship: Exact Window","Aggregate: Exact Window",
					"Local-kinship: Larger Window","Aggregate: Larger Window"};
			plot_folder=data_folder+"_plot/con";
			plot_folder_file=new File(plot_folder);
			if(!plot_folder_file.exists())plot_folder_file.mkdir();
			for(int i=0;i<heritability.length;i++){
				for(int j=0;j<rare_maf.length;j++){
					String title="heritability="+heritability[i]+"maf="+rare_maf[j];
					double[] x_values=concentration_rare.clone();
					double[][] y_values=new double[6][concentration_rare.length];
					String x_lab="Proportion of causal";
					String y_lab="Power";
					for(int m=0;m<concentration_rare.length;m++){
						y_values[0][m]=powers_small_local[i][j][m]*100;
						y_values[1][m]=powers_small_agg[i][j][m]*100;
						y_values[2][m]=powers_exact_local[i][j][m]*100;
						y_values[3][m]=powers_exact_agg[i][j][m]*100;
						y_values[4][m]=powers_big_local[i][j][m]*100;
						y_values[5][m]=powers_big_agg[i][j][m]*100;
					}
					MyLineChart plot=new MyLineChart(title, x_lab, y_lab, legend, x_values, y_values, 
							500, 500, plot_folder+"/"+title+".png");
				}
			}
			plot_folder=data_folder+"_plot/heritability";
			plot_folder_file=new File(plot_folder);
			if(!plot_folder_file.exists())plot_folder_file.mkdir();
			for(int m=0;m<concentration_rare.length;m++){			
				for(int j=0;j<rare_maf.length;j++){
					String title="Proportion="+concentration_rare[m]+"maf="+rare_maf[j];
					double[] x_values=heritability.clone();
					double[][] y_values=new double[6][heritability.length];
					String x_lab="Heritability";
					String y_lab="Power";
					for(int i=0;i<heritability.length;i++){
						y_values[0][i]=powers_small_local[i][j][m]*100;
						y_values[1][i]=powers_small_agg[i][j][m]*100;
						y_values[2][i]=powers_exact_local[i][j][m]*100;
						y_values[3][i]=powers_exact_agg[i][j][m]*100;
						y_values[4][i]=powers_big_local[i][j][m]*100;
						y_values[5][i]=powers_big_agg[i][j][m]*100;
					}
					MyLineChart plot=new MyLineChart(title, x_lab, y_lab, legend, x_values, y_values, 
							500, 500, plot_folder+"/"+title+".png");
				}
			}
			plot_folder=data_folder+"_plot/maf";
			plot_folder_file=new File(plot_folder);
			if(!plot_folder_file.exists())plot_folder_file.mkdir();
			for(int m=0;m<concentration_rare.length;m++){			
				for(int i=0;i<heritability.length;i++){
					String title="Proportion="+concentration_rare[m]+"heritability="+heritability[i];
					double[] x_values=rare_maf.clone();
					double[][] y_values=new double[6][rare_maf.length];
					String x_lab="Maximal MAF";
					String y_lab="Power";					
					for(int j=0;j<rare_maf.length;j++){
						y_values[0][j]=powers_small_local[i][j][m]*100;
						y_values[1][j]=powers_small_agg[i][j][m]*100;
						y_values[2][j]=powers_exact_local[i][j][m]*100;
						y_values[3][j]=powers_exact_agg[i][j][m]*100;
						y_values[4][j]=powers_big_local[i][j][m]*100;
						y_values[5][j]=powers_big_agg[i][j][m]*100;
					}
					MyLineChart plot=new MyLineChart(title, x_lab, y_lab, legend, x_values, y_values, 
							500, 500, plot_folder+"/"+title+".png");
				}
			}
		}catch(Exception e){e.printStackTrace();}
	}
	
	/*
	 * SIM.h0.1.min0.1.max0.2.num1.w10000.chr4.start14067857.end14091171.round0
	 */
	public static void make_plots_common(String data_folder){
		double[][] powers_local=new double[heritability.length][num_of_causal_common.length];
		double[][] powers_emmax=new double[heritability.length][num_of_causal_common.length];		
		double[][] counts=new double[heritability.length][num_of_causal_common.length];		
		try{
			String plot_folder=data_folder+"_plot";
			File plot_folder_file=new File(plot_folder);
			if(!plot_folder_file.exists())plot_folder_file.mkdir();
			File root=new File(data_folder);
			File[] subfolders=root.listFiles();
			for(int k=0;k<subfolders.length;k++){
				if(k%1000==0)System.out.println(k/1000+"K");
				String name=subfolders[k].toString().split("/")[subfolders[k].toString().split("/").length-1];
				//System.out.println(k+":"+name);
				double h=Double.parseDouble(name.substring(5, 8));				
				int max_index_name=name.indexOf("max");
				int num_index_name=name.indexOf("num");
				int w_index_name=name.indexOf("w25");
				//double maf=Double.parseDouble(name.substring(max_index_name+3,con_index_name-1));				
				double number=Double.parseDouble(name.substring(num_index_name+3,w_index_name-2));
//				System.out.println(h+"/"+number);
				int h_index=-1, num_index=-1;
				for(int p=0;p<heritability.length;p++)if(heritability[p]==h)h_index=p;
				for(int p=0;p<num_of_causal_common.length;p++)if(num_of_causal_common[p]==number)num_index=p;
				if(h_index==-1||num_index==-1){
					System.out.println("=====return====="+ number);
					return;
				}
				//System.out.println(h_index+"/"+maf_index+"/"+con_index);
				if(new File(subfolders[k]+"/summary.txt").exists()){
					BufferedReader br=new BufferedReader(new FileReader(subfolders[k]+"/summary.txt"));
					String line=br.readLine();
					if(line!=null){
						counts[h_index][num_index]++;
						if(!(Double.parseDouble(line.split("\t")[1])==0.0))
							powers_emmax[h_index][num_index]++;
						if(!(Double.parseDouble(br.readLine().split("\t")[1])==0.0))
							powers_local[h_index][num_index]++;
					}	
				}							
			}
			for(int i=0;i<heritability.length;i++){
				for(int j=0;j<num_of_causal_common.length;j++){					
						powers_local[i][j]/=counts[i][j];
						powers_emmax[i][j]/=counts[i][j];				
				}		
			}
			String[] legend={"Local-kinship","EMMAX"};
			plot_folder=data_folder+"_plot/num";
			plot_folder_file=new File(plot_folder);
			if(!plot_folder_file.exists())plot_folder_file.mkdir();
			for(int i=0;i<heritability.length;i++){				
				String title="Heritability="+heritability[i];
				double[] x_values=num_of_causal_common.clone();
				double[][] y_values=new double[2][num_of_causal_common.length];
				String x_lab="Number of causal";
				String y_lab="Power";
				for(int j=0;j<num_of_causal_common.length;j++){
					y_values[0][j]=powers_local[i][j]*100;
					y_values[1][j]=powers_emmax[i][j]*100;
				}
				MyLineChart plot=new MyLineChart(title, x_lab, y_lab, legend, x_values, y_values, 
						500, 500, plot_folder+"/"+title+".png");				
			}
			plot_folder=data_folder+"_plot/heritability";
			plot_folder_file=new File(plot_folder);
			if(!plot_folder_file.exists())plot_folder_file.mkdir();
			for(int j=0;j<num_of_causal_common.length;j++){					
				String title="Number of causal="+num_of_causal_common[j];
				double[] x_values=heritability.clone();
				double[][] y_values=new double[2][heritability.length];
				String x_lab="Heritability";
				String y_lab="Power";
				for(int i=0;i<heritability.length;i++){
					y_values[0][i]=powers_local[i][j]*100;
					y_values[1][i]=powers_emmax[i][j]*100;
				}
				MyLineChart plot=new MyLineChart(title, x_lab, y_lab, legend, x_values, y_values, 
						500, 500, plot_folder+"/"+title+".png");
				
			}			
		}catch(Exception e){e.printStackTrace();}
	}
	
	
	public static void main(String[] args) {
		int chr=2; // stars from 0
		int start_loc=4005001;
		int end_loc=4005000+5000;
		
		int[] total_num_of_causal={2, 3, 5, 10, 20, 50, 100, 150, 100};
		int[] num_allele_needed_to_have_efffet={1, 1, 1, 4, 4, 4, 4, 4, 1};
		double[] min_maf={0.3, 0.2, 0.15, 0.1, 0.04, 0, 0, 0, 0}, max_maf={0.35, 0.25, 0.2, 0.15, 0.1, 0.1, 0.1, 0.1, 0.01};
		 
		String regions_file="/Volumes/Projects/DATA/Annotations/arabidopsis/TAIR10_GFF3_genes.gene2Konly.tsv";
		String input_genotype="/Volumes/Projects/DATA/arabidopsis_genomes/swe180_ecker171_removebads.csv.num2.hdf5";
//		String global_kinship_file="/Volumes/Projects/DATA/arabidopsis_genomes/swe180_ecker171_removebads.csv.num2.hdf5.kinship.rescaled.ibs";
		
		String global_kinship_file="/Volumes/Projects/DATA/arabidopsis_genomes/kinship_files_swe180_ecker171_removebads/swe180_ecker171_removebads.2.rescaled.IBS";
//		String local_kinship_files="/Volumes/Projects/DATA/arabidopsis_genomes"+
//				"kinship_files_swe180_ecker171_removebads/100K/";
		String true_signal_file="";
		String temp_pheno_file="/Volumes/Projects/Local-kinship/simulations/main";
		String phenotype_file_folder_recessive="/Volumes/Projects/Local-kinship/simulations/sim_rare_recessive/";
		
		String phenotype_file_folder3="/Volumes/Projects/Local-kinship/simulations/sim_rare_additive_bidir";
		String phenotype_file_folder4="/Volumes/Projects/Local-kinship/simulations/sim_common_additive_bidir";
		String phenotype_file_folder5="/Volumes/Projects/Local-kinship/simulations/sim_rare_liability";
		String phenotype_file_folder6="/Volumes/Projects/Local-kinship/simulations/sim_common_liability";
		
		String phenotype_file_folder1="/Volumes/Projects/Local-kinship/simulations/sim_rare_dominant2";
		String phenotype_file_folder2="/Volumes/Projects/Local-kinship/simulations/sim_rare_recessive2";
		String phenotype_file_folder_common1="/Volumes/Projects/Local-kinship/simulations/additive_apr23";
		String phenotype_file_folder_common2="/Volumes/Projects/Local-kinship/simulations/epistasis_apr23";
		String phenotype_file_folder_common3="/Volumes/Projects/Local-kinship/simulations/dorminant_apr23";
		String phenotype_file_folder_common4="/Volumes/Projects/Local-kinship/simulations/liability_apr23";
		
		int num_replicates=1;
//		double heritability=0.9;
		int analysis_tiling_win_size=50000;
		int sim_win=10000;
		String tmep_output_folder="/Volumes/Projects/Local-kinship/simulations/have2look/";
		String typeIerror_folder="/Volumes/Projects/Local-kinship/simulations/typeIerror/";
		
		try{
//			int round=total_num_of_causal.length;
//			if(round!=num_allele_needed_to_have_efffet.length||round!=min_maf.length||round!=max_maf.length){
//				System.out.println("Wrong");
//				System.exit(0);
//			}
			//for(int r=round-1;r>=0;r--){
//			for(int r=0;r<round;r++){
//				sim_three(input_genotype, true_signal_file, temp_pheno_file+"."+r+".tsv", num_replicates, heritability, 
//						min_maf[r], max_maf[r], total_num_of_causal[r],	num_allele_needed_to_have_efffet[r], 
//						win_size, tmep_output_folder+"r"+r+"/", 
//						global_kinship_file, regions_file, chr, start_loc, end_loc);
//				
//			}
		}catch(Exception e){e.printStackTrace();}
		
//		simulate_phenotypes(phenotype_file_folder, input_genotype, regions_file, 10000);
		
		
//		calculate_typeIerror(regions_file, win_size, input_genotype, 0.05, 
//				typeIerror_folder+"random_phenotype.tsv", 1, typeIerror_folder, global_kinship_file);
		
		//int round_index=0;
		//phenotype_file_folder=args[0];
		
//		for(int round_index=0;round_index<500;round_index++){
//			RARE_power_coparion(phenotype_file_folder5, input_genotype, regions_file, global_kinship_file, sim_win, analysis_tiling_win_size, 
//				round_index);
//		}
//		
//		System.out.println("PLOT RARE:");
//		make_plots_rare(phenotype_file_folder1);
//		make_plots_rare(phenotype_file_folder2);
				
//		if(args.length<7){
//			System.out.println("<output_folder> <input_genotype> <regions_file> <global_kinship_file> <sim_win> " +
//					"<analysis_tiling_win_size>, <analysis_type>");
//			System.exit(0);
//		}
//		for(int round_index=0;round_index<500;round_index++){
//			COMMON_power_coparion(args[0], args[1], args[2], args[3], Integer.parseInt(args[4]), 
//					Integer.parseInt(args[5]), round_index, args[6]);
//		}
		
		
//		for(int round_index=0;round_index<500;round_index++){
//			COMMON_power_coparion(phenotype_file_folder_common4, input_genotype, regions_file, global_kinship_file, sim_win*50, 
//					analysis_tiling_win_size, round_index, analysis_type);
//		}

		System.out.println("PLOT COMMON:");
		make_plots_common(phenotype_file_folder_common1);
		make_plots_common(phenotype_file_folder_common2);
		make_plots_common(phenotype_file_folder_common3);
		make_plots_common(phenotype_file_folder_common4);
//		if(args.length!=1){
//			System.out.println("<folder>");
//			System.exit(0);
//		}make_plots_common(args[0]);
	
// under construction yet		
//		for(int round_index=0;round_index<500;round_index++){
//			COMMON_power_coparion_bk(phenotype_file_folder5, input_genotype, regions_file, global_kinship_file, sim_win*50, analysis_tiling_win_size, 
//					round_index);
//		}
		
//	FOLLOWING ARE DEBUG PROGRAMS 		
//		String debug="/Volumes/Projects/Local-kinship/simulations/sim_common_apr16_liability/SIM.h0.9.min0.05.max0.1.num3.0.w500000.chr1.start17066511.end18069297.round187/";
//		String results_output_folder=debug+"debug3/";
//		String simulated_phenotype_file=debug+"simulated_phenotype.tsv";
//		String summary =debug+"summary2.txt";
//		int chr0=0, start=17066511, end=18069297;
//		running_single_round_common_QUICK(regions_file,  analysis_tiling_win_size, input_genotype, 
//				 simulated_phenotype_file,  results_output_folder,  global_kinship_file,
//				 summary,	chr0, start, end);
	}

}
