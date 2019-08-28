package real_data_analysis;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Random;

import mixedmodel.EMMAX;
import mixedmodel.FaSTLMM_LocalKinship;
import mixedmodel.GBLUP;
import mixedmodel.KinshipMatrix;
import mixedmodel.MultiPhenotype;
import mixedmodel.Phenotype;
import mixedmodel.RareAnalyzerAggregate;
import mixedmodel.Regions;
import mixedmodel.SearchAndCalculateRegionalKinships;
import mixedmodel.VariantsDouble;
import myMathLib.StatFuncs;
import myMathLib.Test;
import simulations.Simulator;

public class PredictionGWAS {
	
	// genetic models
	public static double[] heritability={0.7,0.5,0.3,0.25,0.2,0.15,0.1,0.05};
//	public static double[] concentration_rare={0.03, 0.05, 0.08, 0.12, 0.2, 0.4};
//	public static double[] rare_maf={0.01,0.02,0.03,0.04,0.05};	//for old plots
//	public static double[] rare_maf={0.002,0.005};	
//	public static double[] common_maf={0.05,0.2}; // input as args[] in this class
//	public static double[] num_of_causal_common={1,2,3,5,8,11,15,20,25,30};
	public static double[] num_of_causal_common={2,3,4,5,6,7,8,9,10};

	// GWAS
	//public static double emmax_arabidopsis_threshold=0.019679907/1416064;
	//public static double aggregate_arabidopsis_threshold=0.01212773/12929;
	//public static double local_gene_arabidopsis_threshold=0.052253921/12929;
	//public static double local_win50k_arabidopsis_threshold=0.058879085/4800;
	
	public static double emmax_1kg_threshold=0.0644534742932854/39375579;
	public static double local_win200k_1kg_threshold=0.142046626/28564;
	public static double local_win100k_1kg_threshold=0.142713988/57006;
	public static double local_win50k_1kg_threshold= 0.143789665/113748;
	public static double aggre_win200k_1kg_threshold=0.052616435/28564;
	public static double aggre_win100k_1kg_threshold=0.054369476/57006;
	public static double aggre_win50k_1kg_threshold= 0.0558693/113748;

	public static void running_single_round_common_QUICK(int analysis_win_size, String genotype_file, 
			String simulated_phenotype_file, String results_output_folder, String global_kinship_file,
			String output_summary_file, int chr, int start, int end){
		try{
			//VariantsDouble genotype=new VariantsDouble(genotype_file);
			Phenotype phenotype=new MultiPhenotype(simulated_phenotype_file).phenotypes[0];
			FaSTLMM_LocalKinship calculator=new FaSTLMM_LocalKinship(phenotype, genotype_file, 
					global_kinship_file);
			int num_chr=23;
			//int[][][] regions=(new Regions(regions_file, calculator.genotype)).region_coordinates;
			int[][][] region_only1=new int[num_chr][][];
			String[][] region_name=new String[num_chr][];
			for(int chr_index=0;chr_index<num_chr;chr_index++){
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
			bw.write("EMMAX."+start+"."+end+"\t"+EntranceSim.check_power_specific_region_only(emmax_outputfile, emmax_1kg_threshold)
					+"\t"+(-Math.log10(emmax_1kg_threshold))+"\n");
			bw.write("Local_gene."+start+"."+end+"\t"+EntranceSim.check_power(local_gene_outfile, chr, start, end, local_win100k_1kg_threshold)
					+"\t"+(-Math.log10(local_win100k_1kg_threshold))+"\n");
			bw.close();			
		}catch(Exception e){e.printStackTrace();}
	}
	
	/*
	 * MAIN Rare function for power
	 */
	public static void rare_power_comp(String working_folder, String input_genotype, String global_kinship_file,
			int sim_win, int round_index, String analysis_type, int obs_sample_size, int first_sample_size, 
			double[] rare_mafs, double[] concentration_rare){		
		VariantsDouble genotype=new VariantsDouble(input_genotype);
		Random rg = new Random();		//
		try{			
			for(int h_index=0;h_index<heritability.length;h_index++){
				double max_maf=rare_mafs[1];
				double min_maf=rare_mafs[0];
				for(int concentration_index=0;concentration_index<concentration_rare.length;concentration_index++){
						int chr=(int)(rg.nextDouble()*genotype.num_chrs);
						if(chr==genotype.num_chrs)chr--;
						int max_loc=genotype.locations[chr][genotype.locations[chr].length-1];
						int start_loc=(int)(max_loc*rg.nextDouble());
						if(start_loc>max_loc-2*sim_win)start_loc=start_loc-sim_win;
						int end_loc=start_loc+sim_win;						
						String simulated_phenotype_file_folder=working_folder+"/SIM.h"+heritability[h_index]+".min"+min_maf+".max"+max_maf+
								".con"+concentration_rare[concentration_index]+".w"+sim_win+
								".chr"+(chr+1)+".start"+start_loc+".end"+end_loc+".round"+round_index+"/";
						File the_folder=new File(simulated_phenotype_file_folder);
						if(!the_folder.exists())the_folder.mkdirs();
						Simulator sim=new Simulator(input_genotype, heritability[h_index], min_maf, max_maf, chr, start_loc, end_loc);
						String phenotype_file_simulated=simulated_phenotype_file_folder+"simulated_phenotype.tsv";
						String causal_file_prefix=simulated_phenotype_file_folder+"/causal";
						if(analysis_type.equals("dorminant")){
							sim.heterogeneity2(Simulator.heterogeneityDorminant, concentration_rare[concentration_index], 
									1, 1, false, phenotype_file_simulated, causal_file_prefix);
						}else if(analysis_type.equals("epistasis")){ 
							sim.heterogeneity2(Simulator.heterogeneityEpistasis, concentration_rare[concentration_index], 
									1, 4, false, phenotype_file_simulated, causal_file_prefix);
						}else if(analysis_type.equals("additive")){ 
							sim.additive2(Simulator.additiveBidirection, concentration_rare[concentration_index], 
									-1, 1, false, phenotype_file_simulated, causal_file_prefix);
						}else if(analysis_type.equals("liability")){ 
							sim.liabilityThreshold2(Simulator.liabilityThreshold, concentration_rare[concentration_index], 
									1, false, phenotype_file_simulated, causal_file_prefix);
						}System.out.println("Simulating pheno done.");
						String causal_file_simulated=causal_file_prefix+".0.txt";										
						if(EntranceSim.phenotype_check_passed(phenotype_file_simulated)){
							BufferedWriter bw=new BufferedWriter(new FileWriter(simulated_phenotype_file_folder+"/summary.txt"));		
							running_single_round_rare_QUICK(input_genotype,  
									phenotype_file_simulated, simulated_phenotype_file_folder, global_kinship_file, 
									bw, chr, start_loc-sim_win/2,end_loc+sim_win/2, 
									emmax_1kg_threshold, local_win200k_1kg_threshold, aggre_win200k_1kg_threshold);							
							running_single_round_rare_QUICK(input_genotype,  
									phenotype_file_simulated, simulated_phenotype_file_folder, global_kinship_file, 
									bw, chr, start_loc,end_loc, 
									emmax_1kg_threshold, local_win100k_1kg_threshold, aggre_win100k_1kg_threshold);							
							running_single_round_rare_QUICK(input_genotype,  
									phenotype_file_simulated, simulated_phenotype_file_folder, global_kinship_file, 
									bw, chr, start_loc+sim_win/4,end_loc-sim_win/4, 
									emmax_1kg_threshold, local_win50k_1kg_threshold, aggre_win50k_1kg_threshold);
							bw.close();	
							test_corr(input_genotype, phenotype_file_simulated, global_kinship_file, simulated_phenotype_file_folder, 
									simulated_phenotype_file_folder+"/summary.txt", causal_file_simulated,
									obs_sample_size, first_sample_size, heritability[h_index], chr, start_loc, end_loc);	
						}
										
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
	 * Rare, EMMAX, Local-window, local-region for specific regions ONLY. QUICK!
	 */
	public static void running_single_round_rare_QUICK(String genotype_file, String simulated_phenotype_file, String results_output_folder, 
			String global_kinship_file,	BufferedWriter bw,	int chr, int start, int end, 
			double emmx_threshold, double local_threshold, double aggregate_threshold){
		try{
			//VariantsDouble genotype=new VariantsDouble(genotype_file);
			Phenotype phenotype=new MultiPhenotype(simulated_phenotype_file).phenotypes[0];
			FaSTLMM_LocalKinship calculator=new FaSTLMM_LocalKinship(phenotype, genotype_file, 
					global_kinship_file);
			int num_chr=23;
			//int[][][] regions=(new Regions(regions_file, calculator.genotype)).region_coordinates;
			int[][][] region_only1=new int[num_chr][][];
			String[][] region_name=new String[num_chr][];
			for(int chr_index=0;chr_index<num_chr;chr_index++){
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
			String emmax_outputfile =results_output_folder+"/EMMAX.0_"+phenotype.phe_id+".top";
			EMMAX.emmax_analysis_regions(genotype_file, simulated_phenotype_file, global_kinship_file, results_output_folder,  
					1000000, 10, 0, 0.05, false, region_only1);
			//local
			String local_gene_outfile=results_output_folder+"/Local_VO."+phenotype.phe_id+".gene"+".csv";								
			calculator.analysis_specified_region(local_gene_outfile, region_only1, false);				
			bw.write("EMMAX."+start+"."+end+"\t"+EntranceSim.check_power_specific_region_only(emmax_outputfile, 
					emmx_threshold)+"\t"+(-Math.log10(emmx_threshold))+"\n");			
			bw.write("Local_gene."+start+"."+end+"\t"+EntranceSim.check_power(local_gene_outfile, chr, start, end, local_threshold)
					+"\t"+(-Math.log10(local_threshold))+"\n");
			bw.write("Aggregate."+start+"."+end+"\t"+EntranceSim.check_power(aggregate_outfile, chr, start, end, aggregate_threshold)
					+"\t"+(-Math.log10(aggregate_threshold))+"\n");
//			bw.close(); DO NOT CLOSE BW!!
			
		}catch(Exception e){e.printStackTrace();}
	}
	
	
	public static void common_power_comp(String working_folder, String input_genotype, String global_kinship_file,
			int sim_win, int round_index, String analysis_type, int obs_sample_size, int first_sample_size, 
			double[] common_maf){		
		VariantsDouble genotype=new VariantsDouble(input_genotype);		
		Random rg = new Random();
		try{			
			for(int h_index=0;h_index<heritability.length;h_index++){				
				double max_maf=common_maf[1];
				double min_maf=common_maf[0];				
				for(int num_index=0;num_index<num_of_causal_common.length;num_index++){
					int chr=(int)(rg.nextDouble()*genotype.num_chrs);
					if(chr==genotype.num_chrs)chr--;
					int max_loc=genotype.locations[chr][genotype.locations[chr].length-1];
					int start_loc=(int)(max_loc*rg.nextDouble());
					if(start_loc>max_loc-2*sim_win)start_loc=start_loc-sim_win;
					int end_loc=start_loc+sim_win;
					String simulated_phenotype_file_folder=working_folder+"/SIM.h"+heritability[h_index]+".min"+min_maf+".max"+max_maf+
							".num"+num_of_causal_common[num_index]+".w"+sim_win+
							".chr"+(chr+1)+".start"+start_loc+".end"+end_loc+".round"+round_index+"/";
					File the_folder=new File(simulated_phenotype_file_folder);
					if(!the_folder.exists())the_folder.mkdirs();
					Simulator sim=new Simulator(input_genotype, heritability[h_index], min_maf, max_maf, 
							chr, start_loc, end_loc);
					String phenotype_file_simulated=simulated_phenotype_file_folder+"simulated_phenotype.tsv";
					String causal_file_prefix=simulated_phenotype_file_folder+"/causal";
					if(analysis_type.equals("dorminant")){
						sim.heterogeneity(Simulator.heterogeneityDorminant, (int)num_of_causal_common[num_index], 
								1, 1, false, phenotype_file_simulated, causal_file_prefix);
					}else if(analysis_type.equals("epistasis")){ 
						sim.heterogeneity(Simulator.heterogeneityEpistasis, (int)num_of_causal_common[num_index], 
								1, 4, false, phenotype_file_simulated, causal_file_prefix);
					}else if(analysis_type.equals("additive")){ 
						sim.additive(Simulator.additiveBidirection, (int)num_of_causal_common[num_index], 
								-1, 1, false, phenotype_file_simulated, causal_file_prefix);
					}else if(analysis_type.equals("liability")){ 
						sim.liabilityThreshold(Simulator.liabilityThreshold, (int)num_of_causal_common[num_index], 
								1, false, phenotype_file_simulated, causal_file_prefix);
					}else if(analysis_type.equals("canalization")){
						sim.local_interactions(Simulator.canalization, 1, phenotype_file_simulated, causal_file_prefix);
					}else if(analysis_type.equals("signepistasis")){
						sim.local_interactions(Simulator.signepistasis, 1, phenotype_file_simulated, causal_file_prefix);
					}else if(analysis_type.equals("compensation")){
						sim.local_interactions(Simulator.compensation, 1, phenotype_file_simulated, causal_file_prefix);
					}
					System.out.println("Simulating pheno done.");
					int analysis_tiling_win_size=sim_win;
					String causal_file_simulated=causal_file_prefix+".0.txt";
					if(EntranceSim.phenotype_check_passed(phenotype_file_simulated)){
						running_single_round_common_QUICK(analysis_tiling_win_size, input_genotype,  
										phenotype_file_simulated, simulated_phenotype_file_folder, global_kinship_file, 
										simulated_phenotype_file_folder+"/summary.txt", chr, start_loc,end_loc);
						test_corr(input_genotype, phenotype_file_simulated, global_kinship_file, simulated_phenotype_file_folder, 
								simulated_phenotype_file_folder+"/summary.txt", causal_file_simulated,
								obs_sample_size, first_sample_size, heritability[h_index], chr, start_loc, end_loc);						
					}else{System.out.println("phenotype check NOT passed");}
				}
			}			
		}catch(Exception e){e.printStackTrace();}
	}
	
	/*
	 * causal_file is formated as: "Chr(FromOne)\tIndexInData(FromZero)\tLocation\n"
	 */
	public static double calcualte_distance_LD(String causal_file, VariantsDouble geno_full, int chr, int start_loc, int end_loc){
		double ld=Double.NaN;
		try{
			BufferedReader br=new BufferedReader(new FileReader(causal_file));
			String line=br.readLine();line=br.readLine();// skip header
			ArrayList<double[]> casual=new ArrayList<double[]>();
			while(line!=null){
				String[] tmp=line.split("\t");
				casual.add(geno_full.load_one_variant_by_index(Integer.parseInt(tmp[0])-1, Integer.parseInt(tmp[1])));
				line=br.readLine();
			}
			double[][] X=new double[casual.size()][];
			for(int i=0;i<X.length;i++)X[i]=casual.get(i);
			double[][] Y=geno_full.load_variants_in_region(chr, start_loc, end_loc);
			ld=StatFuncs.distance_correlation(X, Y);
		}catch(Exception e){e.printStackTrace();}
		return ld;
	}
	
	/*
	 * the summary.txt file generated by the comparison between local-k and emmax. 
	 */
	public static ArrayList<int[]> extract_the_marker(String summary_file){
		ArrayList<int[]> output=new ArrayList<int[]>();
		try{
			BufferedReader br=new BufferedReader(new FileReader(summary_file));
			String line[]=br.readLine().split("\t")[1].split("/");
			int[] the_marker=new int[2];
			the_marker[0]=Integer.parseInt(line[1]);
			the_marker[1]=Integer.parseInt(line[2]);
			output.add(the_marker);
		}catch(Exception e){e.printStackTrace();}
		return output;
	}
	
	public static void test_corr(String genotye_file, String phenotype_file, String full_kinship_file, 
			String local_K_folder, String summary_file, String causal_file,
			int obs_sample_size, int first_sample_size,	double h2, int chr, int start_loc, int end_loc){
		try{
			String output_file=phenotype_file+".corr.txt";
			BufferedWriter bw=new BufferedWriter(new FileWriter(output_file));
			bw.write("reg\tgblup\trgblup_sum\trgblup_beta\n");
						
			//String genotye_file=data_folder+"NFBC.num.ch1-22.csv.hdf5";
			//String phenotype_file=data_folder+"NFBC.phenot.tsv";
			//String full_kinship_file=genotye_file+".RRM";
			//String analysis_folder="/Users/quanlong/Documents/projects2/predictions/LDL/";
			
			Phenotype full_pheno=new MultiPhenotype(phenotype_file).phenotypes[0];
			KinshipMatrix kinship_full=new KinshipMatrix(full_kinship_file);
			VariantsDouble geno_full=new VariantsDouble(genotye_file);
			//String local_K_folder="/Users/quanlong/Documents/projects2/predictions/LDL/local_Ks/";
			
			double[] var_comp=new double[1];
			var_comp[0]=h2;
			int[][] regions=new int[1][3];
			regions[0][0]=chr; regions[0][1]=start_loc; regions[0][2]=end_loc;
	//		int obs_sample_size=280;
			double g_record=0, l_record=0, ll_record=0, reg_record=0;
			double dist_LD=calcualte_distance_LD(causal_file, geno_full, chr, start_loc, end_loc);
			for(int round=0;round<10;round++){
				System.out.println("==== Round: "+(round+1)+"======");
				String[] obs_ids=GBLUP.select_ids_random(obs_sample_size, full_pheno, geno_full);
				SearchAndCalculateRegionalKinships assigned_regions=new SearchAndCalculateRegionalKinships(regions, var_comp, 
						local_K_folder, new VariantsDouble(genotye_file));
				System.out.println("Assigning region done.");
				double[] corrs=GBLUP.test_batch(full_pheno, geno_full, kinship_full, obs_ids, 
						h2, assigned_regions.regions, assigned_regions.var_comp, 
					assigned_regions.reginal_full_kinship_files, null, first_sample_size);
				double reg_corr=GBLUP.multiple_regression_prediction(extract_the_marker(summary_file), full_pheno, geno_full, obs_ids)[0];
				g_record+=corrs[0];
				l_record+=corrs[1];
				ll_record+=corrs[2];
				reg_record+=reg_corr;
				System.out.println("===="+reg_record+"/"+g_record+"/"+l_record+"/"+ll_record+"/"+(round+1)+"======");
				bw.write(reg_corr+"\t"+corrs[0]+"\t"+corrs[1]+"\t"+corrs[2]+"\t"+dist_LD+"\n");
			}bw.close();
			//GBLUP.test_batch(full_pheno, kinship_full, obs_ids, h2, null);
		}catch(Exception e){e.printStackTrace();}
		
	}
	
	public static void main(String[] args) {
		if(args.length<9){
			System.out.println("<genotype> <g-kinship> <model> <working_folder> \n"
					+ "<obs_sample_size> <first_sample_size> <sim_win (100000)> "
					+ "<min_MAF> <max_MAF> [con-reare]");
			System.exit(0);
		}
		String genotype_file=args[0];
		String global_kinship_file=args[1];
		String model=args[2];
		String working_folder=args[3];
		int obs_sample_size=Integer.parseInt(args[4]);
		int first_sample_size=Integer.parseInt(args[5]);
		int sim_win=Integer.parseInt(args[6]); //100000
		double[] maf_cutoffs=new double[2];
		maf_cutoffs[0]=Double.parseDouble(args[7]);
		maf_cutoffs[1]=Double.parseDouble(args[8]);
		
		String genotype_file2="/Volumes/Projects/DATA/1000G/1000gdata/g1k_all.hdf5";
		String global_kinship_file2="/Volumes/Projects/DATA/1000G/1000gdata/g1k_all.K.RRM";
		String working_folder_top2="/Users/quanlong/Documents/projects2/predictions/simulation/1000g/";
//		int obs_sample_size2=500; 
//		int first_sample_size2 =300;
		
		for(int round_index=1;round_index<=100;round_index++){
			common_power_comp(working_folder, genotype_file, global_kinship_file, sim_win, round_index, 
						model, obs_sample_size, first_sample_size, maf_cutoffs);
			
		}

//		if(args.length==10){
//			double[] con_rare=new double[1];
//			con_rare[0]=Double.parseDouble(args[9]);
//			for(int round_index=1;round_index<=100;round_index++){
//				rare_power_comp(working_folder, genotype_file, global_kinship_file, sim_win, round_index, 
//							model, obs_sample_size, first_sample_size, maf_cutoffs, con_rare);				
//			}
//		}
			
		
	}

}
