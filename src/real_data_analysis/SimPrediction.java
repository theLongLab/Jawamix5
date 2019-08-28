package real_data_analysis;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.Random;

import mixedmodel.GBLUP;
import mixedmodel.KinshipMatrix;
import mixedmodel.MultiPhenotype;
import mixedmodel.Phenotype;
import mixedmodel.Regions;
import mixedmodel.SearchAndCalculateRegionalKinships;
import mixedmodel.VariantsDouble;
import myMathLib.Test;
import simulations.Simulator;

public class SimPrediction {

	public static double[] heritability={0.9,0.7,0.5,0.3,0.1};
	public static double[] concentration_rare={0.1,0.2,0.3,0.4,0.5};
	public static double[] rare_maf={0.01,0.02,0.03,0.04,0.05};	//for old plots
//	public static double[] rare_maf={0.01};	
	public static double[] common_maf={0.05,0.10};
//	public static double[] num_of_causal_common={1,2,3,5,8,11,15,20,25,30};
	public static double[] num_of_causal_common={2,3,4,5};
	
	public static void test_corr(String genotye_file, String phenotype_file, String full_kinship_file, 
			String local_K_folder, int obs_sample_size, int first_sample_size,
			double h2, int chr, int start_loc, int end_loc){
		try{
			String output_file=phenotype_file+".corr.txt";
			BufferedWriter bw=new BufferedWriter(new FileWriter(output_file));
			bw.write("globle\tlocal_sum\tlocal_beta\n");
						
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
			double g_record=0, l_record=0, ll_record=0;
			for(int round=0;round<100;round++){
				System.out.println("==== Round: "+(round+1)+"======");
				String[] obs_ids=GBLUP.select_ids_random(obs_sample_size, full_pheno, geno_full);
				SearchAndCalculateRegionalKinships assigned_regions=new SearchAndCalculateRegionalKinships(regions, var_comp, 
						local_K_folder, new VariantsDouble(genotye_file));
				System.out.println("Assigning region done.");
				double[] corrs=GBLUP.test_batch(full_pheno, geno_full, kinship_full, obs_ids, 
						h2, assigned_regions.regions, assigned_regions.var_comp, 
					assigned_regions.reginal_full_kinship_files, null, first_sample_size);
				g_record+=corrs[0];
				l_record+=corrs[1];
				ll_record+=corrs[2];
				System.out.println("===="+g_record+"/"+l_record+"/"+ll_record+"/"+(round+1)+"======");
				bw.write(corrs[0]+"\t"+corrs[1]+"\t"+corrs[2]+"\n");
			}bw.close();
			//GBLUP.test_batch(full_pheno, kinship_full, obs_ids, h2, null);
		}catch(Exception e){e.printStackTrace();}
		
	}
	
	public static void COMMON_power_coparion(String working_folder, String input_genotype, String global_kinship_file,
			int sim_win, int round_index, String analysis_type, int obs_sample_size, int first_sample_size){		
		
		//double[][][] correlations=new double[heritability.length][num_of_causal_common.length][];
		VariantsDouble genotype=new VariantsDouble(input_genotype);
		
		//Random rg = new Random();
		//
		try{			
			for(int h_index=0;h_index<heritability.length;h_index++){				
				double max_maf=common_maf[1];
				double min_maf=common_maf[0];				
				for(int num_index=0;num_index<num_of_causal_common.length;num_index++){
					int chr=(int)(Test.randomNumber()*genotype.num_chrs);
					if(chr==genotype.num_chrs)chr--;
					int max_loc=genotype.locations[chr][genotype.locations[chr].length-1];
					int start_loc=(int)(max_loc*Test.randomNumber());
					if(start_loc>max_loc-2*sim_win)start_loc=start_loc-sim_win;
					int end_loc=start_loc+sim_win;

					String simulated_phenotype_file_folder=working_folder+"/SIM.h"+heritability[h_index]+".min"+min_maf+".max"+max_maf+
							".num"+num_of_causal_common[num_index]+".w"+sim_win+
							".chr"+(chr+1)+".start"+start_loc+".end"+end_loc+".round"+round_index+"/";
					File the_folder=new File(simulated_phenotype_file_folder);
					if(!the_folder.exists())the_folder.mkdir();
					Simulator sim=new Simulator(input_genotype, heritability[h_index], min_maf, max_maf, 
							chr, start_loc, end_loc);
					String phenotype_file_simulated=simulated_phenotype_file_folder+"simulated_phenotype.tsv";
					String causal_file_simulated=simulated_phenotype_file_folder+"causal.tsv";
					if(analysis_type.equals("dorminant"))
						sim.heterogeneity(Simulator.heterogeneityDorminant, (int)num_of_causal_common[num_index], 
								1, 1, false, phenotype_file_simulated, causal_file_simulated);
					else if(analysis_type.equals("epistasis")) 
						sim.heterogeneity(Simulator.heterogeneityEpistasis, (int)num_of_causal_common[num_index], 
								1, 4, false, phenotype_file_simulated, causal_file_simulated);
					else if(analysis_type.equals("additive")) 
						sim.additive(Simulator.additiveBidirection, (int)num_of_causal_common[num_index], 
								-1, 1, false, phenotype_file_simulated, causal_file_simulated);
					else if(analysis_type.equals("liability")) 
						sim.liabilityThreshold(Simulator.liabilityThreshold, (int)num_of_causal_common[num_index], 
								1, false, phenotype_file_simulated, causal_file_simulated);
					System.out.println("Simulating pheno done.");
					test_corr(input_genotype, phenotype_file_simulated, global_kinship_file, simulated_phenotype_file_folder, 
							obs_sample_size, first_sample_size, heritability[h_index], chr, start_loc, end_loc);
				}
			}			
		}catch(Exception e){e.printStackTrace();}
	}
	
	public static void summarize_plot(String working_folder){
		try{
			
		}catch(Exception e){e.printStackTrace();}
	}
	public static void main(String[] args) {
		String genotype_file2="/Volumes/Projects/DATA/1000G/1000gdata/g1k_all.hdf5";
		String global_kinship_file2="/Volumes/Projects/DATA/1000G/1000gdata/g1k_all.K.RRM";
		String working_folder_top2="/Users/quanlong/Documents/projects2/predictions/simulation/1000g/";
		int obs_sample_size2=500; 
		int first_sample_size2 =300;
		
//		String regions_file="/Volumes/Projects/DATA/Annotations/arabidopsis/TAIR10_GFF3_genes.gene2Konly.tsv";
		String genotype_file="/Volumes/Projects/DATA/arabidopsis_genomes/swe180_ecker171_removebads.csv.num2.hdf5";
		String working_folder_top="/Users/quanlong/Documents/projects2/predictions/simulation/at/";
		
		String global_kinship_file="/Volumes/Projects/DATA/arabidopsis_genomes/kinship_files_swe180_ecker171_removebads/swe180_ecker171_removebads.2.rescaled.IBS";
//	
		int sim_win=100000;
		String[] types={"additive", "dorminant", "epistasis", "liability"};
		for(int round_index=1;round_index<=10;round_index++){
			for(int type_index=0;type_index<types.length;type_index++){
				String analysis_type=types[type_index];
				String working_folder2=working_folder_top2+analysis_type+"/";			
				COMMON_power_coparion(working_folder2, genotype_file2, global_kinship_file2, sim_win, round_index, 
						analysis_type, obs_sample_size2, first_sample_size2);
			}
		}
		

	}

}
