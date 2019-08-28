package real_data_analysis;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import mixedmodel.KinshipMatrix;
import mixedmodel.MultiPhenotype;
import mixedmodel.VarComp_Result;
import mixedmodel.VarianceComponent;
import mixedmodel.VariantsDouble;
import myMathLib.StatFuncs;
import myMathLib.Test;
import simulations.Simulator;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;

public class Heritability {
	
	public static double[] heritability={0.9,0.7,0.5,0.3,0.1};
	public static double[] concentration_rare={0.1,0.2,0.3,0.4,0.5};
	public static double[] rare_maf={0.01,0.02,0.03,0.04,0.05};	//for old plots
//	public static double[] rare_maf={0.01};	
	public static double[] common_maf={0.3,0.5};
	public static double[] num_of_causal_common={1,2,3,5,8,11,15,20,25,30,50,100};
//	public static double[] num_of_causal_common={2,3,4,5};
	

	public static void COMMON_interval(String working_folder, String input_genotype, String global_kinship_file,
			int round_index, String analysis_type, int num_rep, String summary_file){		
		
		//double[][][] correlations=new double[heritability.length][num_of_causal_common.length][];
		VariantsDouble genotype=new VariantsDouble(input_genotype);
		
		//Random rg = new Random();
		//
		try{			
			BufferedWriter bw=new BufferedWriter(new FileWriter(summary_file));
			bw.write("round\th2\tnum\tPearsonsCorrelation\tSpearmansCorrelation\n");
			for(int h_index=0;h_index<heritability.length;h_index++){				
				double max_maf=common_maf[1];
				double min_maf=common_maf[0];				
				for(int num_index=num_of_causal_common.length-1;num_index>=0;num_index--){
					String simulated_phenotype_file_folder=working_folder+"/SIM.h"+heritability[h_index]+".min"+min_maf+".max"+max_maf+
							".num"+num_of_causal_common[num_index]+".round"+round_index+"/";
					File the_folder=new File(simulated_phenotype_file_folder);
					if(!the_folder.exists())the_folder.mkdir();
					Simulator sim=new Simulator(input_genotype, heritability[h_index], min_maf, max_maf);
					String phenotype_file_simulated=simulated_phenotype_file_folder+"simulated_phenotype.tsv";
					if(analysis_type.equals("dorminant"))
						sim.heterogeneity(Simulator.heterogeneityDorminant, (int)num_of_causal_common[num_index], 
								num_rep, 1, false, phenotype_file_simulated);
					else if(analysis_type.equals("epistasis")) 
						sim.heterogeneity(Simulator.heterogeneityEpistasis, (int)num_of_causal_common[num_index], 
								num_rep, 4, false, phenotype_file_simulated);
					else if(analysis_type.equals("additive")) 
						sim.additive(Simulator.additiveBidirection, (int)num_of_causal_common[num_index], -1, num_rep, false, phenotype_file_simulated);
					else if(analysis_type.equals("liability")) 
						sim.liabilityThreshold(Simulator.liabilityThreshold, (int)num_of_causal_common[num_index], num_rep, false, phenotype_file_simulated);
					System.out.println("Simulating pheno done.");
					VarianceComponent.heritability_fullrank(new KinshipMatrix(global_kinship_file), 
							new MultiPhenotype(phenotype_file_simulated), phenotype_file_simulated+".h2.csv");
					double[] corrs=calcualte_correlation(phenotype_file_simulated+".h2.csv", num_rep,heritability[h_index]);
					bw.write(round_index+"\t"+heritability[h_index]+"\t"+num_of_causal_common[num_index]+"\t"+corrs[0]+"\t"+corrs[1]+"\n");
					bw.flush();
				}
			}bw.close();			
		}catch(Exception e){e.printStackTrace();}
	}
	
	public static double[] calcualte_correlation(String file, int num, double real_h2){
		double[] corrs=new double[2];
		try{
			PearsonsCorrelation c0=new PearsonsCorrelation();
			SpearmansCorrelation c1=new SpearmansCorrelation();
			double[] h2=new double[num];
			double[] fi=new double[num];
			BufferedReader br=new BufferedReader(new FileReader(file));
			String line=br.readLine();line=br.readLine();
			int index=0;
			while(line!=null){
				String[] tmp=line.split(",");
				h2[index]=Math.abs(Double.parseDouble(tmp[1])-real_h2);
				fi[index]=Math.abs(Double.parseDouble(tmp[3])-Double.parseDouble(tmp[4]));
				index++;
				line=br.readLine();
			}
			//System.out.println("===="+c0.correlation(h2, fi)+"===="+c1.correlation(h2, fi)+"====");
			corrs[0]=c0.correlation(h2, fi);
			corrs[1]=c1.correlation(h2, fi);
		}catch(Exception e){e.printStackTrace();}
		return corrs;
	}
	
	public static void main(String[] args) {
		
		String folder="/Users/quanlong/Documents/projects2/interval/";
		
		String genotype_file="/Volumes/Projects/DATA/arabidopsis_genomes/swe180_ecker171_removebads.csv.num2.hdf5";
		String global_kinship_file="/Volumes/Projects/DATA/arabidopsis_genomes/kinship_files_swe180_ecker171_removebads/swe180_ecker171_removebads.2.rescaled.IBS";
		String global_kinship_file_common="/Volumes/Projects/DATA/arabidopsis_genomes/kinship_files_swe180_ecker171_removebads/swe180_ecker171_removebads.K.05.RRM";
		String working_folder_top=folder+"at2/";
		
		
		String genotype_file2="/Volumes/Projects/DATA/1000G/1000gdata/g1k_all.hdf5";
		String global_kinship_file2="/Volumes/Projects/DATA/1000G/1000gdata/g1k_all.K.RRM";
		String global_kinship_file2_common="/Volumes/Projects/DATA/1000G/1000gdata/g1k_all.K.05.RRM";
		String working_folder_top2=folder+"1kg2/";
		
		String[] types={"additive", "dorminant", "epistasis", "liability"};
		int num_rep=1000;
		for(int round_index=1;round_index<=10;round_index++){
			for(int type_index=0;type_index<types.length;type_index++){
				String analysis_type=types[type_index];
				String working_folder=working_folder_top2+analysis_type+"/";	
				String summary_file=working_folder+analysis_type+".summary.r"+round_index+".txt";
				COMMON_interval(working_folder, genotype_file2, global_kinship_file2_common, round_index, 
						analysis_type, num_rep,summary_file);
			}
		}
		
		
		
		
		

	}

}
