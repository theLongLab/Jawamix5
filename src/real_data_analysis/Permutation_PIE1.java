package real_data_analysis;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.Random;

import org.apache.commons.math3.linear.FieldLUDecomposition;

import mixedmodel.FaSTLMM_LocalKinship;
import mixedmodel.MultiPhenotype;
import mixedmodel.Phenotype;
import mixedmodel.Regions;
import mixedmodel.VariantsDouble;
import simulations.Simulator;

public class Permutation_PIE1 {
	public static void permutation(String input_genotype, String phenotype_full, String global_kinship_file,
			int phe_index, int[][][] regions, String local_gene_outfile){		
		
		Phenotype phenotype=new MultiPhenotype(phenotype_full).phenotypes[phe_index];
		
		int round=100000;
		double[] all_p=new double[round];
		try{			
			FaSTLMM_LocalKinship calculator=new FaSTLMM_LocalKinship(phenotype, input_genotype, 
					global_kinship_file);
			calculator.analysis_specified_region(local_gene_outfile, regions, false);
			BufferedReader br=new BufferedReader(new FileReader(local_gene_outfile));	
			br.readLine();br.readLine();
			double pvalue=Double.parseDouble(br.readLine().split(", ")[2]);
			System.out.println("Real-P: "+pvalue);
			int smaller=0;
			for(int i=0;i<round;i++){
				if(i%100==0){
					System.out.println(i+"/"+smaller);
				}
				phenotype=new MultiPhenotype(phenotype_full).phenotypes[phe_index];
				phenotype.permute();
				calculator=new FaSTLMM_LocalKinship(phenotype, input_genotype, 
						global_kinship_file);
				calculator.analysis_specified_region(local_gene_outfile, regions, false);
				br=new BufferedReader(new FileReader(local_gene_outfile));	
				br.readLine();br.readLine();
				double permute_pvalue=Double.parseDouble(br.readLine().split(", ")[2]);
				if(permute_pvalue<pvalue){
					smaller++;
					System.out.println("Smaller: "+pvalue+"/"+i+"th/"+smaller);
				}else{
					System.out.println("Nothing: "+permute_pvalue+"/"+i+"th");
				}
				
			}
		}catch(Exception e){e.printStackTrace();}
	}
	
	public static void main(String[] args){
		String local_gene_outfile="/Volumes/Projects/Local-kinship/germ/tmp.csv";
		String phenotype_full="/Volumes/Projects/DATA/arabidopsis_phenotypes/big_table_gmi/gen_arch_traits_70_plus.useful.tsv";
		int phe_indxe=29;
		String input_genotype="/Volumes/Projects/DATA/arabidopsis_genomes/swe180_ecker171_removebads.csv.num2.hdf5";
		String global_kinship_file="/Volumes/Projects/DATA/arabidopsis_genomes/kinship_files_swe180_ecker171_removebads/swe180_ecker171_removebads.2.rescaled.IBS";
		int[][][] region_only1={{},{},{{3950001,3950001+100000}},{},{}};
		
		permutation(input_genotype, phenotype_full, global_kinship_file, phe_indxe, region_only1, local_gene_outfile);
		
	}
}
