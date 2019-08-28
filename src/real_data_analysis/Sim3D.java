package real_data_analysis;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.HashSet;

import mixedmodel.CompoundAnalyzer;
import mixedmodel.VariantsDouble;
import myMathLib.Test;

public class Sim3D {

	public static double pheno_value(double x1, double x2, double effect_size, String model){
		HashSet<String> allowed_models=new HashSet<String>();
		allowed_models.add("additive");allowed_models.add("epistasis");
		if(!allowed_models.contains(model)){
			System.out.println(model+" is not allowed");return Double.NaN;
		}if(model.equals("additive")){
			double e=Test.randomNumber()*(1-effect_size);
			return (x1+x2)*effect_size+e;
		}else return Double.NaN;
	}
	public static void sim_pheno(CompoundAnalyzer ca, double[] effect_size, int[] compound_index, String pheno_output){
		if(compound_index.length!=effect_size.length){
			System.out.println("compound_index.length!=effect_size.length. RETURN!");
			return;
		}
		int num_pheno=compound_index.length;
		int sample_size=ca.variants.sample_size;
		double[][] pheno=new double[sample_size][num_pheno];
		for(int k=0;k<num_pheno;k++){
			int[] snp_r1=ca.variants.find_vars_in_maf_range_in_region(ca.compound_coordinates[compound_index[k]][0][0]-1,
					ca.compound_coordinates[compound_index[k]][0][1], ca.compound_coordinates[compound_index[k]][0][2],0.2,0.5);
			int[] snp_r2=ca.variants.find_vars_in_maf_range_in_region(ca.compound_coordinates[compound_index[k]][1][0]-1,
					ca.compound_coordinates[compound_index[k]][1][1], ca.compound_coordinates[compound_index[k]][1][2],0.2,0.5);
			System.out.println(snp_r1.length+ " SNPs in Compound"+compound_index[k]+".1");
			System.out.println(snp_r2.length+ " SNPs in Compound"+compound_index[k]+".2");
			double[] snp1=ca.variants.load_one_variant_by_index(ca.compound_coordinates[compound_index[k]][0][0]-1,snp_r1[snp_r1.length/2]);
			double[] snp2=ca.variants.load_one_variant_by_index(ca.compound_coordinates[compound_index[k]][1][0]-1,snp_r2[snp_r2.length/2]);
			for(int i=0;i<sample_size;i++)
				pheno[i][k]=pheno_value(snp1[i],snp2[i],effect_size[k],"additive");
		}
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(pheno_output));
			bw.write("ID");
			for(int k=0;k<num_pheno;k++)bw.write("\tP"+k);
			bw.write("\n");
			for(int i=0;i<sample_size;i++){
				bw.write(ca.variants.sample_ids[i]);
				for(int k=0;k<num_pheno;k++)bw.write("\t"+pheno[i][k]);
				bw.write("\n");
			}bw.close();
		}catch(Exception e){e.printStackTrace();}
	}
	
	 
	
	public static void main(String[] args) {
		String data_folder="/Users/quanlong/Documents/projects/prediction/3D/";
		String compound_coordinates_file=data_folder+"HiC/tmp.txt";
		String hdf5_1000g=data_folder+"1000g/g1k_all.hdf5";
		String hdf5_gtex=data_folder+"gtex/genotype/GTEx_imputed_dosage.hdf5";
		
		CompoundAnalyzer ca=new CompoundAnalyzer(hdf5_1000g, compound_coordinates_file);
		String pheno_output=data_folder+"sim/pheno.tsv";
		int[] compound_index={0, 1, 2, 4, 15};
		double[] effect_size={0.9,0.7,0.5,0.3,0.1};
		Sim3D.sim_pheno(ca, effect_size, compound_index, pheno_output);

	}

}
