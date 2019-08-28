package real_data_analysis;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashMap;

import mixedmodel.GBLUP;
import mixedmodel.MultiPhenotype;
import mixedmodel.Phenotype;
import mixedmodel.VariantsDouble;

public class Multiregression {

	public static void predict_multiregression_snps(String lipid_meta_snps, String map_file, String genotype, String phenotype){
		try{
			ArrayList<int[]> LDL_snps=new ArrayList<int[]>();
			ArrayList<int[]> HDL_snps=new ArrayList<int[]>();
			ArrayList<int[]> TG_snps=new ArrayList<int[]>();
			BufferedReader br_map=new BufferedReader(new FileReader(map_file));
			HashMap<String, int[]> map=new HashMap<String, int[]>();
			String line=br_map.readLine();
			while(line!=null){
				String[] temp=line.split("\t");
				int[] chr_location=new int[2];
				chr_location[0]=Integer.parseInt(temp[0]);
				chr_location[1]=Integer.parseInt(temp[3]);
				map.put(temp[1], chr_location);
				line=br_map.readLine();
			}
			BufferedReader br_lipid=new BufferedReader(new FileReader(lipid_meta_snps));
			line=br_lipid.readLine();
			while(line!=null){
				String[] temp=line.split("\t");
				if(map.containsKey(temp[1])){
					int[] loc=map.get(temp[1]);
					if(Integer.parseInt(temp[2])==loc[0] && Math.abs(Double.parseDouble(temp[3])-loc[1]/1000000.0)<3){
						String[] traits=temp[4].split(",");
						for(int i=0;i<traits.length;i++){
							if(traits[i].equals("LDL"))LDL_snps.add(loc);
							else if(traits[i].equals("HDL"))HDL_snps.add(loc);
							else if(traits[i].equals("TG"))TG_snps.add(loc);
							else if(traits[i].equals("TC"));
							else System.out.println(line);
						}
					}else{
						System.out.println("Location Wrong! "+ Integer.parseInt(temp[2])+"/"+loc[0] +";  " +
					         Double.parseDouble(temp[3]) +"/"+loc[1]/1000000);
					}					
				}
				line=br_lipid.readLine();
			}
			System.out.println(LDL_snps.size());
			System.out.println(HDL_snps.size());
			System.out.println(TG_snps.size());
			VariantsDouble snp_data=new VariantsDouble(genotype);
			MultiPhenotype pheno_data=new MultiPhenotype(phenotype);
			multiple_regression_prediction(LDL_snps, pheno_data.phenotypes[0], snp_data, 4000, 1000);
		}catch(Exception e){e.printStackTrace();}
	}
	
	public static void multiple_regression_prediction(ArrayList<int[]> markers, Phenotype pheno, 
			VariantsDouble geno, int training_sample_size, int rounds){
		double counts=0, total=0;
		double[][] corrs=new double[rounds][];
		for(int r=0;r<rounds;r++){
			corrs[r]=GBLUP.multiple_regression_prediction(markers, pheno, geno, training_sample_size);
			if(!Double.isNaN(corrs[r][0])){
				counts++;
				total+=corrs[r][0];
			}			
		}
		System.out.println(counts+": "+total/counts);
	}
	
	public static void multiple_regression_prediction(ArrayList<int[]> markers, Phenotype pheno, 
			VariantsDouble geno, int training_sample_size, int rounds, String output_file){
		
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(output_file));
			bw.write("Corr\n");
			double counts=0, total=0;
			double[][] corrs=new double[rounds][];
			for(int r=0;r<rounds;r++){
				corrs[r]=GBLUP.multiple_regression_prediction(markers, pheno, geno, training_sample_size);
				if(!Double.isNaN(corrs[r][0])){
					counts++;
					total+=corrs[r][0];
					System.out.println(counts+": "+corrs[r][0]+" || "+total/counts);
					bw.write(corrs[r][0]+"\n");
				}			
			}
			bw.close();
		}catch(Exception e){e.printStackTrace();}
		
	}
	
	public static void analyze_multi(int[][] markers_nfbc_array, String genotype_nfbc, 
			String phenotype_nfbc, String nfbc_out, int rounds, int training_sample_size, int phe_index){
		try{
			ArrayList<int[]> markers_nfbc =new ArrayList<int[]>();
			for(int k=0;k<markers_nfbc_array.length;k++)markers_nfbc.add(markers_nfbc_array[k]);			
			double counts=0, total=0;
			double[][] corrs=new double[rounds][];
			Phenotype nfbc_p=new MultiPhenotype(phenotype_nfbc).phenotypes[phe_index];
			VariantsDouble nfbc_g=new VariantsDouble(genotype_nfbc);
			BufferedWriter bw=new BufferedWriter(new FileWriter(nfbc_out));
			bw.write("Corr\n");
			for(int r=0;r<rounds;r++){
				corrs[r]=GBLUP.multiple_regression_prediction(markers_nfbc,  nfbc_p,
						nfbc_g, training_sample_size);
				if(!Double.isNaN(corrs[r][0])){
					counts++;
					total+=corrs[r][0];
					System.out.println(counts+": "+corrs[r][0]+" || "+total/counts);
					bw.write(corrs[r][0]+"\n");
				}			
			}
			bw.close();			
		}catch(Exception e){e.printStackTrace();}
	}
	
	public static void multi_reg_3data(){
		try{
			int rounds=500;
			// NFBC
//			String genotype_nfbc="/Volumes/Projects/DATA/Human_GWAS/NFBC/NFBC.num.ch1-22.csv.hdf5";
//			String phenotype_nfbc="/Volumes/Projects/DATA/Human_GWAS/NFBC/NFBC.phenot.tsv";
//			String nfbc_out="/Users/quanlong/Documents/projects2/predictions/corr/nfbc_multireg.txt";
//			int training_sample_size_nfbc=4000;
//
//			int[][] markers_nfbc_array={{1,109620053},{1,55579053},{2,21085700},
//					{11,61337211},{19,11056030},{19,50087106}};
//			analyze_multi(markers_nfbc_array, genotype_nfbc, phenotype_nfbc, nfbc_out, rounds, 
//					training_sample_size_nfbc, 0);
			//RA response 
			String ra_data_folder="/Users/quanlong/Documents/projects2/predictions/RA_data/";
			String genotye_ra=ra_data_folder+"genotypes/all.train.hdf5";
			String phenotype_ra=ra_data_folder+"phenotypes/all.tsv";
			String ra_response_out="/Users/quanlong/Documents/projects2/predictions/corr/ra_response_multireg.txt";
			int training_sample_size_ra=1600;
			int[][] markers_ra_response_array={{1,160513456}};
			analyze_multi(markers_ra_response_array, genotye_ra, phenotype_ra, ra_response_out, 
					rounds, training_sample_size_ra, 0);
			//RA severity  
			//int[][] markers_ra_severity_array={{13,110019399},{21,42692862},{6,38606611}};
			int[][] markers_ra_severity_array={{10,119613089},{11,1906815},{14,107036692}};
			String ra_severity_out="/Users/quanlong/Documents/projects2/predictions/corr/ra_severity_multireg10.11.14.txt";
			analyze_multi(markers_ra_severity_array, genotye_ra, phenotype_ra, ra_severity_out, 
					rounds,	training_sample_size_ra, 5);
		}catch(Exception e){e.printStackTrace();}
		
	}
	public static void main(String[] args) {
		String lipid_meta_snps="/Users/quanlong/Documents/projects2/predictions/LDL/meta_result_table.csv";
		String map_file="/Volumes/Projects/DATA/Human_GWAS/NFBC/NFBC.map";
		String genotype="/Volumes/Projects/DATA/Human_GWAS/NFBC/NFBC.num.ch1-22.csv.hdf5";
		String phenotype="/Volumes/Projects/DATA/Human_GWAS/NFBC/NFBC.phenot.tsv";
//		predict_multiregression_snps(lipid_meta_snps, map_file, genotype, phenotype);
		
		multi_reg_3data();
	}

}
