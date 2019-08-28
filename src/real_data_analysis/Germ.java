package real_data_analysis;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;

import mixedmodel.KinshipMatrix;
import mixedmodel.MultiPhenotype;
import mixedmodel.Phenotype;
import mixedmodel.VariantsDouble;

public class Germ {

	public static void generate_pheno(String full, String t0, String t7, String t21){
		try{
			MultiPhenotype all=new MultiPhenotype(full);
			for(int k=0;k<all.num_of_pheno;k++){
				Phenotype p=all.phenotypes[k];
				if(p.phe_id.equals("2052_germ_284_t0")){
					p.write2file(t0);
				}else if(p.phe_id.equals("2053_germ_284_t7")){
					p.write2file(t7);
				}else if(p.phe_id.equals("2054_germ_284_t21")){
					p.write2file(t21);
				}
			}
		}catch(Exception e){e.printStackTrace();}
	}

	public static void jknf(String file, String folder, String genotype, 
			String kinship_file, String command_file){
		Phenotype p=new Phenotype(file);
		int[] nums={1,5,10};
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(command_file));
			for(int k=0;k<nums.length;k++){
				String sub_folder=folder+"/j"+nums[k]+"/";
				File the_sub_folder=new File(sub_folder);
				if(!the_sub_folder.exists())the_sub_folder.mkdir();
				p.jackknife(nums[k], sub_folder);
				//Runtime aRT = Runtime.getRuntime();
				for(int round=0;round<p.sample_ids.length;round++){
					bw.write("java -Xmx4g -jar /Users/quanlong/Documents/programs/my_jar/jawamix5.jar local " +
							" -ig "+ genotype+
							" -ip "+ sub_folder+"jackknife"+round+".tsv"+
							" -o " + sub_folder+
							" -w " +100000+
							" -ik_g "+kinship_file+";\n");
				}	
				
			}
			bw.close();
		}catch(Exception e){e.printStackTrace();}		
	}
	
	/*
	 * adding co-factor or stratifying samples using SNP 5-18570773
	 */
	public static void co_factor(String genotype_file){
		String working_folder="/Volumes/Projects/Local-kinship/germ/";
		String germ_t0=working_folder+"t0.tsv";
		String germ_t7=working_folder+"t7.tsv";
		String germ_t21=working_folder+"t21.tsv";
		VariantsDouble genotype = new VariantsDouble(genotype_file);
		int chr=5-1, loc=18570773;
		String output_folder="/Volumes/Projects/Local-kinship/germ/dog1_controled";
		try{
			String germ_file=germ_t0;
			Phenotype pheno=new Phenotype(germ_file);
			pheno.preprocess_cofactor(genotype, chr, loc, germ_file+".5-18570773.tsv");
			
			germ_file=germ_t7;
			pheno=new Phenotype(germ_file);
			pheno.preprocess_cofactor(genotype, chr, loc, germ_file+".5-18570773.tsv");
			
			germ_file=germ_t21;
			pheno=new Phenotype(germ_file);
			pheno.preprocess_cofactor(genotype, chr, loc, germ_file+".5-18570773.tsv");
		}catch(Exception e){e.printStackTrace();}	
	}
	public static void main(String[] args) {
		String working_folder="/Volumes/Projects/Local-kinship/germ/";
		String full_talbe=working_folder+"gen_arch_traits_70_plus.useful.tsv";
		String germ_t0=working_folder+"t0.tsv";
		String germ_t7=working_folder+"t7.tsv";
		String germ_t21=working_folder+"t21.tsv";
		
		String j_folder_t0=working_folder+"jackknife/t0";
		String j_folder_t7=working_folder+"jackknife/t7";
		String j_folder_t21=working_folder+"jackknife/t21";
//		generate_pheno(full_talbe, germ_t0, germ_t7, germ_t21);
		
		String genotype_file="/Volumes/Projects/DATA/arabidopsis_genomes/swe180_ecker171_removebads.csv.num2.hdf5";		
		String kinship_file="/Volumes/Projects/DATA/arabidopsis_genomes/kinship_files_swe180_ecker171_removebads/swe180_ecker171_removebads.2.rescaled.IBS";

//		jknf(germ_t0, j_folder_t0, genotype_file, kinship_file, j_folder_t0+"/command_file.txt");
//		jknf(germ_t7, j_folder_t7, genotype_file, kinship_file, j_folder_t7+"/command_file.txt");
//		jknf(germ_t21, j_folder_t21, genotype_file, kinship_file, j_folder_t21+"/command_file.txt");
		
//		co_factor(genotype_file);
		
		String all_K="/Users/quanlong/Dropbox (SH Corp)/QuanOnlyDocs/tmp/tree_local_K/kinship.324.chr3.3950001.4050000.raw.ibs";
		String sub_K=all_K+".t0.germ.raw.ibs";
		KinshipMatrix all=new KinshipMatrix(all_K);
		KinshipMatrix sub=all.submatrix(all, new Phenotype(germ_t0).sample_ids);
		sub.write2file(sub_K);
	}

}
