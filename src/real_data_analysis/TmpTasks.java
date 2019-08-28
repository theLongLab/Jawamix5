package real_data_analysis;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import mixedmodel.EMMAX;
import mixedmodel.MultiPhenotype;
import mixedmodel.RelationMatrix;
import mixedmodel.VariantsDouble;

public class TmpTasks {

	public static void correlations(String grid_folder, String transform_folder, String report_file){
		int total=0;
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(report_file));
			bw.write("Phenotype_name,Pearson,Spearman\n");
			File[] grid_files=new File(grid_folder).listFiles();
			File[] transform_files=new File(transform_folder).listFiles();
			HashSet<String> grid_files_in_set=new HashSet<String>();
			for(int i=0;i<grid_files.length;i++){
				String the_file=grid_files[i].toString();
				if(the_file.endsWith("w100000.csv"))grid_files_in_set.add(the_file);
			}
			for(int i=0;i<transform_files.length;i++){
				String the_file=transform_files[i].toString();
				if(the_file.endsWith("w100000.csv")){
					String the_corresponding=grid_folder+the_file.split("/")[the_file.split("/").length-1];
					if(grid_files_in_set.contains(the_corresponding)){
						//System.out.println(the_corresponding);
						double[] corr=corr(the_corresponding,the_file);
						String name0=the_file.split("/")[the_file.split("/").length-1];
						String name=myFileFunctions.FileFunc.mysplit(name0, '.')[2];
						bw.write(name+","+corr[0]+","+corr[1]+"\n");
						System.out.println(name+","+corr[0]+","+corr[1]);
					}else{
						System.out.println("NO corresponding file");
					}
				}
			}bw.close();
		}catch(Exception e){e.printStackTrace();}
		//System.out.println(total);
	}
	
	public static double[] corr(String the_corresponding, String the_file){
		double[] corr=new double[2];
		HashMap<String, Double> data1=load_file(the_corresponding);
		HashMap<String, Double> data2=load_file(the_file);
		ArrayList<Double> p1=new ArrayList<Double>();
		ArrayList<Double> p2=new ArrayList<Double>();
		for(String loc:data1.keySet()){
			if(data2.containsKey(loc)){
				p1.add(data1.get(loc));
				p2.add(data2.get(loc));
			}
		}
		System.out.println(p1.size()+","+p2.size()+","+data1.size()+","+data2.size());
		double[] p1array=myFileFunctions.FileFunc.arraylist2array(p1);
		double[] p2array=myFileFunctions.FileFunc.arraylist2array(p2);
		corr[0]=myMathLib.StatFuncs.correlationPearsons(p1array, p2array);
		corr[1]=myMathLib.StatFuncs.correlationSpearmans(p1array, p2array);
		return corr;
	}
	
	public static HashMap<String, Double> load_file(String file){
		HashMap<String, Double> data=new HashMap<String, Double>();
		try{
			BufferedReader br=new BufferedReader(new FileReader(file));
			String line=br.readLine();line=br.readLine(); line=br.readLine();
			while(line!=null){
				String[] temp=line.split(", ");
				if(Double.parseDouble(temp[2])<0.1)
					data.put(temp[0]+"_"+temp[1], Double.parseDouble(temp[2]));
				line=br.readLine();
			}			
		}catch(Exception e){e.printStackTrace();}
		return data;
	}
	
	public static void RA_responder_genotype_processing(String genotype_folder, String clinic){
		String[] drugs={"adalimumab","etanercept","infliximab","all"};
		ArrayList<Integer>[] indi_array=new ArrayList[drugs.length*2];
		HashMap<Integer,String>[] indi_set=new HashMap[drugs.length*2];
		for(int i=0;i<drugs.length*2;i++){
			indi_array[i]=new ArrayList<Integer>();
			indi_set[i]=new HashMap<Integer,String>();
		}
		try{
			BufferedReader br=new BufferedReader(new FileReader(clinic));
			String line=br.readLine();line=br.readLine();
			int index=0;
			while(line!=null){
				String[] tmp=line.split("\t");
				if(!tmp[1].equals("NA")){ //training
					indi_array[drugs.length-1].add(index);
					indi_set[drugs.length-1].put(index,tmp[0]);
					boolean found=false;
					for(int i=0;i<drugs.length;i++){
						if(drugs[i].equals(tmp[6])){
							indi_array[i].add(index);
							indi_set[i].put(index,tmp[0]);
							found=true;
						}						
					}if(!found && !tmp[6].equals("NA"))System.out.println("Wrong!tmp[6].equals()"+tmp[6]);
				}else{ //testing
					indi_array[2*drugs.length-1].add(index);
					indi_set[2*drugs.length-1].put(index,tmp[0]);
					boolean found=false;
					for(int i=0;i<drugs.length;i++){
						if(drugs[i].equals(tmp[6])){
							indi_array[i+drugs.length].add(index);
							indi_set[i+drugs.length].put(index,tmp[0]);
							found=true;
						}						
					}if(!found && !tmp[6].equals("NA"))System.out.println("Wrong!tmp[6].equals()"+tmp[6]);
				}
				line=br.readLine();
				index++;
			}
			BufferedWriter[] bws=new BufferedWriter[drugs.length*2];
			for(int i=0;i<drugs.length;i++){
				bws[i]=new BufferedWriter(new FileWriter(genotype_folder+drugs[i]+".train.csv"));
				bws[i+drugs.length]=new BufferedWriter(new FileWriter(genotype_folder+drugs[i]+".test.csv"));
			}for(int i=0;i<2*drugs.length;i++){
				bws[i].write("CHR,LOC");
				for(int k=0;k<indi_array[i].size();k++)
					bws[i].write(","+indi_set[i].get(indi_array[i].get(k)));
				bws[i].write("\n");
			}
			for(int chr=1;chr<=22;chr++){
				System.out.println("Chr"+chr);
				BufferedReader br_chr=new BufferedReader(new FileReader(genotype_folder+"Training_chr"+chr+".dos"));
				line=br_chr.readLine();
				while(line!=null){
					String[] tmp=line.split(" ");
					if(tmp.length-6!=2706)System.out.println("Wrong 2706!");
					if(Integer.parseInt(tmp[0])!=chr)System.out.println("Wrong CHR!");
					for(int i=0;i<2*drugs.length;i++){
						bws[i].write(chr+","+tmp[2]);
						for(int k=0;k<indi_array[i].size();k++)
							bws[i].write(","+togenotype(tmp[(indi_array[i].get(k))+6]));
						bws[i].write("\n");
					}
					line=br_chr.readLine();
				}
			}	
			for(int i=0;i<2*drugs.length;i++)bws[i].close();
		}catch(Exception e){e.printStackTrace();}
	}
	public static int togenotype(String x){
		double data=Double.parseDouble(x);
		if(data<0.5)return 0;
		else if(data>=0.5 && data < 1.5)return 1;
		else if(data >=1.5) return 2;
		else {System.out.println("Wrong: genotype NaN");return -100;}
	}
	
	public static void RA_responder_genotype_processing_final(String genotype_folder, String clinic){
		int sample_size=723;
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(genotype_folder+"geno.final.csv"));
			bw.write("CHR,LOC");
			BufferedReader br=new BufferedReader(new FileReader(clinic));
			String line=br.readLine();line=br.readLine();
			while(line!=null){
				String[] tmp=line.split(" ");
				bw.write(","+tmp[0]);
				line=br.readLine();
			}bw.write("\n");
			for(int chr=1;chr<=22;chr++){
				System.out.println("Chr"+chr);
				BufferedReader br_chr=new BufferedReader(new FileReader(genotype_folder+"Testdata_chr"+chr+".dos"));
				line=br_chr.readLine();
				while(line!=null){
					String[] tmp=line.split(" ");
					if(tmp.length-6!=sample_size)System.out.println("Wrong sample_size 723!");
					if(Integer.parseInt(tmp[0])!=chr)System.out.println("Wrong CHR!");
					bw.write(chr+","+tmp[2]);
					for(int i=0;i<sample_size;i++)
						bw.write(","+togenotype(tmp[i+6]));
					bw.write("\n");					
					line=br_chr.readLine();
				}
			}	
			bw.close();
		}catch(Exception e){e.printStackTrace();}
	}
	
	public static void RA_responder_phenotype_processing(String clinic, String output, String output_adj){
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(output));
			BufferedReader br=new BufferedReader(new FileReader(clinic));
			String line=br.readLine();
			String[] header=line.split("\t");
			bw.write(header[0]+"\t"+header[1]+"\t"+header[2]+"\t"+header[3]+"\t"+header[5]+"\t"+
					header[6]+"\t"+header[7]+"\t"+header[8]+"\t"+header[9]+"\t"+header[10]+"\t"+"\n");
			line=br.readLine();			
			while(line!=null){
				String[] tmp=line.split("\t");
				bw.write(tmp[0]+"\t"+tmp[1]+"\t"+tmp[2]);
				//Good, Intermediate,	NA, Non, Supernon
				if(tmp[3].equals("Good"))bw.write("\t"+4.0);
				else if(tmp[3].equals("Intermediate"))bw.write("\t"+3.0);
				else if(tmp[3].equals("Non"))bw.write("\t"+2.0);
				else if(tmp[3].equals("Supernon"))bw.write("\t"+1.0);
				else if(tmp[3].equals("NA"))bw.write("\tNA");
				else System.out.println("tmp[3] Wrong!");
				//DREAM BRAGGS EIRA new react TEAR	ABCoN Immunex BRASS 
				if(tmp[5].equals("DREAM"))bw.write("\t"+1.0);
				else if(tmp[5].equals("BRAGGS"))bw.write("\t"+2.0);
				else if(tmp[5].equals("EIRA"))bw.write("\t"+3.0);
				else if(tmp[5].equals("new"))bw.write("\t"+4.0);
				else if(tmp[5].equals("react"))bw.write("\t"+5.0);
				else if(tmp[5].equals("TEAR"))bw.write("\t"+6.0);
				else if(tmp[5].equals("ABCoN"))bw.write("\t"+7.0);
				else if(tmp[5].equals("Immunex"))bw.write("\t"+8.0);
				else if(tmp[5].equals("BRASS"))bw.write("\t"+9.0);
				else System.out.println("tmp[5] Wrong!");
				//adalimumab infliximab	etanercept NA
				if(tmp[6].equals("adalimumab"))bw.write("\t"+1.0);
				else if(tmp[6].equals("infliximab"))bw.write("\t"+2.0);
				else if(tmp[6].equals("etanercept"))bw.write("\t"+3.0);
				else if(tmp[6].equals("NA"))bw.write("\tNA");
				else System.out.println("tmp[6] Wrong!");
				bw.write("\t"+tmp[7]+"\t"+tmp[8]+"\t"+tmp[9]+"\t"+tmp[10]+"\n");
				line=br.readLine();
			}bw.close();			
			MultiPhenotype data=new MultiPhenotype(output);
			String tobe_adj="Response.deltaDAS";
			String[] adj0={"Age","Gender"};
			String[] adj1={"baselineDAS"};
			String[] adj2={"baselineDAS","Mtx"};
			String[] adj3={"Mtx"};
			String[] adj4={"Mtx","Drug"};
			String[] adj5={"Drug","baselineDAS"};
			String[] adj6={"Age","Gender","baselineDAS"};
			String[] adj7={"Age","Gender","baselineDAS","Mtx"};
			String[] adj8={"Age","Gender","Mtx"};
			String[] adj9={"Age","Gender","Mtx","Drug"};
			String[] adj10={"Age","Gender","Drug","baselineDAS"};
			data.adjust_cofactors(tobe_adj, adj0);
			data.adjust_cofactors(tobe_adj, adj1);
			data.adjust_cofactors(tobe_adj, adj2);
			data.adjust_cofactors(tobe_adj, adj3);
			data.adjust_cofactors(tobe_adj, adj4);
			data.adjust_cofactors(tobe_adj, adj5);
			data.adjust_cofactors(tobe_adj, adj6);
			data.adjust_cofactors(tobe_adj, adj7);
			data.adjust_cofactors(tobe_adj, adj8);
			data.adjust_cofactors(tobe_adj, adj9);
			data.adjust_cofactors(tobe_adj, adj10);
			data.write2file(output_adj);
		}catch(Exception e){e.printStackTrace();}
	}
	
	public static void Eli_processing(){
		try{
			String eli_file="/Users/quanlong/Documents/projects2/predictions/RA_data/fromEli/aTNF/rdrw-SAMs_cept_mab_ifx_ada.pheno";
			String folder="/Users/quanlong/Documents/projects2/predictions/RA_data/genotypes/";
			String pheno_folder="/Users/quanlong/Documents/projects2/predictions/RA_data/phenotypes/";
			String[] types={"etanercept","infliximab","adalimumab","all"};
			for(int k=0;k<types.length;k++){				
				String data_pred_hdf5=folder+types[k]+".test.hdf5";
				String output=pheno_folder+"eli/"+types[k]+".eli.test.tsv"; 
				VariantsDouble geno=new VariantsDouble(data_pred_hdf5);
				String[] ids=geno.sample_ids;
				HashMap<String, Integer> ids_map=geno.sample_id2index;
				BufferedWriter bw=new BufferedWriter(new FileWriter(output));
				bw.write("ID\t"+types[k]+".eli\n");
				BufferedReader br=new BufferedReader(new FileReader(eli_file));
				String line=br.readLine();
				while(line!=null){
					String tmp[] =line.split(" ");
					if(ids_map.containsKey(tmp[0])){
						bw.write(tmp[0]+"\t"+tmp[2]+"\n");
					}
					line=br.readLine();
				}bw.close();
			}
		}catch(Exception e){e.printStackTrace();}
	}
	
	public static void remove_K(String data_folder){
		try{
			File root=new File(data_folder);
			File[] subfolders=root.listFiles();
			for(int k=0;k<subfolders.length;k++){
				File[] the_files=subfolders[k].listFiles();
				for(int i=0;i<the_files.length;i++){
					String name=the_files[i].toString().split("/")[the_files[i].toString().split("/").length-1];
					if(!(name.equals("causal.0.txt")||name.equals("summary.txt")
							||name.equals("simulated_phenotype.tsv.corr.txt")))
						the_files[i].delete();
				}
			}
		}catch(Exception e){e.printStackTrace();}
	}
	
	public static void main(String[] args) {
		String in ="/Volumes/Projects/DATA/GETx/december_release/JAWAMix5_format/GTEx_5M_191_Dec2012.vcf.csv.all.IBS.kinship.raw.corrected.ibs";
		String out="/Volumes/Projects/DATA/GETx/december_release/JAWAMix5_format/GTEx_5M_191_Dec2012.vcf.csv.all.IBS.rescaled.kinship";
//		VariantsDouble.re_scale_kinship_matrix(in, out);
//		
//		double[] x={2,0,1,0,1,0};
//		System.out.println(EMMAX.mafc(x));
//		
//		String folder="/Volumes/Projects/DATA/Seaside/";
//		String input_geno=folder+"Total_152samples_gt_report.txt.csv.hdf5";
//		String input_pheno=folder+"AS208_EFF_MTS.adjusted.durg.tsv";
//		String kinship=folder+"kinship_files/Total_152samples_gt_report.txt.csv.hdf5.rescaled.IBS";
//		String output_folder="/Volumes/Projects/seaside/analysis/drug_adj/";
//		int round=1;
//		int min_sample_size=10;
//		int phe_index=7;
//		double p_after_corr=1000, maf_threshold_plot=0.05;
//		EMMAX.emmax_analysis(input_geno, input_pheno, kinship, output_folder, round, p_after_corr, min_sample_size, 
//				phe_index, maf_threshold_plot, true);
		
		String grid_folder="/Volumes/Projects/backup_Nov2012/Projects/GeneMappingMethodology/LocalGlobal/local_results/DNA_local/100k/";
		String transform_folder="/Volumes/Projects/Local-kinship/Phenotype_table_Seq4/";
		String report="/Volumes/Projects/Local-kinship/stat/sig_corr.csv";
//		correlations(grid_folder, transform_folder,report);
		
		String folder="/Users/quanlong/Documents/projects2/predictions/RA_data/";
		String genotype_folder=folder+	"DREAM_RA_Responders_DosageData/";
		String pheno_folder=folder+"phenotypes/";
		String clinic= folder+"phenotypes/RA_Pheno_All.tsv";
		String pheno_all_num= folder+"phenotypes/RA_Pheno.num.tsv";
		String pheno_all_num_adj= folder+"phenotypes/RA_Pheno.ADJ.tsv";
//		RA_responder_genotype_processing(genotype_folder, clinic);
//		RA_responder_phenotype_processing(clinic, pheno_all_num, pheno_all_num_adj);
//		Eli_processing();
		
		String final_folder="/Users/quanlong/Documents/projects2/predictions/RA_data/test_final/";
//		RA_responder_genotype_processing_final(final_folder, final_folder+"TestData_Cov_Release.txt");
		remove_K(args[0]);
	}

}
