package real_data_analysis;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.LineNumberReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import myFileFunctions.FileFunc;
import myPlotLab.General_fucntions;
import myPlotLab.MyLineChart;

public class SummarizePlots {
	public static HashMap<Double, Integer> heritability_m=new HashMap<Double, Integer>();
	public static HashMap<Double, Integer> num_of_causal_common_m =new HashMap<Double, Integer>();
	public static HashMap<Double, Integer> concentration_rare_m=new HashMap<Double, Integer>();
	public static HashMap<Double, Integer> min_mafs_m=new HashMap<Double, Integer>();
	public static HashMap<Double, Integer> max_mafs_m=new HashMap<Double, Integer>();
	
	public static double[] heritability_a;
	public static double[] num_of_causal_common_a;
	public static double[] concentration_rare_a;
	public static double[] min_mafs_a;
	public static double[] max_mafs_a;
	
	public static boolean common=true; // if rare, then this value is false;
	
	public static void get_structure(String data_folder){
		int common_count=0, rare_count=0;
		try{
			File root=new File(data_folder);
			File[] subfolders=root.listFiles();
			System.out.println("Analyzing the structrue of the folder...");
			for(int k=0;k<subfolders.length;k++){
				if(k%1000==0)System.out.println(k/1000+"K");
				String name=subfolders[k].toString().split("/")[subfolders[k].toString().split("/").length-1];
				//System.out.println(k+":"+name);
				heritability_m.put(Double.parseDouble(name.substring(5, 8)),-1);
				int min_index=name.indexOf("min");
				int max_index=name.indexOf("max");
				min_mafs_m.put(Double.parseDouble(name.substring(min_index+3,max_index-1)),-1);
				int num_index=name.indexOf("num");
				int con_index=name.indexOf("con");
				int w_index_name=name.indexOf("w");
				if(num_index!=-1){
					common_count++;
					max_mafs_m.put(Double.parseDouble(name.substring(max_index+3,num_index-1)),-1);
					num_of_causal_common_m.put(Double.parseDouble(name.substring(num_index+3,w_index_name-1)),-1);
				}else if(con_index!=-1){
					rare_count++;
					max_mafs_m.put(Double.parseDouble(name.substring(max_index+3,con_index-1)),-1);
					concentration_rare_m.put(Double.parseDouble(name.substring(con_index+3,w_index_name-1)),-1);
				}else System.out.println("WRONG: No con or num contained in the folder name!");
			}
	
			heritability_a=FileFunc.hashmap2array_sort(heritability_m);
			num_of_causal_common_a=FileFunc.hashmap2array_sort(num_of_causal_common_m);
			concentration_rare_a=FileFunc.hashmap2array_sort(concentration_rare_m);
			min_mafs_a=FileFunc.hashmap2array_sort(min_mafs_m);
			max_mafs_a=FileFunc.hashmap2array_sort(max_mafs_m);
			for(int i=0;i<heritability_a.length;i++){
				heritability_m.put(heritability_a[i], i);
			}for(int i=0;i<num_of_causal_common_a.length;i++){
				num_of_causal_common_m.put(num_of_causal_common_a[i], i);
			}for(int i=0;i<concentration_rare_a.length;i++){
				concentration_rare_m.put(concentration_rare_a[i], i);
			}for(int i=0;i<min_mafs_a.length;i++){
				min_mafs_m.put(min_mafs_a[i], i);
			}for(int i=0;i<max_mafs_a.length;i++){
				max_mafs_m.put(max_mafs_a[i], i);
			}
			if(common_count==0 && rare_count!=0){
				SummarizePlots.common=false;
				System.out.println(common_count+"/"+rare_count);
			}else if(common_count!=0 && rare_count==0){
				SummarizePlots.common=true;
				System.out.println(common_count+"/"+rare_count);
			}else{
				System.out.println("Mixture of common and rare, perhaps WRONG!");
				System.out.println(common_count+"/"+rare_count);
			}
			System.out.println("Finished analyzing the structrue of the folder.");

		}catch(Exception e){e.printStackTrace();}
	}
	
	/*
	 * Generate data ready for plot: Common
	 */
	public static void summarize_common(String data_folder, double p_cut_emmax, double p_cut_local){
		get_structure(data_folder);
		String con_num="num";
		HashMap<Double, Integer> para_causal=num_of_causal_common_m; 
		double[][][] powers_local=new double[max_mafs_m.size()][heritability_m.size()][para_causal.size()];
		double[][][] powers_emmax=new double[max_mafs_m.size()][heritability_m.size()][para_causal.size()];;		
		double[][][] counts=new double[max_mafs_m.size()][heritability_m.size()][para_causal.size()];	
		if(data_folder.endsWith("/")){
			data_folder=data_folder.substring(0,data_folder.length()-1);
		}
		try{
			String plot_folder=data_folder+"_plot2";
			File plot_folder_file=new File(plot_folder);
			if(!plot_folder_file.exists())plot_folder_file.mkdir();
			File root=new File(data_folder);
			File[] subfolders=root.listFiles();
			BufferedWriter[][][] pred_bws=new BufferedWriter[max_mafs_m.size()][heritability_m.size()][num_of_causal_common_m.size()];
			for(int k=0;k<max_mafs_m.size();k++){
				for(int i=0;i<heritability_m.size();i++){
					for(int j=0;j<num_of_causal_common_m.size();j++){
						String the_file_folder=data_folder+"_plot2/pred/";
						File plot_folder_pred=new File(the_file_folder);
						if(!plot_folder_pred.exists())plot_folder_pred.mkdir();
						String out_file_name=the_file_folder+"m="+max_mafs_a[k]+"h="+heritability_a[i]+
								"n="+num_of_causal_common_a[j];
						pred_bws[k][i][j]=new BufferedWriter(new FileWriter(out_file_name));
						pred_bws[k][i][j].write("reg\tgblup\trgblup_sum\trgblup_beta\tEmmax\tLocal_K\tdis_LD\n");
					}		
				}
			}
			for(int k=0;k<subfolders.length;k++){
				if(k%1000==0)System.out.println(k/1000+"K");
				String name=subfolders[k].toString().split("/")[subfolders[k].toString().split("/").length-1];
				//System.out.println(k+":"+name);
				int h_index=heritability_m.get(Double.parseDouble(name.substring(5, 8)));
				int min_index=name.indexOf("min");
				int max_index=name.indexOf("max");
				//int maf_index=min_mafs.put(Double.parseDouble(name.substring(min_index+3,max_index-1)),-1);
				int num_index=name.indexOf(con_num);				
				int w_index_name=name.indexOf("w");				
				int m_index=max_mafs_m.get(Double.parseDouble(name.substring(max_index+3,num_index-1)));
				int n_index	=para_causal.get(Double.parseDouble(name.substring(num_index+3,w_index_name-1)));				
				//System.out.println(h_index+"/"+maf_index+"/"+con_index);
				File the_file_k=new File(subfolders[k]+"/summary.txt");
				if(the_file_k.exists() && (the_file_k.length()!=0)){
					BufferedReader br=new BufferedReader(new FileReader(subfolders[k]+"/summary.txt"));
					String line=br.readLine();
					double emmax_p=0, local_p=0;
					if(line!=null){
						emmax_p=Double.parseDouble(line.split("\t")[1].split("/")[0]);
						local_p=Double.parseDouble(br.readLine().split("\t")[1]);
						counts[m_index][h_index][n_index]++;
						if(emmax_p>=p_cut_emmax)							
							powers_emmax[m_index][h_index][n_index]++;
						if(local_p>=p_cut_local)
							powers_local[m_index][h_index][n_index]++;
					}br.close();
					if(new File(subfolders[k]+"/simulated_phenotype.tsv.corr.txt").exists()){
						br=new BufferedReader(new FileReader(subfolders[k]+"/simulated_phenotype.tsv.corr.txt"));
						line=br.readLine();line=br.readLine();//skip header
						double[] corrs=new double[4];
						double count=0;
						double dist_LD=0;
						if(line!=null){
							double dis_ld=Double.parseDouble(line.split("\t")[4]);
							if(!Double.isNaN(dis_ld))dist_LD=dis_ld;
						}
						while(line!=null){
							String[] tmp=line.split("\t");
							for(int i=0;i<4;i++){
								double the_value=Double.parseDouble(tmp[i]);
								if(!Double.isNaN(the_value)){
									corrs[i]+=(the_value>0)?the_value:0;
								}
							}count++;
							line=br.readLine();
						}
						if(count!=0){
							for(int i=0;i<4;i++)corrs[i]=corrs[i]/count;
							pred_bws[m_index][h_index][n_index].write(corrs[0]+"\t"+corrs[1]+"\t"+corrs[2]+"\t"+corrs[3]
									+"\t"+emmax_p+"\t"+local_p+"\t"+dist_LD+"\n");
							pred_bws[m_index][h_index][n_index].flush();
						}						
					}					
				}							
			}
			for(int k=0;k<max_mafs_m.size();k++){
				for(int i=0;i<heritability_m.size();i++){
					for(int j=0;j<num_of_causal_common_m.size();j++){					
						powers_local[k][i][j]/=counts[k][i][j];
						powers_emmax[k][i][j]/=counts[k][i][j];				
					}		
				}
			}			
			String[] legend={"Local-kinship","EMMAX"};
			plot_folder=data_folder+"_plot2/num";
			plot_folder_file=new File(plot_folder);
			if(!plot_folder_file.exists())plot_folder_file.mkdir();
			for(int k=0;k<max_mafs_a.length;k++){
				for(int i=0;i<heritability_a.length;i++){				
					String title="Heritability="+heritability_a[i]+"Max_maf="+max_mafs_a[k];
					double[] x_values=num_of_causal_common_a.clone();
					double[][] y_values=new double[2][num_of_causal_common_a.length];
					String x_lab="Number of causal";
					String y_lab="Power";
					for(int j=0;j<num_of_causal_common_a.length;j++){
						y_values[0][j]=powers_local[k][i][j]*100;
						y_values[1][j]=powers_emmax[k][i][j]*100;
					}
					General_fucntions.output_source(plot_folder+"/"+title+".png"+".source.csv", legend, x_values, y_values);
					//MyLineChart plot=new MyLineChart(title, x_lab, y_lab, legend, x_values, y_values, 500, 500, plot_folder+"/"+title+".png");
					
				}
			}			
			plot_folder=data_folder+"_plot2/heritability";
			plot_folder_file=new File(plot_folder);
			if(!plot_folder_file.exists())plot_folder_file.mkdir();
			for(int k=0;k<max_mafs_a.length;k++){
				for(int j=0;j<num_of_causal_common_a.length;j++){					
					String title="NumberCausal="+num_of_causal_common_a[j]+"Max_maf="+max_mafs_a[k];
					double[] x_values=heritability_a.clone();
					double[][] y_values=new double[2][heritability_a.length];
					String x_lab="Heritability";
					String y_lab="Power";
					for(int i=0;i<heritability_a.length;i++){
						y_values[0][i]=powers_local[k][i][j]*100;
						y_values[1][i]=powers_emmax[k][i][j]*100;
					}
					General_fucntions.output_source(plot_folder+"/"+title+".png"+".source.csv", legend, x_values, y_values);
					//MyLineChart plot=new MyLineChart(title, x_lab, y_lab, legend, x_values, y_values, 500, 500, plot_folder+"/"+title+".png");
				}		
			}
			for(int k=0;k<max_mafs_m.size();k++){
				for(int i=0;i<heritability_m.size();i++){
					for(int j=0;j<num_of_causal_common_m.size();j++){
						pred_bws[k][i][j].close();				
					}		
				}
			}
		}catch(Exception e){e.printStackTrace();}
	
	}
	
	/*
	 * Generate data ready for plot: Rare
	 */
	public static void summarize_rare(String data_folder, double p_cut_emmax, 
			double p_cut_local_s, double p_cut_agg_s, 
			double p_cut_local_e, double p_cut_agg_e, 
			double p_cut_local_l, double p_cut_agg_l){
		get_structure(data_folder);
		String con_num="con";
		HashMap<Double, Integer> para_causal_m=concentration_rare_m; 
		
		double[][][] powers_local_l=new double[max_mafs_m.size()][heritability_m.size()][para_causal_m.size()];
		double[][][] powers_emmax_l=new double[max_mafs_m.size()][heritability_m.size()][para_causal_m.size()];;		
		double[][][] powers_agg_l=new double[max_mafs_m.size()][heritability_m.size()][para_causal_m.size()];
		
		double[][][] powers_local_e=new double[max_mafs_m.size()][heritability_m.size()][para_causal_m.size()];
		double[][][] powers_emmax_e=new double[max_mafs_m.size()][heritability_m.size()][para_causal_m.size()];;		
		double[][][] powers_agg_e=new double[max_mafs_m.size()][heritability_m.size()][para_causal_m.size()];
		
		double[][][] powers_local_s=new double[max_mafs_m.size()][heritability_m.size()][para_causal_m.size()];
		double[][][] powers_emmax_s=new double[max_mafs_m.size()][heritability_m.size()][para_causal_m.size()];;		
		double[][][] powers_agg_s=new double[max_mafs_m.size()][heritability_m.size()][para_causal_m.size()];
		
		double[][][] counts=new double[max_mafs_m.size()][heritability_m.size()][para_causal_m.size()];	
		if(data_folder.endsWith("/")){
			data_folder=data_folder.substring(0,data_folder.length()-1);
		}
		try{
			String plot_folder=data_folder+"_plot2";
			File plot_folder_file=new File(plot_folder);
			if(!plot_folder_file.exists())plot_folder_file.mkdir();
			File root=new File(data_folder);
			File[] subfolders=root.listFiles();
			BufferedWriter[][][] pred_bws=new BufferedWriter[max_mafs_m.size()][heritability_m.size()][para_causal_m.size()];
			for(int k=0;k<max_mafs_m.size();k++){
				for(int i=0;i<heritability_m.size();i++){
					for(int j=0;j<para_causal_m.size();j++){
						String the_file_folder=data_folder+"_plot2/pred/";
						File plot_folder_pred=new File(the_file_folder);
						if(!plot_folder_pred.exists())plot_folder_pred.mkdir();
						String out_file_name=the_file_folder+"m="+max_mafs_a[k]+"h="+heritability_a[i]+
								"n="+concentration_rare_a[j];
						pred_bws[k][i][j]=new BufferedWriter(new FileWriter(out_file_name));
						pred_bws[k][i][j].write("reg\tgblup\trgblup_sum\trgblup_beta"
								+ "\tEmmax_L\tLocal_K_L\tAggregate_L"
								+ "\tEmmax_E\tLocal_K_E\tAggregate_E"
								+ "\tEmmax_S\tLocal_K_S\tAggregate_S"
								+ "\tdis_LD\n");
						pred_bws[k][i][j].flush();
					}		
				}
			}
			for(int k=0;k<subfolders.length;k++){
				if(k%1000==0)System.out.println(k/1000+"K");
				String name=subfolders[k].toString().split("/")[subfolders[k].toString().split("/").length-1];
				//System.out.println(k+":"+name);
				int h_index=heritability_m.get(Double.parseDouble(name.substring(5, 8)));
				int min_index=name.indexOf("min");
				int max_index=name.indexOf("max");
				//int maf_index=min_mafs.put(Double.parseDouble(name.substring(min_index+3,max_index-1)),-1);
				int num_index=name.indexOf(con_num);				
				int w_index_name=name.indexOf("w");				
				int m_index=max_mafs_m.get(Double.parseDouble(name.substring(max_index+3,num_index-1)));
				int n_index	=para_causal_m.get(Double.parseDouble(name.substring(num_index+3,w_index_name-1)));				
				//System.out.println(h_index+"/"+maf_index+"/"+con_index);
				if(new File(subfolders[k]+"/summary.txt").exists() && 
						countLines(subfolders[k]+"/summary.txt")==9){
					BufferedReader br=new BufferedReader(new FileReader(subfolders[k]+"/summary.txt"));
					String line=br.readLine();
					double emmax_pl=0, local_pl=0, agg_pl=0;
					double emmax_pe=0, local_pe=0, agg_pe=0;
					double emmax_ps=0, local_ps=0, agg_ps=0;
					if(line!=null){
						counts[m_index][h_index][n_index]++;
						// large
						emmax_pl=Double.parseDouble(line.split("\t")[1].split("/")[0]);
						local_pl=Double.parseDouble(br.readLine().split("\t")[1]);
						agg_pl=Double.parseDouble(br.readLine().split("\t")[1]);						
						if(emmax_pl>=p_cut_emmax)							
							powers_emmax_l[m_index][h_index][n_index]++;
						if(local_pl>=p_cut_local_l)
							powers_local_l[m_index][h_index][n_index]++;
						if(agg_pl>=p_cut_agg_l)
							powers_agg_l[m_index][h_index][n_index]++;
						//exact
						line=br.readLine();
						//System.out.println(line);
						//System.out.println(name);
						emmax_pe=Double.parseDouble(line.split("\t")[1].split("/")[0]);
						local_pe=Double.parseDouble(br.readLine().split("\t")[1]);
						agg_pe=Double.parseDouble(br.readLine().split("\t")[1]);						
						if(emmax_pe>=p_cut_emmax)							
							powers_emmax_e[m_index][h_index][n_index]++;
						if(local_pe>=p_cut_local_e)
							powers_local_e[m_index][h_index][n_index]++;
						if(agg_pe>=p_cut_agg_e)
							powers_agg_e[m_index][h_index][n_index]++;
						//small
						emmax_ps=Double.parseDouble(br.readLine().split("\t")[1].split("/")[0]);
						local_ps=Double.parseDouble(br.readLine().split("\t")[1]);
						agg_ps=Double.parseDouble(br.readLine().split("\t")[1]);						
						if(emmax_ps>=p_cut_emmax)							
							powers_emmax_s[m_index][h_index][n_index]++;
						if(local_ps>=p_cut_local_s)
							powers_local_s[m_index][h_index][n_index]++;
						if(agg_ps>=p_cut_agg_s)
							powers_agg_s[m_index][h_index][n_index]++;
					}br.close();
					if(new File(subfolders[k]+"/simulated_phenotype.tsv.corr.txt").exists()){
						br=new BufferedReader(new FileReader(subfolders[k]+"/simulated_phenotype.tsv.corr.txt"));
						line=br.readLine();line=br.readLine();//skip header
						double[] corrs=new double[4];
						double count=0;
						double dist_LD=0;
						if(line!=null){
							double dis_ld=Double.parseDouble(line.split("\t")[4]);
							if(!Double.isNaN(dis_ld))dist_LD=dis_ld;
						}
						while(line!=null){
							String[] tmp=line.split("\t");
							for(int i=0;i<4;i++){
								double the_value=Double.parseDouble(tmp[i]);
								if(!Double.isNaN(the_value)){
									corrs[i]+=(the_value>0)?the_value:0;
								}
							}count++;
							line=br.readLine();
						}
						if(count!=0){
							for(int i=0;i<4;i++)corrs[i]=corrs[i]/count;
							pred_bws[m_index][h_index][n_index].write(corrs[0]+"\t"+corrs[1]+"\t"+corrs[2]+"\t"+corrs[3]
									+"\t"+emmax_pl+"\t"+local_pl+"\t"+agg_pl
									+"\t"+emmax_pe+"\t"+local_pe+"\t"+agg_pe
									+"\t"+emmax_ps+"\t"+local_ps+"\t"+agg_ps
									+"\t"+dist_LD+"\n");
							pred_bws[m_index][h_index][n_index].flush();
						}						
					}					
				}							
			}
			for(int k=0;k<max_mafs_m.size();k++){
				for(int i=0;i<heritability_m.size();i++){
					for(int j=0;j<concentration_rare_m.size();j++){					
						powers_local_l[k][i][j]/=counts[k][i][j];
						powers_emmax_l[k][i][j]/=counts[k][i][j];
						powers_agg_l[k][i][j]/=counts[k][i][j];
						powers_local_e[k][i][j]/=counts[k][i][j];
						powers_emmax_e[k][i][j]/=counts[k][i][j];
						powers_agg_e[k][i][j]/=counts[k][i][j];
						powers_local_s[k][i][j]/=counts[k][i][j];
						powers_emmax_s[k][i][j]/=counts[k][i][j];
						powers_agg_s[k][i][j]/=counts[k][i][j];
					}		
				}
			}			
			String[] legend={"EMMAX_L","Local_L","Aggregate_L",
					"EMMAX_E","Local_E","Aggregate_E",
					"EMMAX_S","Local_S","Aggregate_S",};
			plot_folder=data_folder+"_plot2/con";
			plot_folder_file=new File(plot_folder);
			if(!plot_folder_file.exists())plot_folder_file.mkdir();
			for(int k=0;k<max_mafs_a.length;k++){
				for(int i=0;i<heritability_a.length;i++){				
					String title="Heritability="+heritability_a[i]+"Max_maf="+max_mafs_a[k];
					double[] x_values=concentration_rare_a.clone();
					double[][] y_values=new double[9][concentration_rare_a.length];
					String x_lab="Proportion of causal";
					String y_lab="Power";
					for(int j=0;j<concentration_rare_a.length;j++){
						y_values[0][j]=powers_emmax_l[k][i][j]*100;
						y_values[1][j]=powers_local_l[k][i][j]*100;
						y_values[2][j]=powers_agg_l[k][i][j]*100;
						y_values[3][j]=powers_emmax_e[k][i][j]*100;
						y_values[4][j]=powers_local_e[k][i][j]*100;
						y_values[5][j]=powers_agg_e[k][i][j]*100;
						y_values[6][j]=powers_emmax_s[k][i][j]*100;
						y_values[7][j]=powers_local_s[k][i][j]*100;
						y_values[8][j]=powers_agg_s[k][i][j]*100;
					}
					General_fucntions.output_source(plot_folder+"/"+title+".png"+".source.csv", legend, x_values, y_values);
					//MyLineChart plot=new MyLineChart(title, x_lab, y_lab, legend, x_values, y_values, 500, 500, plot_folder+"/"+title+".png");
					
				}
			}			
			plot_folder=data_folder+"_plot2/heritability";
			plot_folder_file=new File(plot_folder);
			if(!plot_folder_file.exists())plot_folder_file.mkdir();
			for(int k=0;k<max_mafs_a.length;k++){
				for(int j=0;j<concentration_rare_a.length;j++){					
					String title="PropCausal="+concentration_rare_a[j]+"Max_maf="+max_mafs_a[k];
					double[] x_values=heritability_a.clone();
					double[][] y_values=new double[9][heritability_a.length];
					String x_lab="Heritability";
					String y_lab="Power";
					for(int i=0;i<heritability_a.length;i++){
						y_values[0][i]=powers_emmax_l[k][i][j]*100;
						y_values[1][i]=powers_local_l[k][i][j]*100;
						y_values[2][i]=powers_agg_l[k][i][j]*100;
						y_values[3][i]=powers_emmax_e[k][i][j]*100;
						y_values[4][i]=powers_local_e[k][i][j]*100;
						y_values[5][i]=powers_agg_e[k][i][j]*100;
						y_values[6][i]=powers_emmax_s[k][i][j]*100;
						y_values[7][i]=powers_local_s[k][i][j]*100;
						y_values[8][i]=powers_agg_s[k][i][j]*100;
					}
					General_fucntions.output_source(plot_folder+"/"+title+".png"+".source.csv", legend, x_values, y_values);
					//MyLineChart plot=new MyLineChart(title, x_lab, y_lab, legend, x_values, y_values, 500, 500, plot_folder+"/"+title+".png");
				}		
			}
			for(int k=0;k<max_mafs_m.size();k++){
				for(int i=0;i<heritability_m.size();i++){
					for(int j=0;j<concentration_rare_m.size();j++){
						pred_bws[k][i][j].close();				
					}		
				}
			}
		}catch(Exception e){e.printStackTrace();}
	
	}
	
	/*
	 * Generate data ready for plot: Rare, when EMMAX doesn't work
	 */
	public static void summarize_rare_con_emmax(String data_folder, double p_cut_emmax, 
			double p_cut_local_s, double p_cut_agg_s, 
			double p_cut_local_e, double p_cut_agg_e, 
			double p_cut_local_l, double p_cut_agg_l){
		get_structure(data_folder);
		String con_num="con";
		HashMap<Double, Integer> para_causal_m=concentration_rare_m; 
		
		double[][][] powers_local_l=new double[max_mafs_m.size()][heritability_m.size()][para_causal_m.size()];
		double[][][] powers_emmax_l=new double[max_mafs_m.size()][heritability_m.size()][para_causal_m.size()];;		
		double[][][] powers_agg_l=new double[max_mafs_m.size()][heritability_m.size()][para_causal_m.size()];
		
		double[][][] powers_local_e=new double[max_mafs_m.size()][heritability_m.size()][para_causal_m.size()];
		double[][][] powers_emmax_e=new double[max_mafs_m.size()][heritability_m.size()][para_causal_m.size()];;		
		double[][][] powers_agg_e=new double[max_mafs_m.size()][heritability_m.size()][para_causal_m.size()];
		
		double[][][] powers_local_s=new double[max_mafs_m.size()][heritability_m.size()][para_causal_m.size()];
		double[][][] powers_emmax_s=new double[max_mafs_m.size()][heritability_m.size()][para_causal_m.size()];;		
		double[][][] powers_agg_s=new double[max_mafs_m.size()][heritability_m.size()][para_causal_m.size()];
		
		double[][][] counts=new double[max_mafs_m.size()][heritability_m.size()][para_causal_m.size()];	
		if(data_folder.endsWith("/")){
			data_folder=data_folder.substring(0,data_folder.length()-1);
		}
		try{
			String plot_folder=data_folder+"_plot_CE";
			File plot_folder_file=new File(plot_folder);
			if(!plot_folder_file.exists())plot_folder_file.mkdir();
			File root=new File(data_folder);
			File[] subfolders=root.listFiles();
			BufferedWriter[][][] pred_bws=new BufferedWriter[max_mafs_m.size()][heritability_m.size()][para_causal_m.size()];
			for(int k=0;k<max_mafs_m.size();k++){
				for(int i=0;i<heritability_m.size();i++){
					for(int j=0;j<para_causal_m.size();j++){
						String the_file_folder=data_folder+"_plot_CE/pred/";
						File plot_folder_pred=new File(the_file_folder);
						if(!plot_folder_pred.exists())plot_folder_pred.mkdir();
						String out_file_name=the_file_folder+"m="+max_mafs_a[k]+"h="+heritability_a[i]+
								"n="+concentration_rare_a[j];
						pred_bws[k][i][j]=new BufferedWriter(new FileWriter(out_file_name));
						pred_bws[k][i][j].write("reg\tgblup\trgblup_sum\trgblup_beta"
								+ "\tEmmax_L\tLocal_K_L\tAggregate_L"
								+ "\tEmmax_E\tLocal_K_E\tAggregate_E"
								+ "\tEmmax_S\tLocal_K_S\tAggregate_S"
								+ "\tdis_LD\n");
						pred_bws[k][i][j].flush();
					}		
				}
			}
			for(int k=0;k<subfolders.length;k++){
				if(k%1000==0)System.out.println(k/1000+"K");
				String name=subfolders[k].toString().split("/")[subfolders[k].toString().split("/").length-1];
				//System.out.println(k+":"+name);
				int h_index=heritability_m.get(Double.parseDouble(name.substring(5, 8)));
				int min_index=name.indexOf("min");
				int max_index=name.indexOf("max");
				//int maf_index=min_mafs.put(Double.parseDouble(name.substring(min_index+3,max_index-1)),-1);
				int num_index=name.indexOf(con_num);				
				int w_index_name=name.indexOf("w");				
				int m_index=max_mafs_m.get(Double.parseDouble(name.substring(max_index+3,num_index-1)));
				int n_index	=para_causal_m.get(Double.parseDouble(name.substring(num_index+3,w_index_name-1)));				
				//System.out.println(h_index+"/"+maf_index+"/"+con_index);
				if(new File(subfolders[k]+"/summary.txt").exists() && 
						countLines(subfolders[k]+"/summary.txt")==9){
					BufferedReader br=new BufferedReader(new FileReader(subfolders[k]+"/summary.txt"));
					String line=br.readLine();
					double emmax_pl=0, local_pl=0, agg_pl=0;
					double emmax_pe=0, local_pe=0, agg_pe=0;
					double emmax_ps=0, local_ps=0, agg_ps=0;
					if(line!=null){
						counts[m_index][h_index][n_index]++;
						// large
						emmax_pl=Double.parseDouble(line.split("\t")[1].split("/")[0]);
						local_pl=Double.parseDouble(br.readLine().split("\t")[1]);
						agg_pl=Double.parseDouble(br.readLine().split("\t")[1]);						
						if(emmax_pl>=p_cut_emmax)							
							powers_emmax_l[m_index][h_index][n_index]++;
						else{ // emmax doesn't work
							if(local_pl>=p_cut_local_l)
								powers_local_l[m_index][h_index][n_index]++;
							if(agg_pl>=p_cut_agg_l)
								powers_agg_l[m_index][h_index][n_index]++;
						}
						//exact
						line=br.readLine();
						//System.out.println(line);
						//System.out.println(name);
						emmax_pe=Double.parseDouble(line.split("\t")[1].split("/")[0]);
						local_pe=Double.parseDouble(br.readLine().split("\t")[1]);
						agg_pe=Double.parseDouble(br.readLine().split("\t")[1]);						
						if(emmax_pe>=p_cut_emmax)							
							powers_emmax_e[m_index][h_index][n_index]++;
						else{
							if(local_pe>=p_cut_local_e)
								powers_local_e[m_index][h_index][n_index]++;
							if(agg_pe>=p_cut_agg_e)
								powers_agg_e[m_index][h_index][n_index]++;
						}
						//small
						emmax_ps=Double.parseDouble(br.readLine().split("\t")[1].split("/")[0]);
						local_ps=Double.parseDouble(br.readLine().split("\t")[1]);
						agg_ps=Double.parseDouble(br.readLine().split("\t")[1]);						
						if(emmax_ps>=p_cut_emmax)							
							powers_emmax_s[m_index][h_index][n_index]++;
						else{	
							if(local_ps>=p_cut_local_s)
								powers_local_s[m_index][h_index][n_index]++;
							if(agg_ps>=p_cut_agg_s)
								powers_agg_s[m_index][h_index][n_index]++;
						}
					}br.close();
					if(new File(subfolders[k]+"/simulated_phenotype.tsv.corr.txt").exists()){
						br=new BufferedReader(new FileReader(subfolders[k]+"/simulated_phenotype.tsv.corr.txt"));
						line=br.readLine();line=br.readLine();//skip header
						double[] corrs=new double[4];
						double count=0;
						double dist_LD=0;
						if(line!=null){
							double dis_ld=Double.parseDouble(line.split("\t")[4]);
							if(!Double.isNaN(dis_ld))dist_LD=dis_ld;
						}
						while(line!=null){
							String[] tmp=line.split("\t");
							for(int i=0;i<4;i++){
								double the_value=Double.parseDouble(tmp[i]);
								if(!Double.isNaN(the_value)){
									corrs[i]+=(the_value>0)?the_value:0;
								}
							}count++;
							line=br.readLine();
						}
						if(count!=0){
							for(int i=0;i<4;i++)corrs[i]=corrs[i]/count;
							pred_bws[m_index][h_index][n_index].write(corrs[0]+"\t"+corrs[1]+"\t"+corrs[2]+"\t"+corrs[3]
									+"\t"+emmax_pl+"\t"+local_pl+"\t"+agg_pl
									+"\t"+emmax_pe+"\t"+local_pe+"\t"+agg_pe
									+"\t"+emmax_ps+"\t"+local_ps+"\t"+agg_ps
									+"\t"+dist_LD+"\n");
							pred_bws[m_index][h_index][n_index].flush();
						}						
					}					
				}							
			}
			for(int k=0;k<max_mafs_m.size();k++){
				for(int i=0;i<heritability_m.size();i++){
					for(int j=0;j<concentration_rare_m.size();j++){					
						powers_local_l[k][i][j]/=counts[k][i][j];
						powers_emmax_l[k][i][j]/=counts[k][i][j];
						powers_agg_l[k][i][j]/=counts[k][i][j];
						powers_local_e[k][i][j]/=counts[k][i][j];
						powers_emmax_e[k][i][j]/=counts[k][i][j];
						powers_agg_e[k][i][j]/=counts[k][i][j];
						powers_local_s[k][i][j]/=counts[k][i][j];
						powers_emmax_s[k][i][j]/=counts[k][i][j];
						powers_agg_s[k][i][j]/=counts[k][i][j];
					}		
				}
			}			
			String[] legend={"EMMAX_L","Local_L","Aggregate_L",
					"EMMAX_E","Local_E","Aggregate_E",
					"EMMAX_S","Local_S","Aggregate_S",};
			plot_folder=data_folder+"_plot2/con";
			plot_folder_file=new File(plot_folder);
			if(!plot_folder_file.exists())plot_folder_file.mkdir();
			for(int k=0;k<max_mafs_a.length;k++){
				for(int i=0;i<heritability_a.length;i++){				
					String title="Heritability="+heritability_a[i]+"Max_maf="+max_mafs_a[k];
					double[] x_values=concentration_rare_a.clone();
					double[][] y_values=new double[9][concentration_rare_a.length];
					String x_lab="Proportion of causal";
					String y_lab="Power";
					for(int j=0;j<concentration_rare_a.length;j++){
						y_values[0][j]=powers_emmax_l[k][i][j]*100;
						y_values[1][j]=powers_local_l[k][i][j]*100;
						y_values[2][j]=powers_agg_l[k][i][j]*100;
						y_values[3][j]=powers_emmax_e[k][i][j]*100;
						y_values[4][j]=powers_local_e[k][i][j]*100;
						y_values[5][j]=powers_agg_e[k][i][j]*100;
						y_values[6][j]=powers_emmax_s[k][i][j]*100;
						y_values[7][j]=powers_local_s[k][i][j]*100;
						y_values[8][j]=powers_agg_s[k][i][j]*100;
					}
					General_fucntions.output_source(plot_folder+"/"+title+".png"+".source.csv", legend, x_values, y_values);
					//MyLineChart plot=new MyLineChart(title, x_lab, y_lab, legend, x_values, y_values, 500, 500, plot_folder+"/"+title+".png");
					
				}
			}			
			plot_folder=data_folder+"_plot_CE/heritability";
			plot_folder_file=new File(plot_folder);
			if(!plot_folder_file.exists())plot_folder_file.mkdir();
			for(int k=0;k<max_mafs_a.length;k++){
				for(int j=0;j<concentration_rare_a.length;j++){					
					String title="PropCausal="+concentration_rare_a[j]+"Max_maf="+max_mafs_a[k];
					double[] x_values=heritability_a.clone();
					double[][] y_values=new double[9][heritability_a.length];
					String x_lab="Heritability";
					String y_lab="Power";
					for(int i=0;i<heritability_a.length;i++){
						y_values[0][i]=powers_emmax_l[k][i][j]*100;
						y_values[1][i]=powers_local_l[k][i][j]*100;
						y_values[2][i]=powers_agg_l[k][i][j]*100;
						y_values[3][i]=powers_emmax_e[k][i][j]*100;
						y_values[4][i]=powers_local_e[k][i][j]*100;
						y_values[5][i]=powers_agg_e[k][i][j]*100;
						y_values[6][i]=powers_emmax_s[k][i][j]*100;
						y_values[7][i]=powers_local_s[k][i][j]*100;
						y_values[8][i]=powers_agg_s[k][i][j]*100;
					}
					General_fucntions.output_source(plot_folder+"/"+title+".png"+".source.csv", legend, x_values, y_values);
					//MyLineChart plot=new MyLineChart(title, x_lab, y_lab, legend, x_values, y_values, 500, 500, plot_folder+"/"+title+".png");
				}		
			}
			for(int k=0;k<max_mafs_m.size();k++){
				for(int i=0;i<heritability_m.size();i++){
					for(int j=0;j<concentration_rare_m.size();j++){
						pred_bws[k][i][j].close();				
					}		
				}
			}
		}catch(Exception e){e.printStackTrace();}
	
	}
	
	
	public static int countLines(String filename) throws IOException {
	    LineNumberReader reader  = new LineNumberReader(new FileReader(filename));
	    int cnt = 0;
	    String lineRead = "";
	    while ((lineRead = reader.readLine()) != null) {}
	    cnt = reader.getLineNumber(); 
	    reader.close();
	    return cnt;
	}
	
	public static void make_ready4R(String root_folder){
		String[] folders={"additive_plot2", "canalization_plot2", "compensation_plot2",
				"dorminant_plot2", "epistasis_plot2",  "liability_plot2", "signepistasis_plot2"}; 
		try{
			for(int f=0;f<folders.length;f++){
				String the_folder=root_folder+folders+"/";
				if(new File(the_folder).exists()){
					
				}
			}
		}catch(Exception e){e.printStackTrace();}
	}
	
	public static void main(String[] args) {
		String data_folder="/Users/quanlong/Dropbox (SH Corp)/QuanOnlyDocs/QuanCoding/data/tmp_summary";
		String data_folder2="/Users/quanlong/Documents/projects2/predictions/simulation/1000g/2014-11-17/additive/";
		//get_structure(data_folder2);
		String data_folder3="/Users/quanlong/Documents/projects2/predictions/simulation/1000g/"
				+ "minerva/rgblup-2014-nov-common/dorminant/";
		String data_folder4="/Users/quanlong/Dropbox (SH Corp)/QuanOnlyDocs/QuanCoding/data/tmp_summary/epistasis";
		double p_cut_emmax=-Math.log10(0.00181/39000000), p_cut_local=-Math.log10(0.130715527/60000);
		summarize_common(args[0], PredictionGWAS.emmax_1kg_threshold, PredictionGWAS.local_win100k_1kg_threshold);		
//		summarize_rare(args[0], PredictionGWAS.emmax_1kg_threshold, 
//				PredictionGWAS.local_win50k_1kg_threshold, PredictionGWAS.aggre_win50k_1kg_threshold,
//				PredictionGWAS.local_win100k_1kg_threshold, PredictionGWAS.aggre_win100k_1kg_threshold,
//				PredictionGWAS.local_win200k_1kg_threshold, PredictionGWAS.aggre_win200k_1kg_threshold
//				);
//		summarize_rare_con_emmax(args[0], PredictionGWAS.emmax_1kg_threshold, 
//				PredictionGWAS.local_win50k_1kg_threshold, PredictionGWAS.aggre_win50k_1kg_threshold,
//				PredictionGWAS.local_win100k_1kg_threshold, PredictionGWAS.aggre_win100k_1kg_threshold,
//				PredictionGWAS.local_win200k_1kg_threshold, PredictionGWAS.aggre_win200k_1kg_threshold
//				);
		String plot_folder="/Users/quanlong/Dropbox (SH Corp)/QuanOnlyDocs/QuanCoding/data/plot2/rgblup-2014-nov-common/";
	}
}
