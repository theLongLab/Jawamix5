package real_data_analysis;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.Arrays;

import mixedmodel.AssociationResults;
import mixedmodel.CausalRegions;
import mixedmodel.GBLUP;
import mixedmodel.KinshipMatrix;
import mixedmodel.MultiPhenotype;
import mixedmodel.Phenotype;
import mixedmodel.RelationMatrix;
import mixedmodel.SearchAndCalculateRegionalKinships;
import mixedmodel.VariantsDouble;

public class Prediction_RA {

	public static double[][] var_comp_22=new double[4][22];
	public static double[][] var_comp_22_g=new double[4][23];
	
	public static void old_codes0(){
		int win=100000;
		int[][][] regions={
				//"etanercept"
				{{0,160500001,160500001+win},{10,24350001,24350001+win}},
				//"infliximab"
				{{15,26150001,26150001+win},{11,19950001,19950001+win}},
				// "adalimumab"
				{{16,25500001,25500001+win},{2,82650001,82650001+win}},
				// "all"
				{{12,95250001,95250001+win},{13,101500001,101500001+1}}
		};
		double[][] var_comp={//"etanercept"
				{0.0764,0.0754},
				//"infliximab"
				{0.0630,0.0667},
				// "adalimumab"
				{0.0946,0.0373},
				// "all"
				{0.0516,0.0470}
		};
		String eli_file="/Users/quanlong/Documents/projects2/predictions/RA_data/fromEli/aTNF/rdrw-SAMs_cept_mab_ifx_ada.pheno";
		String geno_folder="/Users/quanlong/Documents/projects2/predictions/RA_data/genotypes/";
		String pheno_folder="/Users/quanlong/Documents/projects2/predictions/RA_data/phenotypes/";
		String local_K_folder=geno_folder+"local_RRM/";
		
		String[] types={"etanercept","infliximab","adalimumab","all"};
		double[] h2={0.4, 0.4, 0.4, 0.4};
		
		
		for(int k=0;k<types.length;k++){	
			System.out.println("==============="+types[k]+"=================");
			Phenotype obs_pheno=new MultiPhenotype(pheno_folder+types[k]+".tsv").phenotypes[0];
			Phenotype comparison=new MultiPhenotype(pheno_folder+"eli/"+types[k]+".eli.test.tsv").phenotypes[0]; 
			VariantsDouble obs_geno=new VariantsDouble(geno_folder+types[k]+".train.hdf5");
			VariantsDouble pred_geno=new VariantsDouble(geno_folder+types[k]+".test.hdf5");
			KinshipMatrix obs_matrix=new KinshipMatrix(geno_folder+types[k]+".train.K.RRM");
			RelationMatrix relation_obs_pred=new RelationMatrix(geno_folder+types[k]+".train.test.K.RRM");
			int[][] the_regions=regions[k];
			String[] reginal_kinship_files=new String[the_regions.length];
			String[] relation_matrix_files=new String[the_regions.length];
			for(int i=0;i<the_regions.length;i++){
				reginal_kinship_files[i]=local_K_folder+types[k]+".train."+(the_regions[i][0]+1)+"."+
						the_regions[i][1]+"."+the_regions[i][2]+".K.RRM";
				relation_matrix_files[i]=local_K_folder+types[k]+".train.test."+(the_regions[i][0]+1)+"."+
						the_regions[i][1]+"."+the_regions[i][2]+".K.RRM";
			}
			
			double[] beta0_array={0,1}, beta1_array={0,1};
			for(int j=0;j<beta1_array.length;j++){
				for(int i=0;i<beta0_array.length;i++){
					double beta0=beta0_array[i];				
					double beta1=beta1_array[j];
					double[] the_new_comp=var_comp[k].clone();
					the_new_comp[0]=the_new_comp[0]*beta0;
					the_new_comp[1]=the_new_comp[1]*beta1;
					CausalRegions causal=new CausalRegions(the_regions, null, reginal_kinship_files, relation_matrix_files,  
							obs_matrix, relation_obs_pred);
					GBLUP gblup=new GBLUP(obs_geno, pred_geno, obs_pheno, obs_matrix, relation_obs_pred, h2[k]);
					//gblup.estimate_local(causal, the_new_comp);
					gblup.estimate_global_local_sum(causal, the_new_comp);
					double[] corr=MultiPhenotype.correlation(comparison, gblup.pred_pheno);
					System.out.println(beta0+":"+var_comp[k][0]+"/"+beta1+":"+var_comp[k][1]+" = "+corr[0]+"\t"+corr[1]+"\t"+corr[2]+"\t"+corr[3]);
				}
			}			
		}
	}
	
	public static void old_RA(){
		int win=100000;
		int[][][] regions={
				//"etanercept"
				{{0,160500001,160500001+win},{10,24350001,24350001+win}},
				//"infliximab"
				{{15,26150001,26150001+win},{11,19950001,19950001+win}},
				// "adalimumab"
				{{16,25500001,25500001+win},{2,82650001,82650001+win}},
				// "all"
				{{12,95250001,95250001+win},{13,101500001,101500001+1}}
		};
//		double[][] var_comp={//"etanercept"
//				{0.0764,0.0754},
//				//"infliximab"
//				{0.0630,0.0667},
//				// "adalimumab"
//				{0.0946,0.0373},
//				// "all"
//				{0.0516,0.0470}
//		};
		String eli_file="/Users/quanlong/Documents/projects2/predictions/RA_data/fromEli/aTNF/rdrw-SAMs_cept_mab_ifx_ada.pheno";
		String geno_folder="/Users/quanlong/Documents/projects2/predictions/RA_data/genotypes/";
		String pheno_folder="/Users/quanlong/Documents/projects2/predictions/RA_data/phenotypes/";
		String result_folder="/Users/quanlong/Documents/projects2/predictions/RA_data/result/";
		String local_K_folder=geno_folder+"local_RRM/";
		
		String[] types={"etanercept","infliximab","adalimumab","all"};
		double[] h2={0.4, 0.4, 0.4, 0.4};
		
		for(int k=0;k<types.length;k++){	
			System.out.println("==============="+types[k]+"=================");
			String local_kinship_result_file=result_folder+types[k]+"2/Local_VO.0.Response.deltaDAS.w100000.csv";
			regions[k]=generate_regions(local_kinship_result_file, var_comp_22[k]);
			var_comp_22_g[k][0]=h2[k];
			for(int i=0;i<22;i++)var_comp_22_g[k][i+1]=var_comp_22[k][i];
			Phenotype obs_pheno=new MultiPhenotype(pheno_folder+types[k]+".tsv").phenotypes[0];
			
			Phenotype comparison=new MultiPhenotype(pheno_folder+"eli/"+types[k]+".eli.test.tsv").phenotypes[0]; 
			VariantsDouble obs_geno=new VariantsDouble(geno_folder+types[k]+".train.hdf5");
			VariantsDouble pred_geno=new VariantsDouble(geno_folder+types[k]+".test.hdf5");
			KinshipMatrix obs_matrix=new KinshipMatrix(geno_folder+types[k]+".train.K.RRM");
			RelationMatrix relation_obs_pred=new RelationMatrix(geno_folder+types[k]+".train.test.K.RRM");
			int[][] the_regions=regions[k];
			String[] reginal_kinship_files=new String[the_regions.length];
			String[] relation_matrix_files=new String[the_regions.length];
			for(int i=0;i<the_regions.length;i++){
				reginal_kinship_files[i]=local_K_folder+types[k]+".train."+(the_regions[i][0]+1)+"."+
						the_regions[i][1]+"."+the_regions[i][2]+".K.RRM";
				relation_matrix_files[i]=local_K_folder+types[k]+".train.test."+(the_regions[i][0]+1)+"."+
						the_regions[i][1]+"."+the_regions[i][2]+".K.RRM";
			}			
		
			CausalRegions causal=new CausalRegions(the_regions, null, reginal_kinship_files, relation_matrix_files,  
							obs_matrix, relation_obs_pred);//, obs_geno, pred_geno);
			GBLUP gblup=new GBLUP(obs_geno, pred_geno, obs_pheno, obs_matrix, relation_obs_pred, h2[k]);
			gblup.estimate_local_reg(causal, var_comp_22[k]);
			//gblup.estimate_local(causal, var_comp_22[k]);
			//gblup.estimate_local_no_global(causal, var_comp_22[k]);
			//gblup.estimate_global_local_sum(causal, var_comp_22[k]);
			//gblup.estimate_local_only_sum(causal, var_comp_22[k]);
			double[] corr=MultiPhenotype.correlation(comparison, gblup.pred_pheno);
			System.out.println(corr[0]+"\t"+corr[1]+"\t"+corr[2]+"\t"+corr[3]);
//			for(int i=0;i<gblup.pred_pheno.sample_ids.length;i++){
//				//System.out.println(gblup.pred_pheno.sample_ids[i]+","+gblup.pred_pheno.values[i]+","
//							//+gblup.pred_pheno.values[i]+"\n");
//			}
			
			
		}
	}
	public static int[][] generate_regions(String local_kinship_result_file, double[] var_comp){
		int num_chr=22,  win=100000;
		AssociationResults local=new AssociationResults(local_kinship_result_file, 0.9);
		int[][] regions=new int[num_chr][3];
		double[] ps=new double[num_chr];
		//double[] r2s=new double[num_chr];
		Arrays.fill(ps, 1); Arrays.fill(var_comp, 1); 
		for(int chr=0;chr<num_chr;chr++){
			regions[chr][0]=chr;
		}
		for(int index=0;index<local.num_of_var;index++){
			if(local.pvalue[index]<ps[local.chr[index]-1]){
				ps[local.chr[index]-1]=local.pvalue[index];
				var_comp[local.chr[index]-1]=local.AdjustedR2[index];//*(-Math.log(local.pvalue[index]));
				regions[local.chr[index]-1][1]=local.location[index];
				regions[local.chr[index]-1][2]=local.location[index]+win;
			}
		}
		for(int chr=0;chr<num_chr;chr++){
			System.out.println(regions[chr][0]+":"+regions[chr][1]+":"+regions[chr][2]+":"+var_comp[chr]+":");
		}
		return regions;
	}
	
	public static void nfbc(){
		try{
			String output_file="/Users/quanlong/Documents/projects2/predictions/LDL/gblup_localblup.txt";
			BufferedWriter bw=new BufferedWriter(new FileWriter(output_file));
			bw.write("globle\tlocal_sum\tlocal_beta\n");
			int win=100000, peak_width=10000000;
			double logged_pcutpoff=5.0;
			String data_folder="/Volumes/Projects/DATA/Human_GWAS/NFBC/";
			String genotye_file=data_folder+"NFBC.num.ch1-22.csv.hdf5";
			String phenotype_file=data_folder+"NFBC.phenot.tsv";
			String full_kinship_file=genotye_file+".RRM";
			//String analysis_folder="/Users/quanlong/Documents/projects2/predictions/LDL/";
			String local_k_gwas="/Volumes/Projects/Local-kinship/NFBC/Local_VO.0.ldlres.w100000.csv";
			Phenotype full_pheno=new MultiPhenotype(phenotype_file).phenotypes[0];
			KinshipMatrix kinship_full=new KinshipMatrix(full_kinship_file);
			VariantsDouble geno_full=new VariantsDouble(genotye_file);
			String local_K_folder="/Users/quanlong/Documents/projects2/predictions/LDL/local_Ks/";
			
			double h2=0.45;
			int obs_sample_size=4500;
			double g_record=0, l_record=0, ll_record=0;
			for(int round=0;round<100;round++){
				System.out.println("==== Round: "+(round+1)+"======");
				String[] obs_ids=GBLUP.select_ids_random(obs_sample_size, full_pheno, geno_full);
				SearchAndCalculateRegionalKinships search_regions=new SearchAndCalculateRegionalKinships(local_k_gwas, win, 
					peak_width, logged_pcutpoff, local_K_folder,null);//, new VariantsDouble(genotye_file));
				double[] corrs=GBLUP.test_batch(full_pheno, geno_full, kinship_full, obs_ids, 
						h2, search_regions.regions, search_regions.var_comp, 
					search_regions.reginal_full_kinship_files, null, 3000);
				g_record+=corrs[0];
				l_record+=corrs[1];
				ll_record+=corrs[2];
				System.out.println("===="+g_record+"/"+l_record+"/"+ll_record+"/"+(round+1)+"======");
				bw.write(corrs[0]+"\t"+corrs[1]+"\t"+corrs[2]+"\n");
			}bw.close();
			//GBLUP.test_batch(full_pheno, kinship_full, obs_ids, h2, null);
		}catch(Exception e){e.printStackTrace();}
		
	}
	
	public static void seaside(){
		double g_record=0;
		double l_record=0, ll_record=0;
		int win=1000000, peak_width=10000000;
		double logged_pcutpoff=8.5;
		String data_folder="/Volumes/Projects/DATA/Seaside/";
		String genotye_file=data_folder+"Total_152samples_gt_report.txt.csv.hdf5";
		String phenotype_file=data_folder+"AS208_EFF_MTS.durg.tsv";
		String full_kinship_file=genotye_file+".RRM";
		String local_k_gwas="/Volumes/Projects/seaside/analysis/drug/Local_VO.4.EXPRERS.w500000.csv";
		Phenotype full_pheno=new MultiPhenotype(phenotype_file).phenotypes[4];
		KinshipMatrix kinship_full=new KinshipMatrix(full_kinship_file);
		VariantsDouble geno_full=new VariantsDouble(genotye_file);
		String local_K_folder="/Users/quanlong/Documents/projects2/predictions/seaside/local_Ks/";
		double h2=0.5;
		int obs_sample_size=40;
		for(int round=0;round<100;round++){
			String[] obs_ids=GBLUP.select_ids_random(obs_sample_size, full_pheno, geno_full);
			SearchAndCalculateRegionalKinships search_regions=new SearchAndCalculateRegionalKinships(local_k_gwas, win, 
					peak_width, logged_pcutpoff, local_K_folder, geno_full);
			double[] corrs=GBLUP.test_batch(full_pheno, geno_full, kinship_full, obs_ids, h2, search_regions.regions, search_regions.var_comp, 
					search_regions.reginal_full_kinship_files, null,30);
			g_record+=corrs[0];
			l_record+=corrs[1];
			ll_record+=corrs[2];
			System.out.println("===="+g_record+"/"+l_record+"/"+ll_record+"/"+(round+1)+"======");
		}
		System.out.println("====\n"+g_record+"/"+l_record);
		//GBLUP.test_batch(full_pheno, kinship_full, obs_ids, h2, null);
	}
	
	public static void ra_drug(){
		try{
			String output_file="/Users/quanlong/Documents/projects2/predictions/RA_data/gblup_localblup.txt";
			BufferedWriter bw=new BufferedWriter(new FileWriter(output_file));
			bw.write("globle\tlocal_sum\tlocal_beta\n");
			int win=100000, peak_width=10000000;
			double logged_pcutpoff=5.0;
			String data_folder="/Users/quanlong/Documents/projects2/predictions/RA_data/";
			String genotye_file=data_folder+"genotypes/all.train.hdf5";
			String phenotype_file=data_folder+"phenotypes/all.tsv";
			String full_kinship_file=data_folder+"genotypes/all.train.K.RRM";
			//String analysis_folder="/Users/quanlong/Documents/projects2/predictions/LDL/";
			String local_k_gwas=data_folder+"result/etanercept2/Local_VO.0.Response.deltaDAS.w100000.csv";
			Phenotype full_pheno=new MultiPhenotype(phenotype_file).phenotypes[0];
			KinshipMatrix kinship_full=new KinshipMatrix(full_kinship_file);
			VariantsDouble geno_full=new VariantsDouble(genotye_file);
			String local_K_folder=data_folder+"genotypes/local_RRM/";
			
			double h2=0.05;
			int obs_sample_size=1600;
			double g_record=0, l_record=0, ll_record=0;
			for(int round=0;round<100;round++){
				System.out.println("==== Round: "+(round+1)+"======");
				String[] obs_ids=GBLUP.select_ids_random(obs_sample_size, full_pheno, geno_full);
				SearchAndCalculateRegionalKinships search_regions=new SearchAndCalculateRegionalKinships(local_k_gwas, win, 
					peak_width, logged_pcutpoff, local_K_folder,null);//, new VariantsDouble(genotye_file));
				double[] corrs=GBLUP.test_batch(full_pheno, geno_full, kinship_full, obs_ids, 
						h2, search_regions.regions, search_regions.var_comp, 
					search_regions.reginal_full_kinship_files, null, 1000);
				g_record+=corrs[0];
				l_record+=corrs[1];
				ll_record+=corrs[2];
				System.out.println("===="+g_record+"/"+l_record+"/"+ll_record+"/"+(round+1)+"======");
				bw.write(corrs[0]+"\t"+corrs[1]+"\t"+corrs[2]+"\n");
			}bw.close();
			//GBLUP.test_batch(full_pheno, kinship_full, obs_ids, h2, null);
		}catch(Exception e){e.printStackTrace();}
		
	}
	
	public static void ra_severity(double logged_pcutpoff){
		try{
			String output_file="/Users/quanlong/Documents/projects2/predictions/RA_data/severity_gblup_localblup." +
					+logged_pcutpoff+".txt";
			BufferedWriter bw=new BufferedWriter(new FileWriter(output_file));
			bw.write("globle\tlocal_sum\tlocal_beta\n");
			int win=100000, peak_width=10000000;
			//double logged_pcutpoff=22.0;
			String data_folder="/Users/quanlong/Documents/projects2/predictions/RA_data/";
			String genotye_file=data_folder+"genotypes/all.train.hdf5";
			String phenotype_file=data_folder+"phenotypes/all.tsv";
			String full_kinship_file=data_folder+"genotypes/all.train.K.RRM";
			//String analysis_folder="/Users/quanlong/Documents/projects2/predictions/LDL/";
			String local_k_gwas=data_folder+"result/all2/Local_VO.5.baselineDAS.w100000.csv";
			Phenotype full_pheno=new MultiPhenotype(phenotype_file).phenotypes[5]; // baselineDAS
			KinshipMatrix kinship_full=new KinshipMatrix(full_kinship_file);
			VariantsDouble geno_full=new VariantsDouble(genotye_file);
			String local_K_folder=data_folder+"genotypes/local_RRM/";
			
			double h2=0.05;
			int obs_sample_size=1600;
			double g_record=0, l_record=0, ll_record=0;
			for(int round=0;round<100;round++){
				System.out.println("==== Round: "+(round+1)+"======");
				String[] obs_ids=GBLUP.select_ids_random(obs_sample_size, full_pheno, geno_full);
				SearchAndCalculateRegionalKinships search_regions=new SearchAndCalculateRegionalKinships(local_k_gwas, win, 
					peak_width, logged_pcutpoff, local_K_folder, null);//new VariantsDouble(genotye_file));
				double[] corrs=GBLUP.test_batch(full_pheno, geno_full, kinship_full, obs_ids, 
						h2, search_regions.regions, search_regions.var_comp, 
					search_regions.reginal_full_kinship_files, null, 1000);
				g_record+=corrs[0];
				l_record+=corrs[1];
				ll_record+=corrs[2];
				System.out.println("===="+g_record+"/"+l_record+"/"+ll_record+"/"+(round+1)+"======");
				bw.write(corrs[0]+"\t"+corrs[1]+"\t"+corrs[2]+"\n");
			}bw.close();
			//GBLUP.test_batch(full_pheno, kinship_full, obs_ids, h2, null);
		}catch(Exception e){e.printStackTrace();}
		
	}
	
	public static void RA_final_global(String predicting_clinic, String predicting_geno, String predicting_relationship,
			String output_quanti, String output_binary){
		String[][] adj_try={
				{"baselineDAS","Age","Gender","Mtx"},
				{"baselineDAS","Gender","Mtx"},
				{"baselineDAS","Age","Mtx"},
				{"baselineDAS","Age","Gender"},
				{"baselineDAS","Age",},
				{"baselineDAS","Gender"},
				{"baselineDAS","Mtx"},
				{"baselineDAS"},
				{"Age","Gender","Mtx"},
				{"Gender","Mtx"},
				{"Age","Mtx"},
				{"Age","Gender"},
				{"Age"},
				{"Gender"},
				{"Mtx"},
				{}			
		};
		String[][] adj={				
				{"baselineDAS","Age","Gender"}						
		};
		String tobe_adj="Response.deltaDAS";
		
		String folder="/Users/quanlong/Documents/projects2/predictions/RA_data/";
		String training_file=folder+"test_final/all.tsv";
		String training_geno=folder+"genotypes/all.train.hdf5";
		String training_kinship=folder+"genotypes/all.train.K.RRM";
		
		MultiPhenotype training=new MultiPhenotype(training_file); 
		MultiPhenotype predicting=new MultiPhenotype(predicting_clinic);
		VariantsDouble obs_geno=new VariantsDouble(training_geno);
		VariantsDouble pred_geno=new VariantsDouble(predicting_geno);
		Phenotype obs_pheno=training.phenotypes[0];
		KinshipMatrix obs_matrix=new KinshipMatrix(training_kinship);
		RelationMatrix relation_obs_pred=new RelationMatrix(predicting_relationship);
		double h2=0.3;//		

		GBLUP gblup=new GBLUP(obs_geno, pred_geno, obs_pheno, obs_matrix, relation_obs_pred, h2);
		gblup.estimate_global();
		// storing results
		String[] ids0=gblup.pred_pheno.sample_ids.clone();
		double[] values_noclinic=gblup.pred_pheno.values.clone();
		adjust(values_noclinic, obs_pheno.values);

		Phenotype comparison=new MultiPhenotype("/Users/quanlong/Documents/projects2/predictions/RA_data/test_final/all.eli.test.tsv").phenotypes[0];
		for(int adj_index=0; adj_index<adj.length;adj_index++){
			double[] betas=training.adjust_cofactors(tobe_adj, adj[adj_index]);
			obs_pheno=training.phenotypes[training.phenotypes.length-1];
			gblup=new GBLUP(obs_geno, pred_geno, obs_pheno, obs_matrix, relation_obs_pred, h2);
			gblup.estimate_global();
			// again, check results
			String[] ids1=gblup.pred_pheno.sample_ids.clone();
			double[] values_clinic=gblup.pred_pheno.values.clone();
			adjust(values_clinic, obs_pheno.values);
			//add_back
			for(int k=0;k<values_clinic.length;k++){
				values_clinic[k]=values_clinic[k]+betas[0];
				for(int i=0;i<adj[adj_index].length;i++){
					Phenotype the_tmp_pheno=predicting.phenotypes[predicting.phe_id2index.get(adj[adj_index][i])];
					if(the_tmp_pheno.sample_id2index.containsKey(ids1[k])){
						int sample_index=the_tmp_pheno.sample_id2index.get(ids1[k]);
						values_clinic[k]=values_clinic[k]+betas[i+1]*
								the_tmp_pheno.values[sample_index];
					}else{
						values_clinic[k]=values_clinic[k]+betas[i+1]*the_tmp_pheno.sample_mean();
					}				
				}
			}
			try{
				String sep=",";
				BufferedWriter bw1=new BufferedWriter(new FileWriter(output_quanti));
				bw1.write("ID"+sep+"pred_gen_facs"+sep+"pred_clin_gen_facs\n");
				BufferedWriter bw2=new BufferedWriter(new FileWriter(output_binary));
				bw2.write("ID"+sep+"belief_gen"+sep+"belief_clin_gen\n");

				for(int k=0;k<ids0.length;k++){
					if(!ids0[k].equals(ids1[k]))System.out.println("Wrong");
					bw1.write(ids0[k]+sep+values_noclinic[k]+sep+values_clinic[k]+"\n");
					bw2.write(ids0[k]+sep+(nonresponder(values_noclinic[k], predicting.phenotypes[0].values[k])?1:0)
							+sep+(nonresponder(values_clinic[k], predicting.phenotypes[0].values[k])?1:0)+"\n");
				}
				bw1.close();bw2.close();
			}catch(Exception e){e.printStackTrace();}
			MultiPhenotype reuslt=new MultiPhenotype(output_quanti);
			double[] corr=MultiPhenotype.correlation(comparison, reuslt.phenotypes[0]);
			System.out.println("===="+adj_index+"======"+h2);
			System.out.println("non:  "+corr[0]+"\t"+corr[1]+"\t"+corr[2]+"\t"+corr[3]);
			corr=MultiPhenotype.correlation(comparison, reuslt.phenotypes[1]);
			System.out.println("gen:  "+corr[0]+"\t"+corr[1]+"\t"+corr[2]+"\t"+corr[3]);

		}
	}
	
	
	public static void RA_final_local(String predicting_clinic, String predicting_geno, String predicting_relationship,
			String output_quanti, String output_binary){
		
		int win=100000;
		int[][][] regions={
				//set 1
				{{0,160500001,160500001+win},
				{13,105700001,105700001+win},
				{13, 75650001, 75650001+win}, 
				{11, 14400001, 14400001+win}},
				//set 2				
				{{13,105700001,105700001+win},
				{13, 75650001, 75650001+win}, 
				{11, 14400001, 14400001+win}},
				//set 3
				{{0,160500001,160500001+win},				
				{13, 75650001, 75650001+win}, 
				{11, 14400001, 14400001+win}},
				//set 4
				{{0,160500001,160500001+win},
				{13,105700001,105700001+win},				
				{11, 14400001, 14400001+win}},
				//set 5
				{{0,160500001,160500001+win},
				{13,105700001,105700001+win},
				{13, 75650001, 75650001+win}},
				//set 6
				{{13, 75650001, 75650001+win}, 
				{11, 14400001, 14400001+win}},
				//set 7
				{{13,105700001,105700001+win},
				{11, 14400001, 14400001+win}},
				//set 8
				{{13,105700001,105700001+win},
				{13, 75650001, 75650001+win}}, 				
				//set 9
				{{0,160500001,160500001+win},
				{11, 14400001, 14400001+win}},
				//set 10
				{{0,160500001,160500001+win},				
				{13, 75650001, 75650001+win}},
				//set 11
				{{0,160500001,160500001+win},
				{13,105700001,105700001+win}},
				//set 12
				{{0,160500001,160500001+win}},
				//set 13
				{{13,105700001,105700001+win}},
				//set 14
				{{13, 75650001, 75650001+win}}, 
				//set 15
				{{11, 14400001, 14400001+win}},
		};
		//double[][] var_comp={{0.05,0.05,0.05,0.05}};
		String[] adj={"baselineDAS","Age","Gender"};
		String tobe_adj="Response.deltaDAS";
		
		String folder="/Users/quanlong/Documents/projects2/predictions/RA_data/";
		String local_K_folder=folder+"test_final/local_K/";
		String training_file=folder+"test_final/all.tsv";
		String training_geno=folder+"genotypes/all.train.hdf5";
		String training_kinship=folder+"genotypes/all.train.K.RRM";

		MultiPhenotype training=new MultiPhenotype(training_file); 
		MultiPhenotype predicting=new MultiPhenotype(predicting_clinic);
		VariantsDouble obs_geno=new VariantsDouble(training_geno);
		VariantsDouble pred_geno=new VariantsDouble(predicting_geno);
		Phenotype obs_pheno=training.phenotypes[0];
		KinshipMatrix obs_matrix=new KinshipMatrix(training_kinship);
		RelationMatrix relation_obs_pred=new RelationMatrix(predicting_relationship);
		double h2=0.3;// 		
		String local_kinship_result_file="/Users/quanlong/Documents/projects2/predictions/RA_data/result/all2/Local_VO.0.Response.deltaDAS.w100000.csv";
		//double[] var_comp=new double[22];
		//int[][] regions22=generate_regions(local_kinship_result_file, var_comp);
		for(int r=12;r<13;r++){
			System.out.println("===Regions: "+r+"===");
			int[][] the_regions=regions[r];
			String[] reginal_kinship_files=new String[the_regions.length];
			String[] relation_matrix_files=new String[the_regions.length];
			for(int i=0;i<the_regions.length;i++){
				reginal_kinship_files[i]=local_K_folder+"final"+".train."+(the_regions[i][0]+1)+"."+
						the_regions[i][1]+"."+the_regions[i][2]+".K.RRM";
				relation_matrix_files[i]=local_K_folder+"final"+".train.test."+(the_regions[i][0]+1)+"."+
						the_regions[i][1]+"."+the_regions[i][2]+".K.RRM";
			}
			double[] var_comp=new double[the_regions.length];
			for(int i=0;i<var_comp.length;i++)var_comp[i]=h2;
			CausalRegions causal=new CausalRegions(the_regions, var_comp, reginal_kinship_files, relation_matrix_files,  
					obs_matrix, relation_obs_pred, obs_geno, pred_geno);
			GBLUP gblup=new GBLUP(obs_geno, pred_geno, obs_pheno, obs_matrix, relation_obs_pred, h2);
			gblup.estimate_local_reg_(causal, var_comp, obs_geno, 1500);
			
			// storing results
			String[] ids0=gblup.pred_pheno.sample_ids.clone();
			double[] values_noclinic=gblup.pred_pheno.values.clone();
			adjust(values_noclinic, obs_pheno.values);

			Phenotype comparison=new MultiPhenotype("/Users/quanlong/Documents/projects2/predictions/RA_data/test_final/all.eli.test.tsv").phenotypes[0];

			double[] betas=training.adjust_cofactors(tobe_adj, adj);
			obs_pheno=training.phenotypes[training.phenotypes.length-1];
			gblup=new GBLUP(obs_geno, pred_geno, obs_pheno, obs_matrix, relation_obs_pred, h2);
			gblup.estimate_local_reg_(causal, var_comp, obs_geno, 1500);
			// again, check results
			String[] ids1=gblup.pred_pheno.sample_ids.clone();
			double[] values_clinic=gblup.pred_pheno.values.clone();
			adjust(values_clinic, obs_pheno.values);
			//add_back
			for(int k=0;k<values_clinic.length;k++){
				values_clinic[k]=values_clinic[k]+betas[0];
				for(int i=0;i<adj.length;i++){
					Phenotype the_tmp_pheno=predicting.phenotypes[predicting.phe_id2index.get(adj[i])];
					if(the_tmp_pheno.sample_id2index.containsKey(ids1[k])){
						int sample_index=the_tmp_pheno.sample_id2index.get(ids1[k]);
						values_clinic[k]=values_clinic[k]+betas[i+1]*
								the_tmp_pheno.values[sample_index];
					}else{
						values_clinic[k]=values_clinic[k]+betas[i+1]*the_tmp_pheno.sample_mean();
					}				
				}
			}
			try{
				String sep=",";
				BufferedWriter bw1=new BufferedWriter(new FileWriter(output_quanti));
				bw1.write("ID"+sep+"pred_gen_facs"+sep+"pred_clin_gen_facs\n");
				BufferedWriter bw2=new BufferedWriter(new FileWriter(output_binary));
				bw2.write("ID"+sep+"belief_gen"+sep+"belief_clin_gen\n");

				for(int k=0;k<ids0.length;k++){
					if(!ids0[k].equals(ids1[k]))System.out.println("Wrong");
					bw1.write(ids0[k]+sep+values_noclinic[k]+sep+values_clinic[k]+"\n");
					bw2.write(ids0[k]+sep+(nonresponder(values_noclinic[k], predicting.phenotypes[0].values[k])?1:0)
							+sep+(nonresponder(values_clinic[k], predicting.phenotypes[0].values[k])?1:0)+"\n");
				}
				bw1.close();bw2.close();
			}catch(Exception e){e.printStackTrace();}
			MultiPhenotype reuslt=new MultiPhenotype(output_quanti);
			double[] corr=MultiPhenotype.correlation(comparison, reuslt.phenotypes[0]);

			System.out.println("non:  "+corr[0]+"\t"+corr[1]+"\t"+corr[2]+"\t"+corr[3]);
			corr=MultiPhenotype.correlation(comparison, reuslt.phenotypes[1]);
			System.out.println("gen:  "+corr[0]+"\t"+corr[1]+"\t"+corr[2]+"\t"+corr[3]);

		}
		

	}
	
	public static void adjust(double[] data, double[] template){
		double mean0=myMathLib.StatFuncs.mean(data);
		double var0=Math.sqrt(myMathLib.StatFuncs.var(data, mean0));		
		double mean1=myMathLib.StatFuncs.mean(template);
		double var1=Math.sqrt(myMathLib.StatFuncs.var(template,mean1));
		for(int k=0;k<data.length;k++){
			data[k]=((data[k]-mean0)/var0)*var1+mean1;
		}
	}
	
	public static boolean nonresponder(double delta_das, double das){
		return (delta_das<=0.6) || ((das-delta_das)>=5.1 && delta_das>0.6 && delta_das<=1.2);
	}
	public static void main(String[] args) {
		
		//nfbc();
		//seaside();
		//ra_drug();
		ra_severity(22);
		ra_severity(19);
		ra_severity(18);
		String folder="/Users/quanlong/Documents/projects2/predictions/RA_data/";
		String predicting_clinic0=folder+"test_final/all.test.tsv";
		String predicting_geno0=folder+"genotypes/all.test.hdf5";
		String predicting_relationship0=folder+"genotypes/all.train.test.K.RRM";
		String output_binary0=folder+"test_final/blup.b.csv";
		String output_quanti0=folder+"test_final/blup.q.csv";
		
		String output_binary1=folder+"test_final/blup.b1.csv";
		String output_quanti1=folder+"test_final/blup.q1.csv";
		
		String predicting_clinic2=folder+"test_final/final.test.tsv";
		String predicting_geno2=folder+"test_final/geno.final.hdf5";
		String predicting_relationship2=folder+"test_final/all.final.K.4.RRM";
		String output_binary2=folder+"test_final/blup.b2.csv";
		String output_quanti2=folder+"test_final/blup.q2.csv";
		
		
		//RA_final_global(predicting_clinic2, predicting_geno2, predicting_relationship2, output_quanti2, output_binary2);
		//RA_final_local(predicting_clinic2, predicting_geno2, predicting_relationship2, output_quanti1, output_binary1);

	}

}
