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

public class Prediction_RA_Final {

	public static void RA_final_global(String predicting_clinic, String predicting_geno, String predicting_relationship,
			String output_quanti, String output_binary){
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
		double h2=0.3;// from mixed model estimation  		

		GBLUP gblup=new GBLUP(obs_geno, pred_geno, obs_pheno, obs_matrix, relation_obs_pred, h2);
		gblup.estimate_global();
		// storing results
		String[] ids0=gblup.pred_pheno.sample_ids.clone();
		double[] values_noclinic=gblup.pred_pheno.values.clone();
		adjust(values_noclinic, obs_pheno.values);

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
		
		String folder="/Users/quanlong/Documents/projects2/predictions/RA_data/";
//		String predicting_clinic0=folder+"test_final/all.test.tsv";
//		String predicting_geno0=folder+"genotypes/all.test.hdf5";
//		String predicting_relationship0=folder+"genotypes/all.train.test.K.RRM";
//		String output_binary0=folder+"test_final/blup.b.csv";
//		String output_quanti0=folder+"test_final/blup.q.csv";
		
		String output_binary1=folder+"test_final/blup.b1.csv";
		String output_quanti1=folder+"test_final/blup.q1.csv";
		
		String predicting_clinic2=folder+"test_final/final.test.tsv";
		String predicting_geno2=folder+"test_final/geno.final.hdf5";
		String predicting_relationship2=folder+"test_final/all.final.K.4.RRM";
		String output_binary2=folder+"test_final/blup.b2.csv";
		String output_quanti2=folder+"test_final/blup.q2.csv";		
		
		RA_final_global(predicting_clinic2, predicting_geno2, predicting_relationship2, output_quanti2, output_binary2);
		//RA_final_local(predicting_clinic2, predicting_geno2, predicting_relationship2, output_quanti1, output_binary1);

	}

}
