package mixedmodel;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;

import javax.naming.spi.DirStateFactory.Result;

import myMathLib.Test;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularValueDecomposition;
import org.apache.commons.math3.optimization.GoalType;
import org.apache.commons.math3.optimization.univariate.BrentOptimizer;
import org.apache.commons.math3.optimization.univariate.UnivariatePointValuePair;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

public class FaSTLMM_Compound extends FaSTLMM {
	
	public static int min_sample_size=10;
	RealMatrix global_Ut; // eigen vectors 
	final double[] S; //eigen values
	final double[] Y_K_transformed_null;
	final double[] Y_global_transformed;
	
	final double[] intercept;
	double Likelihood_null; // could be REML or ML, depending on the parameter passed to constructor.
	double sigma_g_null;
	double sigma_e_null;
	double delta_null;
		
	final RealMatrix transformation_matrix; // =U'(S+delta*I)^(-1/2) 
	
	/* inherited:	
	public final static double machep=1e-10;//No need to use 2.220446049250313E-16; // Machine epsilon = Math.ulp(1d)
	public final static double brent_rel=Math.sqrt(machep);
	public final static double brent_abs=brent_rel;
	public final int maxEval=100;	
	public final static double uplimit=10;
	public final static double lowlimit=-10;
	public final static double grid_step=0.1;	
	
	String genotype_hdf5_file;
	VariantsDouble genotype;
	int[] indexes_with_phenotype_in_genotype;
	Phenotype phenotype; // all samples with phenotype will have genotype in the object: the constructor will guarantee that by checking the data.
	double[][] global_kinship_matrix; // this is corresponding to the phenotype, but generated from genotype matched kinship file.
	int sample_size;			
	*/
	
	public FaSTLMM_Compound(Phenotype ori_phenotype, String genotype_file_hdf5, String global_kinship_file){
		super(ori_phenotype, genotype_file_hdf5, global_kinship_file);		
		EigenDecomposition eigen=new EigenDecomposition(new Array2DRowRealMatrix(this.global_kinship_matrix, false),0); 
		this.global_Ut=eigen.getVT();
		this.Y_K_transformed_null=this.global_Ut.operate(this.phenotype.values);
		this.S=eigen.getRealEigenvalues();
		eigen=null;
		this.global_kinship_matrix=null; // free memory
		System.out.println("Eigen-decomposition of global kinship done.");
		this.fullRankSolver_null_ML();
		System.out.println("Null ML with global kinship calculated.");
		// set up null-ML in the model that only fixed effect is the transformed intercept. 	    
	    double[][] transformation_matrix_array= new double[this.sample_size][this.sample_size];
		double[] diag= new double[this.sample_size];
		for(int i=0;i<this.sample_size;i++){
	//		if(this.reml_eig_L.eigenvalues[i]>0)
				diag[i]=1.0/Math.sqrt(this.S[i]+this.delta_null); //(S+delta_null)^(-1/2), no initialization of delta_null, delta_null=0
		}
		for(int i=0;i<this.sample_size;i++){
			for(int j=0;j<this.sample_size;j++){
				transformation_matrix_array[i][j]=this.global_Ut.getEntry(i,j)*diag[i]; //Get the entry in the specified row and column.
			}
		}
		this.global_Ut=null;
		this.transformation_matrix=new Array2DRowRealMatrix(transformation_matrix_array, false);
		System.out.println("Transformation Matrix formed.");
		this.Y_global_transformed=this.transformation_matrix.operate(this.phenotype.values);
		double[] all_one=new double[this.sample_size];
		Arrays.fill(all_one, 1);
		this.intercept=this.transformation_matrix.operate(all_one);		
		// remove rescaled intercept: NO NEED TO DO?		
		double mean_y=(new DescriptiveStatistics(this.Y_global_transformed)).getMean();
		double mean_beta0=(new DescriptiveStatistics(this.intercept)).getMean();
		double ratio=mean_y/mean_beta0;
		for(int k=0;k<this.Y_global_transformed.length;k++){
			this.Y_global_transformed[k]-=(ratio*this.intercept[k]);
		}
	}
	
	public void fullRankSolver_null_ML(){
		double[][] transformed_X=new double[1][];
		double[] intercept_1=new double[this.sample_size];
		Arrays.fill(intercept_1, 1);
		transformed_X[0]=this.global_Ut.operate(intercept_1);
		BrentOptimizer bo=new BrentOptimizer(FaSTLMM.brent_rel, FaSTLMM.brent_abs);
		double best_delta=-1;
		double best_ml=-(Double.MAX_VALUE);
		UnivariateFunction ml_function=new MLFullRank4BrentOptimal(transformed_X,Y_K_transformed_null,S);
		for(double min=lowlimit;min<uplimit;min+=grid_step){
			double min_delta=Math.exp(min);
			double max_delta=Math.exp(min+grid_step);
			UnivariatePointValuePair result=bo.optimize(maxEval, ml_function, GoalType.MAXIMIZE, min_delta, max_delta);
			double the_ml=result.getValue();
			if(the_ml>best_ml){
				best_ml=the_ml;
				best_delta=result.getPoint();
			}
		}
		this.Likelihood_null=best_ml;		
		double[] beta_0=MLFullRank4BrentOptimal.calculate_beta(Y_K_transformed_null, transformed_X, S, best_delta);
		this.sigma_g_null=MLFullRank4BrentOptimal.sigma_g2_LL(Y_K_transformed_null, transformed_X, S, best_delta, beta_0);
		this.sigma_e_null=this.sigma_g_null*best_delta;
		this.delta_null=best_delta;
	}
	
	public void analysis_specified_compounds(String mlr_output_file, int[][][] compounds, boolean plot){
		try{
			System.out.println(phenotype.phe_id);		
			BufferedWriter bw=new BufferedWriter(new FileWriter(mlr_output_file));
			bw.write("#SampleSize="+this.sample_size+ "; REML="+this.Likelihood_null+"; " +
					"h0-heritability="+this.sigma_g_null/(this.sigma_e_null+this.sigma_g_null) +"; delta="+this.delta_null+"\n");
			bw.write("compound(chr:start:end_chr:start_end), best_ml, p-value, variance_local, variance_global, variance_e\n");	
			final double[] all_one=new double[this.sample_size];
			Arrays.fill(all_one, 1);
			final double[] all_one_gt=this.transformation_matrix.operate(all_one); //Returns the result of multiplying this by the vector v.
			int total=0;
			int[] low_counts={0};			
			for(int compound_index=0;compound_index<compounds.length;compound_index++){
				String compound_name="";
				total++;
				//int start=1, end=start+tiling_win_size-1;
				ArrayList<double[][]> data_buffer=new ArrayList<double[][]>();
				int num_snps_in_compound=0;
				for(int region_index=0;region_index<compounds[compound_index].length;region_index++){					
					int chr=compounds[compound_index][region_index][0];
					int start=compounds[compound_index][region_index][1];
					int end=compounds[compound_index][region_index][2];
					double[][] snp_in_this_region=genotype.load_variants_in_region(chr-1, start, end);
					data_buffer.add(snp_in_this_region);
					if(region_index==compounds[compound_index].length-1)
						compound_name=compound_name+chr+":"+start+":"+end;
					else compound_name=compound_name+chr+":"+start+":"+end+"_";
					num_snps_in_compound+=snp_in_this_region.length;
				}
				double[][] combined_data=new double[num_snps_in_compound][];
				int combined_index=0;
				for(int region_index=0;region_index<compounds[compound_index].length;region_index++){	
					double[][] snp_in_this_region=data_buffer.get(region_index);
					for(int k=0;k<snp_in_this_region.length;k++){
						combined_data[combined_index++]=snp_in_this_region[k];
					}
				}
				double[][] data_matching_phenotype=this.select_subset_of_data(combined_data);
				if(data_matching_phenotype.length<2){
					//System.out.println("NO SNP: "+chr+"/"+start);
					continue;
				}
				run_a_compound(data_matching_phenotype, all_one_gt, low_counts, bw, compound_name);		
			}bw.close();
			System.out.println("Total #compounds="+total+"; #Low-rank regions="+low_counts[0]);
			//if(plot)LocalKinshipAnalyzer.make_plot_one_phen(mlr_output_file);
		}catch(Exception e){e.printStackTrace();}
	}
	
	public void run_a_compound(double[][] data_matching_phenotype, double[] all_one_gt, int[] low_full_counts, BufferedWriter bw, 
			String compound_name){
		//System.out.println(chr+"/"+start);
		try{
			FaSTLMM.normalize_terms_withSqrtNumVarDivided(data_matching_phenotype);
			final double[] all_one=new double[this.sample_size];
			Arrays.fill(all_one, 1);
			double[][] transformed_data=FaSTLMM.transform_X_by_Ut(transformation_matrix, data_matching_phenotype);
			double[][] X_for_null=new double[1][];
			X_for_null[0]=all_one_gt.clone();
			
			VarComp_Result result_null=FaSTLMM_FullRank.fullRankSolver_FixedEffects_ML(X_for_null, this.Y_global_transformed.clone(), 
					all_one, 0);
			double ml_null=result_null.ml;
			
			
			if(data_matching_phenotype.length<sample_size){ //low rank
				//System.out.println("Low-Rank chr"+(chr+1)+" win "+start);
				low_full_counts[0]++;
				double[] Y_new_original=this.Y_global_transformed;
				//FaSTLMM.normalize_terms_withSqrtNumVarDivided(transformed_data); NO NEED TO NORMALIZE AGAIN!!!
				SingularValueDecomposition local_svd=new SingularValueDecomposition(
						(new Array2DRowRealMatrix(transformed_data)).transpose());	
				
				double[] S=local_svd.getSingularValues();
				for(int k=0;k<S.length;k++)S[k]=(S[k]*S[k]);
				RealMatrix local_U1t= local_svd.getUT();
				double[] transformed_Y=local_U1t.operate(Y_new_original);
				RealMatrix local_IminusU1U1t=local_svd.getU().multiply(local_U1t).scalarMultiply(-1);
				for(int i=0;i<this.sample_size;i++)local_IminusU1U1t.setEntry(i, i, 1+local_IminusU1U1t.getEntry(i, i));
				double[] Y_transformed2=local_IminusU1U1t.operate(Y_new_original);
				double[][] X_new_original=new double[1][];
				X_new_original[0]=all_one_gt;
				//double[][] transformed_X=new double[1][];
				//transformed_X[0]=local_U1t.operate(all_one);						
				VarComp_Result result=FaSTLMM_LowRank.lowRankSolver_FixedEffects(X_new_original, Y_new_original, 
						transformed_Y,Y_transformed2, S, ml_null, local_U1t, local_IminusU1U1t);								
				if(result!=null){
					double total_variance=result.sigma_g+result.sigma_e;//*(1+this.delta_null);
					bw.write(compound_name + ", "+result.ml+", "+result.pvalue+", "+
							result.sigma_g/total_variance+", "+(result.sigma_e/total_variance)*(1/(1+this.delta_null))+", "+ 
							(result.sigma_e/total_variance)*(this.delta_null/(1+this.delta_null))+"\n");
					//System.out.println("Succ finsihed chr"+(chr+1)+" win "+start);
				}
			}else if(data_matching_phenotype.length>=sample_size){	// full-rank						
				RealMatrix local_kinship= VariantsDouble.calculate_RRM_local_kinship_RealMatrix_noNomalization(transformed_data);
				//RealMatrix local_kinship=KinshipMatrix.re_scale_kinship_matrix(local_kinship0, true);
				double[] Y_new_original=this.Y_global_transformed;
				EigenDecomposition local_eigen= null;
				try{local_eigen= new EigenDecomposition((Array2DRowRealMatrix)local_kinship,0);}
				catch(Exception e){}
				if(local_eigen!=null){
					double[] S=local_eigen.getRealEigenvalues();
					RealMatrix local_kinship_Ut= local_eigen.getVT();
					double[] transformed_Y=local_kinship_Ut.operate(Y_new_original);
					double[][] transformed_X=new double[1][];
					transformed_X[0]=local_kinship_Ut.operate(all_one_gt);
					//double ml_null=myMathLib.StatFuncs.null_logL(transformed_Y,transformed_X[0]);
					//System.out.println(ml_null+","+this.Likelihood_null);
					VarComp_Result result=FaSTLMM_FullRank.fullRankSolver_FixedEffects_ML(transformed_X, transformed_Y, 
							S, ml_null);
					
					if(result!=null){
						double total_variance=result.sigma_g+result.sigma_e;//*(1+this.delta_null);
						bw.write(compound_name+", "+result.ml+", "+result.pvalue+", "+
								result.sigma_g/total_variance+", "+(result.sigma_e/total_variance)*(1/(1+this.delta_null))+", "+ 
								(result.sigma_e/total_variance)*(this.delta_null/(1+this.delta_null))+"\n");
				}				
				//	System.out.println("Succ finsihed chr"+(chr+1)+" win "+start);
				}
			}	
		}catch(Exception e){e.printStackTrace();}
	}
	
	
}
