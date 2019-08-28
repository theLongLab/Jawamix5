package simulations;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.Arrays;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression;

import mixedmodel.EMMA;
import mixedmodel.EMMAX;
import mixedmodel.KinshipMatrix;
import mixedmodel.MultiPhenotype;
import mixedmodel.Phenotype;
import mixedmodel.VariantsDouble;
import myMathLib.Test;

public class PowerEWAS {	
	
	/*
	 * This method serves for two functions:
	 * 
	 * (1) Simulate phenotype using (multiple) SNPs
	 * (2) Simulate phenotype using (multiple) methylation sites (or regions)
	 * 
	 * As long as the input format is the same, i.e., our HDF5 file converted from the CSV-file, whether it contains
	 *  SNPs or methylation sites doesn't matter.  
	 * 
	 * The String causal_variants is in the form of ChrIndex_LocIndex;ChrIndex_LocIndex;...;ChrIndex_LocIndex
	 * Please note that all the indexes start from ZERO.
	 * 
	 * It supports the following mechanisms:
	 * (1) Additive  
	 * (2) Genetic Heterogeneity  
	 * (3) Effect size: for case/control, it is penetrance; for quantitative data, it is   
	 * 
	 */
	public static void sim_phenotype_mSNPs(String input_genotype, String causal_variants, String output_pheno_file,
			double heritability){
		int[][] causal_indexes=parse_causal_var_string(causal_variants);
		int num_of_causal=causal_indexes.length;
		// assign effect sizes
		double[] effects=new double[num_of_causal];
		NormalDistribution norm=new NormalDistribution();
		for(int i=0;i<effects.length;i++){
			effects[i]=norm.sample();
		}			
		// assign phenotype based on genotype and effect size
		VariantsDouble genotype=new VariantsDouble(input_genotype);
		double[][] phenotypes=new double[1][genotype.sample_size]; // only one trait will be simulated
		for(int var=0;var<num_of_causal;var++){
			double[] genotypes=genotype.load_one_variant_by_index(causal_indexes[var][0], causal_indexes[var][1]);
			for(int indi=0;indi<genotype.sample_size;indi++){					
				phenotypes[0][indi]+=(effects[var]*genotypes[indi]);											
			}
		}				
		// normalize the phenotype to ensure the specified heritability
		normalization(phenotypes[0], heritability);
		// write the phenotype to a file
		write_simulated_pheno_file(phenotypes, output_pheno_file,  genotype);
		System.out.println("Finsihed simulating phenotype with additive model and writing to "+output_pheno_file);
	}
	
	/*
	 * Test the association of (simulated) causal variants and assess their p-value by permutations
	 * EMMAX algorithm is used to control population structure 
	 * Note that it applies to both DNA=> Phenotype and Methylation=>Phenotype. 
	 */
	public static void test_causal_emmax(String genotype_hdf5_file, String kinship_file, String pheno_file, 
			String causal_variants, int permutation, String output_file){
		try{
			int[][] causal_indexes=parse_causal_var_string(causal_variants);
			int num_of_causal=causal_indexes.length;
			Phenotype phenotype=(new MultiPhenotype(pheno_file)).phenotypes[0];
			VariantsDouble genotype=new VariantsDouble(genotype_hdf5_file);
			EMMA emma=new EMMA(phenotype, genotype,  new KinshipMatrix(kinship_file).kinship_matrix.getData());	
			int[] indexes={}, chrs={};
			emma.REMLE(emma.phenotype.values, emma.cbind(chrs, indexes), null);		
			double[][] decomposed_array = emma.reml_decompositioned_array();				
			emma.phenotype.generate_new_Y_by_multiplying(decomposed_array); 
			double[] intsept=new double[emma.sample_size];
			for(int i=0;i<emma.sample_size;i++){
				for(int j=0;j<emma.sample_size;j++){
					intsept[i]+=decomposed_array[i][j];
				}
			}	
			BufferedWriter bw=new BufferedWriter(new FileWriter(output_file));
			bw.write("ChrIndex\tLocIndex\tPermutedPValue\n");
			for(int var=0;var<num_of_causal;var++){
				double[] X_ori=genotype.load_one_variant_by_index(causal_indexes[var][0], causal_indexes[var][1]);
				double[][] Xs_after=new double[emma.sample_size][2];						
				for(int i=0;i<emma.sample_size;i++){
					Xs_after[i][0]=intsept[i];
					for(int j=0;j<emma.sample_size;j++){
						Xs_after[i][1]+=(decomposed_array[i][j]*X_ori[j]);
					}
				}
				// test the variant.
				double[] Y=emma.phenotype.new_Y_4emmax.clone();		
				double real_P=regression_P(Y, Xs_after, true);
				// use permutation to assess p-value
				double[] permuted_Ps=new double[permutation];
				for(int i_permute=0;i_permute<permutation;i_permute++){
					double[] Yp=emma.phenotype.new_Y_4emmax.clone();
					Test.randomPermute(Yp);
					permuted_Ps[i_permute]=regression_P(Yp, Xs_after, true);
				}
				Arrays.sort(permuted_Ps);
				double rank=0.0;
				for(int i_permute=0;i_permute<permutation;i_permute++){
					if(permuted_Ps[i_permute]<real_P)rank=i_permute+1;
					else break;
				}bw.write(causal_indexes[var][0]+"\t"+causal_indexes[var][1]+"\t"+(rank/permutation)+"\n");
			}bw.close();	
		}catch(Exception e){e.printStackTrace();}
	}
	
	/*
	 * Test the association of (simulated) causal variants and assess their p-value by permutations
	 * Linear Model is used to control population structure 
	 * Note that it applies to both DNA=> Phenotype and Methylation=>Phenotype. 
	 */
	public static void test_causal_lm(String genotype_hdf5_file, String pheno_file, 
			String causal_variants, int permutation, String output_file){
		try{
			int[][] causal_indexes=parse_causal_var_string(causal_variants);
			int num_of_causal=causal_indexes.length;
			Phenotype phenotype=(new MultiPhenotype(pheno_file)).phenotypes[0];
			VariantsDouble genotype=new VariantsDouble(genotype_hdf5_file);
			BufferedWriter bw=new BufferedWriter(new FileWriter(output_file));
			bw.write("ChrIndex\tLocIndex\tPermutedPValue\n");
			for(int var=0;var<num_of_causal;var++){
				double[] X_ori=genotype.load_one_variant_by_index(causal_indexes[var][0], causal_indexes[var][1]);
				double[][] Xs_after=new double[X_ori.length][2];						
				for(int i=0;i<X_ori.length;i++){
					Xs_after[i][0]=1.0;
					Xs_after[i][1]=(X_ori[i]);					
				}
				// test the variant.
				double[] Y=phenotype.values.clone();		
				double real_P=regression_P(Y, Xs_after, true);
				// use permutation to assess p-value
				double[] permuted_Ps=new double[permutation];
				for(int i_permute=0;i_permute<permutation;i_permute++){
					double[] Yp=phenotype.values.clone();
					Test.randomPermute(Yp);
					permuted_Ps[i_permute]=regression_P(Yp, Xs_after, true);
				}
				Arrays.sort(permuted_Ps);
				double rank=0.0;
				for(int i_permute=0;i_permute<permutation;i_permute++){
					if(permuted_Ps[i_permute]<real_P)rank=i_permute+1;
					else break;
				}bw.write(causal_indexes[var][0]+"\t"+causal_indexes[var][1]+"\t"+(rank/permutation)+"\n");
			}bw.close();	
		}catch(Exception e){e.printStackTrace();}
	}
	
	public static double regression_P(double[] Y, double[][] Xs, boolean no_intercept){
		OLSMultipleLinearRegression reg1=new OLSMultipleLinearRegression();	
		if(no_intercept)reg1.setNoIntercept(true);
		reg1.newSampleData(Y, Xs);	
		double pvalue=myMathLib.StatFuncs.multi_reg_pvalues(reg1.estimateRegressionParameters(), 
				reg1.estimateRegressionParametersStandardErrors(), Y.length)[1];
		return pvalue; 
	}

	/*
	 * Estimate the overall R2 of multiple mSNPs that explains the variance of a methylated site (Could be a CpG site or a DMR)
	 * Only LM is provided: this is because that the function that simulates phenotype based on a given heritability only uses LM (not EMMAX).
	 */
	public static void mSNP_aggregated_R2(String genotype_hdf5_file, String pheno_file, 
			String causal_variants, String output_file){
		try{
			int[][] causal_indexes=parse_causal_var_string(causal_variants);
			int num_of_causal=causal_indexes.length;
			Phenotype phenotype=(new MultiPhenotype(pheno_file)).phenotypes[0];
			VariantsDouble genotype=new VariantsDouble(genotype_hdf5_file);
			BufferedWriter bw=new BufferedWriter(new FileWriter(output_file));
			bw.write("The Variance component explained by the input mSNPs ("+causal_variants+") is\n");
			double[] Y=phenotype.values.clone();	
			double[][] Xs=new double[Y.length][num_of_causal];		
			for(int var=0;var<num_of_causal;var++){
				double[] X_ori=genotype.load_one_variant_by_index(causal_indexes[var][0], causal_indexes[var][1]);								
				for(int i=0;i<Y.length;i++)	Xs[i][var]=X_ori[i];				
			}
			// test the regression of all the mSNPs.
			OLSMultipleLinearRegression reg1=new OLSMultipleLinearRegression();	
			reg1.newSampleData(Y, Xs);	
			bw.write(reg1.calculateAdjustedRSquared()+"\n");			
			bw.close();	
		}catch(Exception e){e.printStackTrace();}	
	}
	
	public static void write_simulated_pheno_file(double[][] phenos, String file, VariantsDouble genotype){
		try{
			BufferedWriter bw= new BufferedWriter(new FileWriter(file));
			bw.write("ID");
			for(int i=0;i<phenos.length;i++){
				bw.write("\tSim_Pheno_"+i);
			}bw.write("\n");
			for(int indi=0;indi<genotype.sample_size;indi++){
				bw.write(genotype.sample_ids[indi]);
				for(int i=0;i<phenos.length;i++){
					bw.write("\t"+phenos[i][indi]);
				}bw.write("\n");
			}bw.close();
		}catch(Exception e){e.printStackTrace();}
	}
	
	/*
	 * adding the variance of error term
	 * change will be made to the double[] original 
	 */
	public static void normalization(double[] original, double heritability){
		NormalDistribution norm=new NormalDistribution();
		double sigma_g2=StatUtils.populationVariance(original);
		double sigma_e=Math.sqrt(sigma_g2*(1-heritability)/heritability); // heritability=sigma_g2/(sigma_g2+sigma_e2)
		for(int indi=0;indi<original.length;indi++){
			original[indi]+=(sigma_e*norm.sample());
		}
	}
	
	/*
	 * parse the causal_variants string into an index array.
	 */
	public static int[][] parse_causal_var_string(String causal_variants){
		String[] indexes_array=causal_variants.split(";");
		int num_of_causal=indexes_array.length;
		int[][] casual_indexes=new int[num_of_causal][2];
		for(int var=0;var<num_of_causal;var++){
			casual_indexes[var][0]=Integer.parseInt(indexes_array[var].split("_")[0]);  // ChrIndex
			casual_indexes[var][1]=Integer.parseInt(indexes_array[var].split("_")[1]);	// LocIndex
		}
		return casual_indexes;
	}	
	
	public static void main(String[] args) {
		if(args.length==0){
			System.out.println("EWAS/GWAS Power estimator. Four functions: \n"
					+ "\tsim-phenotype\n"
					+ "\tmSNPsR2\n"
					+ "\tpvalue-lm\n"	
					+ "\tpvalue-emmax\n"	);
			System.exit(0);
		}
		String function=args[0];
		if(args[0].equals("pvalue-lm")){
			if(args.length==1){
				System.out.println("Conduct association test using linear model and calcualte the P-value using permutation.\n"
						+ "Usage: <genotype_file> <causal_variants_string> <phenotype_file> <#permutations> <output_file>\n"
						+ "\tNote: \n\tgenotype_file must be an HDF5 file. \n"
						+ "\tcausal_variants_string is in the format of: ChrIndex_LocIndex;ChrIndex_LocIndex;...;ChrIndex_LocIndex. \n"
						+ "\t#permutations can be an integer between 200 and 2000.");
			}else{
				String genotype_hdf5_file=args[1];
				String causal_variants=args[2];
				String pheno_file=args[3];
				int permutation=Integer.parseInt(args[4]);
				String output_file=args[5];
				test_causal_lm(genotype_hdf5_file, pheno_file, causal_variants, permutation, output_file);
			}
		}else if(args[0].equals("pvalue-emmax")){
			if(args.length==1){
				System.out.println("Conduct association test using linear mixed model (EMMAX algorithm) and calcualte the P-value using permutation.\n"
						+ "Usage: <genotype_file> <causal_variants_string> <phenotype_file> <#permutations> <kinship_file> <output_file>\n"
						+ "\tNote: \n\tgenotype_file must be an HDF5 file. \n"
						+ "\tcausal_variants_string is in the format of: ChrIndex_LocIndex;ChrIndex_LocIndex;...;ChrIndex_LocIndex. \n"
						+ "\t#permutations can be an integer between 200 and 2000.\n"
						+ "\tkinship_file can be generated using jawamix5.jar.\n");
			}else{
				String genotype_hdf5_file=args[1];
				String causal_variants=args[2];
				String pheno_file=args[3];
				int permutation=Integer.parseInt(args[4]);
				String kinship_file=args[5];
				String output_file=args[6];
				test_causal_emmax(genotype_hdf5_file, kinship_file, pheno_file, causal_variants, permutation, output_file);
			}
		}else if(args[0].equals("sim-phenotype")){
			if(args.length==1){
				System.out.println("Simulate phenotype based on genotype (of DNAm or DNA).\n"
						+ "Usage: <genotype_file> <causal_variants_string> <output_phenotype_file> <heritability>\n"
						+ "\tNote:\n\t genotype_file must be an HDF5 file. \n"
						+ "\tcausal_variants_string is in the format of: ChrIndex_LocIndex;ChrIndex_LocIndex;...;ChrIndex_LocIndex. \n"
						+ "\theritability is a number between 0 and 1.");
			}else{
				String input_genotype=args[1];
				String causal_variants=args[2];
				String output_pheno_file=args[3];
				double heritability=Double.parseDouble(args[4]);
				sim_phenotype_mSNPs(input_genotype, causal_variants, output_pheno_file, heritability);
			}
		}else if(args[0].equals("mSNPsR2")){
			if(args.length==1){
				System.out.println("Calcualte the total methylation variance explained by multiple mSNPs.\n"
						+ "Usage: <genotype_file> <causal_variants_string> <phenotype_file> <output_file>\n"
						+ "\tNote:\n\t genotype_file must be an HDF5 file. \n"
						+ "\tcausal_variants_string is in the format of: ChrIndex_LocIndex;ChrIndex_LocIndex;...;ChrIndex_LocIndex. \n"
						+ "");
			}else{
				String genotype_hdf5_file=args[1];
				String causal_variants=args[2];
				String pheno_file=args[3];
				String output_file=args[4];
				mSNP_aggregated_R2(genotype_hdf5_file, pheno_file, causal_variants, output_file);
			}
		}else{
			System.out.println("The function "+function+" doesn't exist. Typo?\n");
		}

	}

}
