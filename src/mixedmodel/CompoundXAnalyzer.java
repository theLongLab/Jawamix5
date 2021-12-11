package mixedmodel;

public class CompoundXAnalyzer extends CompoundAnalyzer{

	public CompoundXAnalyzer(String hdf5_file, String compound_coordinates_file) {
		super(hdf5_file, compound_coordinates_file);
		// TODO Auto-generated constructor stub
	}
	
	public void compound_VO(Phenotype phenotype, String genotype_hdf5_file, String global_kinship_file, 
			String mlr_output_file, boolean plot, int min_sample_size, String emmax_res, int m){
		int the_sample_size=FaSTLMM.sample_size(phenotype, genotype_hdf5_file);
		if(the_sample_size< min_sample_size){
			System.out.println("Overlap of sample size too small: "+the_sample_size);
			return;
		}
		FaSTLMM_CompoundX calculator=new FaSTLMM_CompoundX(phenotype, genotype_hdf5_file, global_kinship_file);
		calculator.analysis_specified_compounds(mlr_output_file, super.compound_coordinates, plot, emmax_res,m);
	}

}
