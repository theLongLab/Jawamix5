package real_data_analysis;

import mixedmodel.VariantsDouble;
import simulations.AdhocFunctions4Analysis;
import simulations.Simulator;

public class TypeIError {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// analyze random data to convert the type-I-error
		
//		String input_geno="/Users/quan.long/Documents/Projects/GeneMappingMethodology/LocalGlobal/simulations/" +
//				"swe180_ecker171_removebads.csv.double.hdf5";
		String input_geno="/Volumes/Projects/DATA/arabidopsis_genomes/swe180_ecker171_removebads.csv.num.mafc.hdf5";
		VariantsDouble data= new VariantsDouble(input_geno);
		String input_pheno="/Volumes/Projects/Jawamix5/simulations/tmp.phe.txt";
		Simulator.generate_random_phenotype(input_pheno, data.sample_ids);
		String output_folder="/Volumes/Projects/Jawamix5/simulations/type1error/";
//		String local_kinship_files_folder="/Users/quan.long/Documents/Projects/GeneMappingMethodology/LocalGlobal/" +
//				"simulations/kinship_local_50k/";
		String local_kinship_files_folder="/Volumes/Projects/Jawamix5/simulations/kinship_local_100k/";
		String global_kinship_file="/Volumes/Projects/DATA/arabidopsis_genomes/kinship.rescaled.ibs";
		int win_size=100000;
		AdhocFunctions4Analysis.analyze_local(input_geno, input_pheno, output_folder, local_kinship_files_folder, 
				global_kinship_file, win_size);
		
	}
}
