package real_data_analysis;

import mixedmodel.CompoundAnalyzer;
import mixedmodel.EMMAX;
import mixedmodel.FaSTLMM_LocalKinship;
import mixedmodel.LocalKinshipAnalyzer;
import mixedmodel.MultiPhenotype;
import mixedmodel.Phenotype;
import mixedmodel.RareAnalyzerAggregate;
import mixedmodel.Regions;

public class Debug_in_details {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
////		String pheotype="/Volumes/Projects/DATA/arabidopsis_phenotypes/big_table_gmi/sinle_pheno_for_debug/176_Leaf_roll_10.tsv";
////		String input_genotype="/Volumes/Projects/DATA/arabidopsis_genomes/swe180_ecker171_removebads.csv.num2.hdf5";		
//		//String global_kinship_file="/Volumes/Projects/DATA/arabidopsis_genomes/kinship_files_swe180_ecker171_removebads/kinship.swe180_ecker171_removebads.rescaled.ibs";
////		String global_kinship_file="/Volumes/Projects/DATA/arabidopsis_genomes/swe180_ecker171_removebads.csv.num2.hdf5.kinship.rescaled.ibs";
////		int win_size=100000;
////		String out_folder="/Volumes/Projects/Local-kinship/simulations/have2look/"; 
//		
////		MultiPhenotype phenotypeS=new MultiPhenotype(pheotype);
////		for(int phe_index=0;phe_index<phenotypeS.num_of_pheno;phe_index++){
////			String out_phe_file=out_folder+"Local_VO."+phe_index+"."+phenotypeS.phenotypes[phe_index].phe_id+".w"+
////					win_size+".again.csv";
////			FaSTLMM_LocalKinship calculator=new FaSTLMM_LocalKinship(phenotypeS.phenotypes[phe_index], input_genotype, 
////					global_kinship_file);
////			calculator.analysis_tiling_win(out_phe_file, win_size, true);			
////		}
//		
////		String folder="/Volumes/Projects/Local-kinship/simulations/sim_rare_common/SIM.h0.9.min0.0.max0.01.con0.4.w10000.chr3.start1585369.end1608271.round0/";
////		String files[]={"Aggregate.Sim_Pheno_0.csv","EMMAX.0_Sim_Pheno_0.top","Local_VO.Sim_Pheno_0.gene.csv","Local_VO.Sim_Pheno_0.w50000.csv"};
////		double[] cut_off={EntranceSim.aggregate_arabidopsis_threshold,EntranceSim.emmax_arabidopsis_threshold,
////				EntranceSim.local_gene_arabidopsis_threshold,EntranceSim.local_win50k_arabidopsis_threshold};
////		for(int i=0;i<4;i++){
////			System.out.println(EntranceSim.check_power(folder+files[i], 2, 1585369-50000, 1608271, cut_off[i])+"\t"+files[i]);
////		}
//		
//		String results_output_folder="/Volumes/Projects/Local-kinship/simulations/sim_rare_common_debug/";
//		String simulated_phenotype_file=results_output_folder+"SIM.h0.9.min0.0.max0.01.con0.5.w10000.chr2.start12114415.end12140069.round0/simulated_phenotype.tsv";
//		
//		String regions_file="/Volumes/Projects/DATA/Annotations/arabidopsis/TAIR10_GFF3_genes.gene2Konly.tsv";
//		String genotype_file="/Volumes/Projects/DATA/arabidopsis_genomes/swe180_ecker171_removebads.csv.num2.hdf5";
////		String global_kinship_file="/Volumes/Projects/DATA/arabidopsis_genomes/swe180_ecker171_removebads.csv.num2.hdf5.kinship.rescaled.ibs";
//		
//		String global_kinship_file="/Volumes/Projects/DATA/arabidopsis_genomes/kinship_files_swe180_ecker171_removebads/swe180_ecker171_removebads.2.rescaled.IBS";
//
//		int chr=1, start=12114415+10000, end=12140069-10000;
//		
//		Phenotype phenotype=new MultiPhenotype(simulated_phenotype_file).phenotypes[0];
//		FaSTLMM_LocalKinship calculator=new FaSTLMM_LocalKinship(phenotype, genotype_file, 
//				global_kinship_file);
//		int[][][] regions=(new Regions(regions_file, calculator.genotype)).region_coordinates;
//		int[][][] region_only1=new int[regions.length][][];
//		String[][] region_name=new String[regions.length][];
//		for(int chr_index=0;chr_index<regions.length;chr_index++){
//			if(chr_index==chr){
//				region_only1[chr]=new int[1][2];
//				region_only1[chr][0][0]=start;
//				region_only1[chr][0][1]=end;
//				region_name[chr]=new String[1];
//				region_name[chr][0]="The_Gene";
//			}else{
//				region_only1[chr_index]=new int[0][];
//				region_name[chr_index]=new String[0];
//			}
//		}
//		//aggregate
//		String aggregate_outfile=results_output_folder+"/Aggregate."+phenotype.phe_id+".csv";
//		RareAnalyzerAggregate rare=new RareAnalyzerAggregate(genotype_file, region_only1, region_name, phenotype,
//				global_kinship_file, 0.01);
//		rare.rare_association(0.01, results_output_folder, false);	
//		// emmax, region
//		String emmax_outputfile =results_output_folder+"/EMMAX.0_"+phenotype.phe_id+".top";
//		EMMAX.emmax_analysis_regions(genotype_file, simulated_phenotype_file, global_kinship_file, results_output_folder,  
//				1000000, 10, 0, 0.05, false, region_only1);
//		//local
////		String local_win_outfile=results_output_folder+"/Local_VO."+phenotype.phe_id+".w"+
////				analysis_win_size+".csv";
//		String local_gene_outfile=results_output_folder+"/Local_VO."+phenotype.phe_id+".gene"+".csv";								
//		calculator.analysis_specified_region(local_gene_outfile, region_only1, false);
//		System.out.println(Math.log10(7.361884745083473E-15)+"\t"+Math.log10(2.220446049250313E-16));
//		System.out.println("EMMAX"+start+"."+end+"\t"+EntranceSim.check_power(emmax_outputfile, chr, start, end, EntranceSim.emmax_arabidopsis_threshold)+"\t"+(-Math.log10(EntranceSim.emmax_arabidopsis_threshold))+"\n");
//		System.out.println("Local_gene"+start+"."+end+"\t"+EntranceSim.check_power(local_gene_outfile, chr, start, end, EntranceSim.local_gene_arabidopsis_threshold)+"\t"+(-Math.log10(EntranceSim.local_gene_arabidopsis_threshold))+"\n");
//		System.out.println("Aggregate"+start+"."+end+"\t"+EntranceSim.check_power(aggregate_outfile, chr, start, end, EntranceSim.aggregate_arabidopsis_threshold)+"\t"+(-Math.log10(EntranceSim.aggregate_arabidopsis_threshold))+"\n");
		
		String data_folder="/Users/quanlong/Documents/projects/prediction/3D/";
		String compound_coordinates_file=data_folder+"HiC/tmp.txt";
		String hdf5_1000g=data_folder+"1000g/g1k_all.hdf5";
		String kinship_1000g=data_folder+"1000g/g1k_all.K.05.RRM";
		String hdf5_gtex=data_folder+"gtex/genotype/GTEx_imputed_dosage.hdf5";
		
		CompoundAnalyzer ca=new CompoundAnalyzer(hdf5_1000g, compound_coordinates_file);
		String simulated_phenotype_file=data_folder+"sim/pheno.tsv";
		int[] compound_index={0, 1, 2, 4, 15};
		double[] effect_size={0.9,0.7,0.5,0.3,0.1};
		Sim3D.sim_pheno(ca, effect_size, compound_index, simulated_phenotype_file);
		Phenotype phenotype=new MultiPhenotype(simulated_phenotype_file).phenotypes[0];
		FaSTLMM_LocalKinship calculator=new FaSTLMM_LocalKinship(phenotype, hdf5_1000g, 
				global_kinship_file);
		
	}

}
