package mixedmodel;

import java.util.HashMap;

public class CompoudX {

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		
		String input_geno=args[0];
		String input_pheno=args[1];
		String global_kinship_file=args[2];
		String output_folder=args[3];
		String input_compound=args[4];
		String emmax_res=args[5];
		String mString=args[6];
		String the_phe_indexString = args[7];
		int m=Integer.parseInt(mString);
		int the_phe_index=-1;
		if(the_phe_indexString!=null) {
			the_phe_index=Integer.parseInt(the_phe_indexString);
		}
		int min_sample_size=40;
		boolean plot=false;
		
		CompoundXAnalyzer local_k=new CompoundXAnalyzer(input_geno, input_compound);
		MultiPhenotype phenotypeS=new MultiPhenotype(input_pheno);
		if(the_phe_index==-1){
			System.out.println("Running all phenotypes? It is suggested to specify a phenotype index. " +
					"Otherwise it may be slow." +
					"\nLet us have a try!");
			for(int phe_index=0;phe_index<phenotypeS.num_of_pheno;phe_index++){
				String out_phe_file=output_folder+"Compound_VO."+phe_index+"."+phenotypeS.phenotypes[phe_index].phe_id+".r"+".csv";
				local_k.compound_VO(phenotypeS.phenotypes[phe_index], input_geno, global_kinship_file, out_phe_file, plot, min_sample_size,emmax_res, m);
			}
			
		}else{
			String out_phe_file=output_folder+"Compound_VO."+the_phe_index+"."+phenotypeS.phenotypes[the_phe_index].phe_id+".csv";						
			local_k.compound_VO(phenotypeS.phenotypes[the_phe_index], input_geno, global_kinship_file, out_phe_file, plot, min_sample_size, emmax_res, m);								
		
		}

		
	}
}
