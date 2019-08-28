package real_data_analysis;

import java.io.BufferedWriter;
import java.io.FileWriter;

public class GenerateBsubs {

	
	/*
	 * <genotype> <g-kinship> <model> <working_folder> <obs_sample_size> <first_sample_size> <sim_win (100000)> <min_MAF> <max_MAF>
	 */
	public static void generated_commands_common(String command_file){
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(command_file));
			String genotype="/sc/orga/work/longq01/rGBLUP/g1k_all.hdf5";
			String kinship="/sc/orga/work/longq01/rGBLUP/g1k_all.K.RRM";
			String[] models={"dorminant","epistasis","additive","liability","canalization","signepistasis","compensation"};
			String working_folder="/sc/orga/scratch/longq01/projects/prediction/rgblup-2015-feb-common/";
			String obs_sample_size="700";
			String first_sample_size ="500";
			String sim_win="100000";
			double[] maf={0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3};
			for(int k=0;k<4;k++){
				for(int i=0;i<maf.length-1;i++){
					bw.write("bsub -q alloc -W 144:00 -P acc_GTEX -R \"rusage[mem=6000]\" java -Xmx6g -jar "
							+ "/hpc/users/longq01/Documents/programs/my_jar/rgblup_sim_common.jar "+genotype+" "+
							kinship+" "+models[k]+" "+working_folder+models[k]+"/ "+obs_sample_size+" "+first_sample_size+" "+sim_win+" "+
							maf[i]+" "+maf[i+1]+"\n");
				}
			}for(int k=4;k<=6;k++){				
				bw.write("bsub -q alloc -W 144:00 -P acc_GTEX -R \"rusage[mem=6000]\" java -Xmx6g -jar "
							+ "/hpc/users/longq01/Documents/programs/my_jar/rgblup_sim_common.jar "+genotype+" "+
						kinship+" "+models[k]+" "+working_folder+models[k]+"/ "+obs_sample_size+" "+first_sample_size+" "+sim_win+" "+
							"0.2 0.5\n");				
			}
			bw.close();
		}catch(Exception e){e.printStackTrace();}
	}
	
	/*
	 * /Volumes/Projects/DATA/1000G/1000gdata/g1k_all.hdf5 /Volumes/Projects/DATA/1000G/1000gdata/g1k_all.K.RRM additive /Users/quanlong/Documents/projects2/predictions/simulation/1000g/2014-11-18-rare/additive/  700 200 100000 0.002 0.005
	 */
	public static void generated_commands_rare(String command_file){
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(command_file));
			String genotype="/sc/orga/work/longq01/rGBLUP/g1k_all.hdf5";
			String kinship="/sc/orga/work/longq01/rGBLUP/g1k_all.K.RRM";
			String[] models={"dorminant","epistasis","additive","liability","canalization","signepistasis","compensation"};
			String working_folder="/sc/orga/scratch/longq01/projects/prediction/rgblup-2015-feb-rare/";
			String obs_sample_size="700";
			String first_sample_size ="500";
			String sim_win="100000";
			double[] maf={0, 0.002, 0.005, 0.01};
			double[] concentration_rare={0.03, 0.05, 0.08, 0.12, 0.2, 0.4, 0.6, 0.8};
			for(int k=0;k<4;k++){
				for(int i=0;i<maf.length-1;i++){
					for(int j=0;j<concentration_rare.length;j++){
						bw.write("bsub -q alloc -W 144:00 -P acc_GTEX -R \"rusage[mem=6000]\" java -Xmx6g -jar "
								+ "/hpc/users/longq01/Documents/programs/my_jar/rgblup_sim_rare.jar "+genotype+" "+
								kinship+" "+models[k]+" "+working_folder+models[k]+"/ "+obs_sample_size+" "+first_sample_size+" "+sim_win+" "+
								maf[i]+" "+maf[i+1]+" "+concentration_rare[j]+"\n");
					}
					
				}
			}
//			for(int k=4;k<=6;k++){				
//				bw.write("bsub -q scavenger -W 23:59 -R \"rusage[mem=6000]\" java -Xmx6g -jar "
//							+ "/hpc/users/longq01/Documents/programs/my_jar/rgblup_sim.jar "+genotype+" "+
//						kinship+" "+models[k]+" "+working_folder+models[k]+"/ "+obs_sample_size+" "+first_sample_size+" "+sim_win+" "+
//							"0.2 0.5\n");				
//			}
			bw.close();
		}catch(Exception e){e.printStackTrace();}
	}
	
	public static void ipm(String output_file){
		String command="bsub -q alloc -W 144:00 -P acc_GTEX -R \"rusage[mem=12000]\" java -Xmx12g -jar "
				+ "/hpc/users/longq01/Documents/programs/my_jar/jawamix5.jar ";
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(output_file));
			int icd_num=4247, reg_num=506;
			for(int k=0;k<reg_num;k++){
				bw.write(command + "emmax -ig ipm.hdf5 -ip pheno_regular.tsv -o emmax_reg/ -ik ipm.RRM -min_size 1000 -index "+k+"\n");
				bw.write(command + "local -ig ipm.hdf5 -ip pheno_regular.tsv -o local_reg/ -w 100000 -ik_g ipm.RRM -min_size 1000 -index "+k+"\n");
			}
			for(int k=0;k<icd_num;k++){
				bw.write(command + "emmax -ig ipm.hdf5 -ip pheno_icd.tsv -o emmax_icd/ -ik ipm.RRM -min_size 1000 -index "+k+"\n");
				bw.write(command + "local -ig ipm.hdf5 -ip pheno_icd.tsv -o local_icd/ -w 100000 -ik_g ipm.RRM -min_size 1000 -index "+k+"\n");
			}bw.close();
		}catch(Exception e){e.printStackTrace();}
	}
	public static void main(String[] args) {
		 String command_file="/Users/quanlong/Documents/projects2/predictions/simulation/1000g/2014-11-17/command2014-11-17.txt";
		 String rare_command_file="/Users/quanlong/Documents/projects2/predictions/simulation/1000g/2014-11-17/command2015-02.rare.txt";
		 String common_command_file="/Users/quanlong/Documents/projects2/predictions/simulation/1000g/2014-11-17/command2015-02.common.txt";
//		 generated_commands_common(common_command_file);
//		 generated_commands_rare(rare_command_file);
		 
		 ipm("/Users/quanlong/Documents/projects2/predictions/IPM/commands.txt");
	}

}
