package real_data_analysis;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.HashMap;

import mixedmodel.VariantsDouble;

public class Real {

	public static void combine_phenotype(String file_9, String file_height, String combined){
		HashMap<String, String> height=new HashMap<String, String>();
		try{
			BufferedReader br=new BufferedReader(new FileReader(file_height));
			String line=br.readLine();
			while(line!=null){
				String[] tmp=line.split(" ");
				if(tmp[2].equals("-9")){
					height.put(tmp[1], "NA");
				}else{
					height.put(tmp[1], tmp[2]);
				}
				line=br.readLine();
			}
			br=new BufferedReader(new FileReader(file_9));
			line=br.readLine();
			BufferedWriter bw=new BufferedWriter(new FileWriter(combined));
			while(line!=null){
				String[] tmp=line.split(" ");				
				if(tmp[4].equals("-9")){
					bw.write(tmp[1]+"\t"+"NA\t");
				}else{
					bw.write(tmp[1]+"\t"+tmp[4]+"\t");
				}
				for(int k=2;k<tmp.length;k++){
					if(k!=4){
						if(tmp[k].equals("-9")){
							bw.write("NA\t");
						}else{
							bw.write(tmp[k]+"\t");
						}					
					}					
				}bw.write(height.get(tmp[1])+"\n");
				line=br.readLine();
			}
			bw.close();
		}catch(Exception e){e.printStackTrace();}
	}
	public static void main(String[] args) {
		String data_folder_NFBC="/Volumes/Projects/DATA/Human_GWAS/NFBC/";
		String data_folder_WTCCC="/Volumes/Projects/DATA/Human_GWAS/WTCCC_data/";
		
		String tped_nfbc=data_folder_NFBC+"NFBC.tped";
		String tfam_nfbc=data_folder_NFBC+"NFBC.tfam";
		String csv_nfbc=data_folder_NFBC+"NFBC.char.csv";
		String csv_num=data_folder_NFBC+"NFBC.num.csv";
		
//		VariantsDouble.tped2csv_num(tfam_nfbc, tped_nfbc, csv_nfbc, csv_num);
		
		String file_9=data_folder_NFBC+"phenotype/MetaboPheno.txt";
		String file_height=data_folder_NFBC+"phenotype/pheno.Height";
		String combined=data_folder_NFBC+"NFBC.phenot.tsv";
		combine_phenotype(file_9, file_height, combined);
//		VariantsDouble.char2num(csv_nfbc, csv_num);
	}

}
