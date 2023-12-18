package simulations;

import myMathLib.StatFuncs;
import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;
import java.util.HashMap;
import org.apache.commons.math3.distribution.UniformRealDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularValueDecomposition;
import java.io.BufferedWriter;
import java.io.FileWriter;

import mixedmodel.KinshipMatrix;
import mixedmodel.MultiPhenotype;
import mixedmodel.Phenotype;
import mixedmodel.VariantsDouble;

public class SimedLMM
{
    private static String edLMMIntro;
    private static String supported_functions;
    
    static {
        SimedLMM.edLMMIntro = "Simulation programme for edLMM\nDeveloper: Qing Li and Quan Long\nUsage: java -Xmx2g -jar SimedLMM.jar [function]";
        SimedLMM.supported_functions = "\n\tnull_h\n\talter_h";
    }
    
    public static boolean check_control_case_balance(final double[] sim_pheno_one_round) {
        boolean control_case_balance = false;
        int control_count = 0;
        int case_count = 0;
        for (int i = 0; i < sim_pheno_one_round.length; ++i) {
            if (Double.compare(sim_pheno_one_round[i], 1.0) == 0) {
                ++case_count;
            }
            else if (Double.compare(sim_pheno_one_round[i], 2.0) == 0) {
                ++control_count;
            }
            else {
                System.out.println("Error in simulated binary phenotype");
            }
        }
        final double control_case_ratio_lower_bound = 0.9;
        final double control_case_ratio_upper_bound = 1.1;
        final double control_case_ratio = control_count / (double)case_count;
        if (Double.compare(control_case_ratio, control_case_ratio_lower_bound) >= 0 && Double.compare(control_case_ratio, control_case_ratio_upper_bound) <= 0) {
            control_case_balance = true;
        }
        return control_case_balance;
    }
    
    public static int mafc(final double[] one_variant_geno) {
        int zero_count = 0;
        int other_count = 0;
        for (int k = 0; k < one_variant_geno.length; ++k) {
            final int the_allele = (int)(one_variant_geno[k] + 0.1);
            if (the_allele == 0) {
                zero_count += 2;
            }
            else if (the_allele == 1) {
                ++zero_count;
                ++other_count;
            }
            else if (the_allele == 2) {
                other_count += 2;
            }
        }
        return (zero_count < other_count) ? zero_count : other_count;
    }
    
    public static void null_h(final String input, final String output_prefix, final int simulation_round, final String pheno_type, final String pheno_method, final double thvalue) {
        try {
            final VariantsDouble genotype = new VariantsDouble(input);
            final BufferedWriter pheno_h0_bw = new BufferedWriter(new FileWriter(String.valueOf(output_prefix) + "." + pheno_type + ".tsv"));
            final double[][] sim_all_pheno = new double[simulation_round][genotype.sample_size];
            final NormalDistribution normal = new NormalDistribution(0.0, Math.sqrt(1.0));
            final UniformRealDistribution uniq = new UniformRealDistribution(0.0, 1.0);
            int round = 0;
            while (round < simulation_round) {
                final double[] sim_pheno_one_round = new double[genotype.sample_size];
                if (pheno_type.equals("Q")) {
                    for (int sample_index = 0; sample_index < sim_pheno_one_round.length; ++sample_index) {
                        sim_pheno_one_round[sample_index] = normal.sample();
                    }
                    sim_all_pheno[round] = sim_pheno_one_round;
                    ++round;
                }
                else if (pheno_type.equals("B")) {
                    for (int sample_index = 0; sample_index < sim_pheno_one_round.length; ++sample_index) {
                        if (pheno_method.equals("T")) {
                            final double tmp = uniq.sample();
                            if (Double.compare(tmp, thvalue) > 0) {
                                sim_pheno_one_round[sample_index] = 2.0;
                            }
                            else {
                                sim_pheno_one_round[sample_index] = 1.0;
                            }
                        }
                        else if (pheno_method.equals("OR")) {
                            System.out.println("Odds Ratio method to generate simulated phenotype will come soon. Please use Threshold for now");
                            System.exit(0);
                        }
                    }
                    if (check_control_case_balance(sim_pheno_one_round)) {
                        sim_all_pheno[round] = sim_pheno_one_round;
                        ++round;
                    }
                    else {
                        System.out.println("Round " + Integer.toString(round) + " does not have a balanced case and control, simulate again");
                    }
                }
                else {
                    System.out.println("Unrecognized phenotype type. Exit");
                    System.exit(0);
                }
            }
            String header = "#ID";
            for (int round_index = 0; round_index < simulation_round; ++round_index) {
                header = String.valueOf(header) + "\t" + Integer.toString(round_index);
            }
            pheno_h0_bw.write(String.valueOf(header) + "\n");
            for (int sample_index = 0; sample_index < genotype.sample_size; ++sample_index) {
                String sample_ID_pheno = genotype.sample_ids[sample_index];
                for (int round_index2 = 0; round_index2 < simulation_round; ++round_index2) {
                    sample_ID_pheno = String.valueOf(sample_ID_pheno) + "\t" + Double.toString(sim_all_pheno[round_index2][sample_index]);
                }
                pheno_h0_bw.write(String.valueOf(sample_ID_pheno) + "\n");
            }
            pheno_h0_bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public static double[] generate_phenotype_from_random_term(final VariantsDouble genotype, final int random_term_snps_num, final String[] sim_all_x_chr_pos, final int[] sim_x_num_per_chr, final HashMap<String, Double> snps_chr_pos_weights_map, final String output_random_term_prefix, final String verbose) {
        final Random rand = new Random();
        final double[][] z_genotype = new double[random_term_snps_num][genotype.sample_size];
        final String[] z_snps_info = new String[random_term_snps_num];
        final double[] z_weights = new double[random_term_snps_num];
        final double[] phenotype_random_term = new double[genotype.sample_size];
        final NormalDistribution STnormal = new NormalDistribution(0.0, 1.0);
        Arrays.fill(phenotype_random_term, 0.0);
        final ArrayList<String> sim_all_x_chr_pos_arrayArrayList = new ArrayList<String>();
        int random_term_snps_count = 0;
        for (final String tmp : sim_all_x_chr_pos) {
            sim_all_x_chr_pos_arrayArrayList.add(tmp);
        }
        for (int chr_index = 0; chr_index < sim_x_num_per_chr.length; ++chr_index) {
            final ArrayList<String> selected_z_snp_info_per_chr = new ArrayList<String>();
            while (selected_z_snp_info_per_chr.size() < sim_x_num_per_chr[chr_index]) {
                final int z_pos_index = rand.nextInt(genotype.locations[chr_index].length);
                final int z_chr = chr_index + 1;
                final String z_pos_ID = String.valueOf(Integer.toString(z_chr)) + "," + Integer.toString(genotype.locations[chr_index][z_pos_index]);
                if (!sim_all_x_chr_pos_arrayArrayList.contains(z_pos_ID) && snps_chr_pos_weights_map.containsKey(z_pos_ID)) {
                    z_genotype[random_term_snps_count] = genotype.load_one_variant_by_index(chr_index, z_pos_index);
                    z_weights[random_term_snps_count] = snps_chr_pos_weights_map.get(z_pos_ID) + STnormal.sample(); //add noise to gene expression weights
                    z_snps_info[random_term_snps_count] = z_pos_ID;
                    ++random_term_snps_count;
                    selected_z_snp_info_per_chr.add(z_pos_ID);
                }
            }
        }
        if (random_term_snps_count == random_term_snps_num) {
            try {
                if (verbose.equals("true")) {
                    final BufferedWriter geno_z_h1_bw = new BufferedWriter(new FileWriter(String.valueOf(output_random_term_prefix) + ".csv"));
                    final BufferedWriter weights_z_h1_bw = new BufferedWriter(new FileWriter(String.valueOf(output_random_term_prefix) + ".weights.txt"));
                    final BufferedWriter weighted_geno_z_h1_bw = new BufferedWriter(new FileWriter(String.valueOf(output_random_term_prefix) + ".weighted.csv"));
                    String header = "SNP,POS";
                    for (int sample_index = 0; sample_index < genotype.sample_size; ++sample_index) {
                        header = String.valueOf(header) + "," + genotype.sample_ids[sample_index];
                    }
                    geno_z_h1_bw.write(String.valueOf(header) + "\n");
                    weights_z_h1_bw.write("SNP,POS,Weights\n");
                    weighted_geno_z_h1_bw.write(String.valueOf(header) + "\n");
                    for (int snp_index = 0; snp_index < random_term_snps_num; ++snp_index) {
                        String snp_info = z_snps_info[snp_index];
                        String weighted_snp_info = z_snps_info[snp_index];
                        for (int sample_index2 = 0; sample_index2 < genotype.sample_size; ++sample_index2) {
                            snp_info = String.valueOf(snp_info) + "," + z_genotype[snp_index][sample_index2];
                            phenotype_random_term[sample_index2] += z_genotype[snp_index][sample_index2] * z_weights[snp_index];
                            weighted_snp_info = String.valueOf(weighted_snp_info) + "," + z_genotype[snp_index][sample_index2] * z_weights[snp_index];
                        }
                        geno_z_h1_bw.write(String.valueOf(snp_info) + "\n");
                        weights_z_h1_bw.write(String.valueOf(z_snps_info[snp_index]) + "," + z_weights[snp_index] + "\n");
                        weighted_geno_z_h1_bw.write(String.valueOf(weighted_snp_info) + "\n");
                    }
                    geno_z_h1_bw.close();
                    weights_z_h1_bw.close();
                    weighted_geno_z_h1_bw.close();
                }
                else {
                    for (int snp_index2 = 0; snp_index2 < random_term_snps_num; ++snp_index2) {
                        for (int sample_index3 = 0; sample_index3 < genotype.sample_size; ++sample_index3) {
                            final double[] array = phenotype_random_term;
                            final int n = sample_index3;
                            array[n] += z_genotype[snp_index2][sample_index3] * z_weights[snp_index2];
                        }
                    }
                }
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        }
        else {
            System.out.println("Randome term matrix is not full, with SNPs number " + Integer.toString(random_term_snps_count) + "But we need " + Integer.toString(random_term_snps_num));
            System.exit(0);
        }
        return phenotype_random_term;
    }
    
    public static double calculate_sim_pheno_one_round_mean(final double[] sim_pheno_one_round) {
        final int sample_size = sim_pheno_one_round.length;
        double sim_pheno_one_round_mean = 0.0;
        for (int sample_index = 0; sample_index < sample_size; ++sample_index) {
            sim_pheno_one_round_mean += sim_pheno_one_round[sample_index];
        }
        return sim_pheno_one_round_mean / sample_size;
    }
    
    public static void alter_h(final String input, final double maf, final int simulation_round, final String weight_filename, final int random_term_snps_num, final String output_prefix, final String pheno_type, final String pheno_method, final double h2_1, final double h2_2, final String verbose, final double thvalue) {
        try {
            final VariantsDouble genotype = new VariantsDouble(input);
            final BufferedReader pheno_snps_weights_br = new BufferedReader(new FileReader(weight_filename));
            final BufferedWriter geno_h1_bw = new BufferedWriter(new FileWriter(String.valueOf(output_prefix) + "." + pheno_type + "." + random_term_snps_num + ".h2_1." + h2_1 + ".h2_2." + h2_2 + ".csv"));
            final BufferedWriter pheno_h1_bw = new BufferedWriter(new FileWriter(String.valueOf(output_prefix) + "." + pheno_type + "." + random_term_snps_num + ".h2_1." + h2_1 + ".h2_2." + h2_2 + ".tsv"));
            final BufferedWriter geno_h1_log_bw = new BufferedWriter(new FileWriter(String.valueOf(output_prefix) + ".csv.log"));
            final BufferedWriter pheno_h1_log_bw = new BufferedWriter(new FileWriter(String.valueOf(output_prefix) + "." + pheno_type + "." + random_term_snps_num + ".h2_1." + h2_1 + ".h2_2." + h2_2 + ".tsv.log"));
            final String[] sim_all_x_chr_pos = new String[simulation_round];
            final int[] sim_x_num_per_chr = new int[22];
            final double[][] sim_all_x = new double[simulation_round][genotype.sample_size];
            final double[][] sim_all_pheno = new double[genotype.sample_size][simulation_round];
            final HashMap<String, Double> snps_chr_pos_weights_map = new HashMap<String, Double>();
            int round = 0;
            final int[][] genotype_locations = genotype.locations;
            final NormalDistribution weights_error_from_nromal = new NormalDistribution(0.0, Math.sqrt(1.0));
            final Random rand = new Random();
            //Select focal and non-focal SNPs
            while (round < simulation_round) {
                final int random_chr = rand.nextInt(22);
                final int[] random_pos = genotype_locations[random_chr];
                final int random_pos_index = rand.nextInt(random_pos.length);
                final double[] selected_x = genotype.load_one_variant_by_index(random_chr, random_pos_index);
                final double minor_allels_frequency = mafc(selected_x) / (double)(2 * selected_x.length);
                if (Double.compare(minor_allels_frequency, maf) >= 0) {
                    if (selected_x.length == genotype.sample_size) {
                        sim_all_x[round] = selected_x;
                        sim_all_x_chr_pos[round] = String.valueOf(Integer.toString(random_chr + 1)) + "," + Integer.toString(random_pos[random_pos_index]);
                        ++round;
                    }
                    else {
                        System.out.println("Selected SNP does not have the same sample number as genotype.\nPlease Double check your genotype.\nExit");
                        System.exit(0);
                    }
                    geno_h1_log_bw.write(String.valueOf(Integer.toString(1 + random_chr)) + "," + Integer.toString(random_pos[random_pos_index]) + " SNP maf is " + Double.toString(minor_allels_frequency) + "\n");
                }
                else {
                    final int random_chr_value = random_chr + 1;
                    geno_h1_log_bw.write(String.valueOf(Integer.toString(random_chr_value)) + "," + Integer.toString(random_pos[random_pos_index]) + " SNP maf is " + Double.toString(minor_allels_frequency) + ", less than maf threshold " + Double.toString(maf) + ", will select again\n");
                }
            }
            geno_h1_log_bw.close();
            //Output genoytpes used for simulations
            final String[] sampleID = genotype.sample_ids;
            String header = "CHR,POS";
            for (int sample_index = 0; sample_index < sampleID.length; ++sample_index) {
                header = String.valueOf(header) + "," + sampleID[sample_index];
            }
            geno_h1_bw.write(String.valueOf(header) + "\n");
            for (int round_index = 0; round_index < simulation_round; ++round_index) {
                String selected_x_info = sim_all_x_chr_pos[round_index];
                final double[] selected_x_geno = sim_all_x[round_index];
                for (int sample_index2 = 0; sample_index2 < sampleID.length; ++sample_index2) {
                    selected_x_info = String.valueOf(selected_x_info) + "," + Double.toString(selected_x_geno[sample_index2]);
                }
                geno_h1_bw.write(String.valueOf(selected_x_info) + "\n");
            }
            geno_h1_bw.close();
            final int[] genotype_snps_num_per_chr = genotype.num_sites;
            int genotype_snps_num_total = 0;
            for (int i = 0; i < genotype_snps_num_per_chr.length; ++i) {
                genotype_snps_num_total += genotype_snps_num_per_chr[i];
            }
            int sim_x_num_per_chr_total = 0;
            for (int j = 0; j < sim_x_num_per_chr.length; ++j) {
                sim_x_num_per_chr[j] = (int)(random_term_snps_num * (genotype_snps_num_per_chr[j] / (double)genotype_snps_num_total));
                sim_x_num_per_chr_total += sim_x_num_per_chr[j];
            }
            final int diff = random_term_snps_num - sim_x_num_per_chr_total;
            final int[] array = sim_x_num_per_chr;
            final int n = 0;
            array[n] += diff;
            for (int item_index = 0; item_index < sim_x_num_per_chr.length; ++item_index) {
                pheno_h1_log_bw.write("chr" + (item_index + 1) + ":" + sim_x_num_per_chr[item_index] + "\n");
            }
            //load SNPs weights from weights file
            pheno_snps_weights_br.readLine();
            for (String weight_line = pheno_snps_weights_br.readLine(); weight_line != null; weight_line = pheno_snps_weights_br.readLine()) {
                final String[] weight_line_arr = weight_line.split(" ");
                final String weight_key = String.valueOf(weight_line_arr[1]) + "," + weight_line_arr[0];
                final double weight_value = Math.abs(Double.parseDouble(weight_line_arr[3]) + Double.parseDouble(weight_line_arr[4]))+ weights_error_from_nromal.sample(); //add one error from normal distribution 
                if (!snps_chr_pos_weights_map.containsKey(weight_key)) {
                    snps_chr_pos_weights_map.put(weight_key, weight_value);
                }
                else {
                    final double old_weight_value = snps_chr_pos_weights_map.get(weight_key);
                    if (Double.compare(weight_value, old_weight_value) > 0) {
                        snps_chr_pos_weights_map.put(weight_key, weight_value);
                    }
                }
            }
            pheno_snps_weights_br.close();
            //Simulate phenotypes based on weights, selected geneotypes
            int round_count = 0;
            while (round_count < simulation_round) {
                final String output_random_term_prefix = String.valueOf(output_prefix) + ".s" + Integer.toString(round_count) + ".h2_1." + Double.toString(h2_1) + ".h2_2." + Double.toString(h2_2) + ".z";
                final double[] phenotype_random_term = generate_phenotype_from_random_term(genotype, random_term_snps_num, sim_all_x_chr_pos, sim_x_num_per_chr, snps_chr_pos_weights_map, output_random_term_prefix, verbose);
                final double variance_g = StatFuncs.popVar_NaN(sim_all_x[round_count]);
                final double variance_z = StatFuncs.popVar_NaN(phenotype_random_term);
                final double b = Math.sqrt((h2_1 - h2_2) / h2_2 * (variance_g / variance_z));
                final double c = Math.sqrt(variance_g * (1.0 - h2_1) / h2_2);
                final double variance_e = 1.0;
                final double[] sim_pheno_one_round = new double[genotype.sample_size];
                pheno_h1_log_bw.write("Round " + round_count + ", variance_g " + variance_g + ", variance_z " + variance_z + ", variance_e " + variance_e + ", h2_1 " + h2_1 + ", h2_2 " + h2_2 + ", b " + b + ", c " + c + "\n");
                final NormalDistribution normal = new NormalDistribution(0.0, Math.sqrt(variance_e));
                for (int sample_index3 = 0; sample_index3 < genotype.sample_size; ++sample_index3) {
                    sim_pheno_one_round[sample_index3] = 1.0 * sim_all_x[round_count][sample_index3] + b * phenotype_random_term[sample_index3] + c * normal.sample();
                }
                if (pheno_type.equals("B")) {
                    if (pheno_method.equals("T")) {
                        final double sim_pheno_one_round_mean = calculate_sim_pheno_one_round_mean(sim_pheno_one_round);
                        for (int sample_index4 = 0; sample_index4 < genotype.sample_size; ++sample_index4) {
                            if (Double.compare(sim_pheno_one_round[sample_index4], sim_pheno_one_round_mean) > 0) {
                                sim_pheno_one_round[sample_index4] = 2.0;
                            }
                            else {
                                sim_pheno_one_round[sample_index4] = 1.0;
                            }
                        }
                        if (!check_control_case_balance(sim_pheno_one_round)) {
                            continue;
                        }
                        for (int sample_index4 = 0; sample_index4 < genotype.sample_size; ++sample_index4) {
                            sim_all_pheno[sample_index4][round_count] = sim_pheno_one_round[sample_index4];
                        }
                        ++round_count;
                    }
                    else {
                        if (!pheno_method.equals("OR")) {
                            continue;
                        }
                        System.out.println("Odds Ratio method to generate simulated phenotype will come soon. Please use Threshold for now");
                        System.exit(0);
                    }
                }
                else if (pheno_type.equals("Q")) {
                    for (int sample_index3 = 0; sample_index3 < genotype.sample_size; ++sample_index3) {
                        sim_all_pheno[sample_index3][round_count] = sim_pheno_one_round[sample_index3];
                    }
                    ++round_count;
                }
                else {
                    System.out.println("Unrecognized phenotype type. Exit");
                    System.exit(0);
                }
            }
            header = "SampleID";
            for (int round_index2 = 0; round_index2 < simulation_round; ++round_index2) {
                header = String.valueOf(header) + "\tR" + round_index2;
            }
            pheno_h1_bw.write(String.valueOf(header) + "\n");
            for (int sample_index5 = 0; sample_index5 < genotype.sample_size; ++sample_index5) {
                String sample_sim_pheno = genotype.sample_ids[sample_index5];
                for (int round_index3 = 0; round_index3 < simulation_round; ++round_index3) {
                    sample_sim_pheno = String.valueOf(sample_sim_pheno) + "\t" + sim_all_pheno[sample_index5][round_index3];
                }
                pheno_h1_bw.write(String.valueOf(sample_sim_pheno) + "\n");
            }
            pheno_h1_log_bw.close();
            pheno_h1_bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public static void main(final String[] args) {
        if (args.length == 0) {
            System.out.println(SimedLMM.edLMMIntro);
            System.out.println("\nSupported functions:" + SimedLMM.supported_functions);
            System.exit(0);
        }
        final String function = args[0];
        if (function.equals("null_h")) {
            if (args.length == 1) {
                System.out.println("Simulation under null hypothsis to estimate Type I Error");
                System.out.println("Usage: \n\t<-ig\tinput_genotype_file in hdf5>\n\t<-r\tsimulation round (df=1000)>\n\t<-o\toutput_prefix>\n\t<-t\tphenotype type (B:Binary; Q:quantitative, df=Q)>\n\t<-m\tmethod to generate Binary phenotype (T: Threshold; OR: OddsRatio, df=T)>\n\t[-tv\tthreshold value for phenotype to set to 2/control (above tv) or 1/case (below or equal to tv), df=0.5]");
                System.exit(0);
            }
            else {
                String input = null;
                String output_prefix = null;
                int simulation_round = 1000;
                String pheno_type = "Q";
                String pheno_method = "T";
                double thvalue = 0.5;
                for (int k = 1; k < args.length; ++k) {
                    if (args[k].startsWith("-")) {
                        if (args[k].equals("-ig")) {
                            input = args[k + 1];
                        }
                        else if (args[k].equals("-o")) {
                            output_prefix = args[k + 1];
                        }
                        else if (args[k].equals("-r")) {
                            simulation_round = Integer.parseInt(args[k + 1]);
                        }
                        else if (args[k].equals("-t")) {
                            pheno_type = args[k + 1];
                        }
                        else if (args[k].equals("-m")) {
                            pheno_method = args[k + 1];
                        }
                        else if (args[k].equals("-tv")) {
                            thvalue = Double.parseDouble(args[k + 1]);
                        }
                    }
                }
                if (input == null || output_prefix == null) {
                    System.out.println("Input or output-prefix can't be null!");
                }
                else {
                    System.out.println(input);
                    System.out.println(output_prefix);
                    System.out.println(pheno_type);
                    null_h(input, output_prefix, simulation_round, pheno_type, pheno_method, thvalue);
                }
            }
        }
        else if (function.equals("alter_h")) {
            if (args.length == 1) {
                System.out.println("Simulation under alternative hypothsis to calculate power");
                System.out.println("Usage: \n\t<-ig\tinput_genotype_file in hdf5>\n\t<-maf\tminor allele frequency above which SNPs will be selected (df=0.2)>\n\t<-r\tsimulation round (df=1000)>\n\t<-wf\tweights file to generate randome term>\n\t<-wm\tnumber of SNPs used to generate random term (df=15000)>\n\t<-o\toutput_prefix>\n\t<-t\tphenotype type (B:Binary; Q:quantitative, df=Q)\n\t<-m\tmethod to generate Binary phenotype (T: Threshold; OR: OddsRatio, df=T)>\n\t<-h2_1\toverall heritablity ((0-1), df=0.5)>\n\t<-h2_2\theritablity of focal variant x ((0-0.5), df=0.001)>\n\t<-v\t verbose (true OR false, output detailed information. E.g., random term genotype, df=false)>\n\t[-tv\tthreshold value for phenotype to set to 2/control (above tv) or 1/case (below or equal to tv), df=mean(pheno)]\n\t");
                System.exit(0);
            }
            else {
                String input = null;
                String output_prefix = null;
                String weight_filename = null;
                String verbose = "false";
                int simulation_round2 = 1000;
                int random_term_snps_num = 15000;
                double maf = 0.2;
                double h2_1 = 0.5;
                double h2_2 = 0.001;
                double thvalue2 = 0.0;
                String pheno_type2 = "Q";
                String pheno_method2 = "T";
                for (int i = 1; i < args.length; ++i) {
                    if (args[i].startsWith("-")) {
                        if (args[i].equals("-ig")) {
                            input = args[i + 1];
                        }
                        else if (args[i].equals("-maf")) {
                            maf = Double.parseDouble(args[i + 1]);
                        }
                        else if (args[i].equals("-r")) {
                            simulation_round2 = Integer.parseInt(args[i + 1]);
                        }
                        else if (args[i].equals("-wf")) {
                            weight_filename = args[i + 1];
                        }
                        else if (args[i].equals("-wm")) {
                            random_term_snps_num = Integer.parseInt(args[i + 1]);
                        }
                        else if (args[i].equals("-o")) {
                            output_prefix = args[i + 1];
                        }
                        else if (args[i].equals("-t")) {
                            pheno_type2 = args[i + 1];
                        }
                        else if (args[i].equals("-m")) {
                            pheno_method2 = args[i + 1];
                        }
                        else if (args[i].equals("-h2_1")) {
                            h2_1 = Double.parseDouble(args[i + 1]);
                        }
                        else if (args[i].equals("-h2_2")) {
                            h2_2 = Double.parseDouble(args[i + 1]);
                        }
                        else if (args[i].equals("-v")) {
                            verbose = args[i + 1];
                        }
                        else if (args[i].equals("-tv")) {
                            thvalue2 = Double.parseDouble(args[i + 1]);
                        }
                    }
                }
                if (input == null || output_prefix == null) {
                    System.out.println("Input or output-prefix can't be null!");
                }
                else {
                    alter_h(input, maf, simulation_round2, weight_filename, random_term_snps_num, output_prefix, pheno_type2, pheno_method2, h2_1, h2_2, verbose, thvalue2);
                }
            }
        }
        
    }
}
