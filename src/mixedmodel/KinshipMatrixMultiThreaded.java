package mixedmodel;

import ch.systemsx.cisd.base.mdarray.MDDoubleArray;
import ch.systemsx.cisd.hdf5.HDF5MDDataBlock;
import org.jfree.util.Log;

import java.util.Random;

public class KinshipMatrixMultiThreaded extends Thread {
    private VariantsDouble data1;
    private int chr1;
    private double scale1;
    private double min_MAF1;
    private double[][] kinship1;
    private double[][] total_num_var_used1;

    public KinshipMatrixMultiThreaded(VariantsDouble data, int chr, double scale, double min_MAF, double[][] kinship, double[][] total_num_var_used) {
        this.data1 = data;
        this.chr1 = chr;
        this.scale1 = scale;
        this.min_MAF1 = min_MAF;
        this.kinship1 = kinship;
        this.total_num_var_used1 = total_num_var_used;
    }

    public static int number_of_cores() {
        return Runtime.getRuntime().availableProcessors();
    }

    public void run() {
        process_chr(data1, chr1, scale1, min_MAF1, kinship1, total_num_var_used1);
    }

    private synchronized void process_chr(VariantsDouble data, int chr, double scale, double min_MAF, double[][] kinship, double[][] total_num_var_used) {
        System.out.println("Working on Chr" + (chr + 1));
        for (HDF5MDDataBlock<MDDoubleArray> block : data.position_fast_blocks[chr]) {
            double[][] data4thisblock = block.getData().toMatrix();
            for (int var_index = 0; var_index < data4thisblock.length; var_index++) {
                if (KinshipMatrix.maf(data4thisblock[var_index], scale) < min_MAF) continue;
                for (int i = 0; i < data.sample_size; i++) {
                    for (int j = i + 1; j < data.sample_size; j++) {
                        if ((!Double.isNaN(data4thisblock[var_index][i])) && (!Double.isNaN(data4thisblock[var_index][j]))) {
                            kinship[i][j] = kinship[i][j] + (scale - Math.abs(data4thisblock[var_index][i] - data4thisblock[var_index][j]));
                            total_num_var_used[i][j]++;
                        }

                    }
                }
            }
        }
        System.out.println("Finished Chr" + (chr + 1));
    }
}
