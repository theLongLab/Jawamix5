package real_data_analysis;

import myMathLib.Test;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;

public class TryPorgrams2 {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		
		double[][] data={{1,2},{3,4}};
		
		double[][] add={{1,1},{1,1}};
		RealMatrix kinship_matrix=new Array2DRowRealMatrix(data, false);
		RealMatrix reginal_kinship=new Array2DRowRealMatrix(add, false);
		Test.outputArray(kinship_matrix.getData());
		kinship_matrix=kinship_matrix.add(reginal_kinship.scalarMultiply(100));
		Test.outputArray(reginal_kinship.getData());
		Test.outputArray(kinship_matrix.getData());

	}

}
