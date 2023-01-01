package ccg_unidirectional;

import java.io.File;
import java.io.IOException;

import ilog.concert.IloException;
import ilog.concert.IloNumExpr;
import ilog.concert.IloNumVar;
import ilog.concert.IloNumVarType;

import ilog.concert.IloRange;
import ilog.cplex.IloCplex;


public class StochasticProgrammingModel {
	IloCplex upperPro;
	IloCplex[] lowerProList;
	
	//上层问题变量
	IloNumVar[] CmaxKsiUpper;
	double[] CmaxKsiUpperValue;
	IloNumVar[][] Xhk;
	double[][] XhkValue;
	IloNumVar[] alpha;
	double[] alphaValue;
	IloNumVar[] beta;
	double[] betaValue;
	// 下层问题变量
	IloNumVar[] CmaxKsiLower;
	IloNumVar[][][] YKsiHK;
	double[][][] YksiHKValue;
	IloNumVar[][][] QKsiJK;
	double[][][] QKsiJKValue;
	IloNumVar[][] RKsiK;
	double[][] RKsiKValue;
	// 上层问题约束
	IloRange[] upper4_2;
	IloRange[][] upper4_3;
	IloRange[][] upper4_4;
	IloRange[] upper4_8;
	IloRange[] upper4_9;
	// 下层问题约束
	IloRange[][][] lower4_12;
	IloRange[][] lower4_13;
	// 下层问题参数矩阵
	double[][] matrixT;
	double[] vectorH;
	// 下层问题对偶变量
	double[][] pai_value;
	
	long startTime; // 算法运行开始时间
	long endTime; // 算法运行结束时间
	private static final double M = 3000;
	private static final double ACCURANCY = 0.1; // C&CG算法终止的精度
	private static final double FLOATACCURACY = 0.001;
	double terminalSpan = 1800; // 程序最大运行时间，以s为单位
	int Iteration = 1;
	int jValue; // 贝位数减去安全距离再减1
	int numScenarioByBay = 2; // 每一个贝位作业时间的类型数量
	double[] probabilityKsi; // 不同scenario对应的概率
	
	public void upperModel(Input input, boolean variableOutput, boolean iterationDetail) throws IloException, IOException {
		System.out.println("上层模型开始求解");
		int num_ksi = (int) Math.round(Math.pow(numScenarioByBay, input.numBay));
		startTime = System.currentTimeMillis();
		upperPro = new IloCplex();
		// 变量
		Xhk = new IloNumVar[input.numBay][input.numQC];
		CmaxKsiUpper = new IloNumVar[num_ksi];
		alpha = new IloNumVar[input.numQC];
		beta = new IloNumVar[input.numQC];
		// 下层问题对偶变量
		
		// 约束
		upper4_2 = new IloRange[input.numBay];
		upper4_3 = new IloRange[input.numBay][input.numQC];
		upper4_4 = new IloRange[input.numBay][input.numQC];
		upper4_8 = new IloRange[input.numQC];
		upper4_9 = new IloRange[input.numQC];
		
		// 设置变量
		for (int i = 0; i < input.numBay; i++) {
			for (int j = 0; j < input.numQC; j++) {
				Xhk[i][j] = upperPro.numVar(0, 1, IloNumVarType.Int, "X" + i + j);
			}
		}
		
		for (int i = 0; i < num_ksi; i++) {
			CmaxKsiUpper[i] = upperPro.numVar(0, Double.MAX_VALUE, IloNumVarType.Float, "Cmax" + i);
		}
		
		for (int i = 0; i < input.numQC; i++) {
			alpha[i] = upperPro.numVar(1, input.numBay, IloNumVarType.Float, "alpha" + i);
		}

		for (int i = 0; i < input.numQC; i++) {
			beta[i] = upperPro.numVar(1, input.numBay, IloNumVarType.Float, "beta" + i);
		}
		
		// 设置目标函数
		IloNumExpr upperObj = upperPro.numExpr();
		for (int i = 0; i < num_ksi; i++) {
			upperObj = upperPro.sum(upperObj, upperPro.prod(CmaxKsiUpper[i], 1));
		}
		upperPro.addMinimize(upperObj);
		
		// 上层问题约束
		// 约束4-2
		for (int i = 0; i < input.numBay; i++) {
			if (input.processTime[i] != 0) {
				IloNumExpr expr = upperPro.numExpr();
				for (int j = 0; j < input.numQC; j++) {
					expr = upperPro.sum(expr, upperPro.prod(Xhk[i][j], 1));
				}
				upper4_2[i] = upperPro.addEq(expr, 1, "constraint4_2 " + i);
			}
		}

		// 约束4-3
		for (int i = 0; i < input.numBay; i++) {
			for (int j = 0; j < input.numQC; j++) {
				IloNumExpr expr = upperPro.numExpr();
				expr = upperPro.sum(expr, upperPro.prod(Xhk[i][j], M - i - 1));
				expr = upperPro.sum(expr, upperPro.prod(alpha[j], 1));
				upper4_3[i][j] = upperPro.addLe(expr, M, "constraint4_3 " + i + "," + j);
			}
		}

		// 约束4-4
		for (int i = 0; i < input.numBay; i++) {
			for (int j = 0; j < input.numQC; j++) {
				IloNumExpr expr = upperPro.numExpr();
				expr = upperPro.sum(expr, upperPro.prod(Xhk[i][j], -M - i - 1));
				expr = upperPro.sum(expr, upperPro.prod(beta[j], 1));
				upper4_4[i][j] = upperPro.addGe(expr, -M, "constraint4_4 " + i + "," + j);
			}
		}

		// 约束4-8
		for (int i = 0; i < input.numQC - 1; i++) {
			IloNumExpr expr = upperPro.numExpr();
			expr = upperPro.sum(expr, upperPro.prod(alpha[i], 1));
			expr = upperPro.sum(expr, upperPro.prod(alpha[i + 1], -1));
			upper4_8[i] = upperPro.addLe(expr, -input.safetyMargin - 1, "constraint4_8 " + i);
		}

		// 约束4-9
		for (int i = 0; i < input.numQC - 1; i++) {
			IloNumExpr expr = upperPro.numExpr();
			expr = upperPro.sum(expr, upperPro.prod(beta[i], 1));
			expr = upperPro.sum(expr, upperPro.prod(beta[i + 1], -1));
			upper4_9[i] = upperPro.addLe(expr, -input.safetyMargin - 1, "constraint4_9 " + i);
		}
		
		//进行求解
		if (upperPro.solve()) {
			System.out.println("问题求解成功");
			XhkValue = new double[input.numBay][input.numQC];
			CmaxKsiUpperValue = new double[num_ksi];
			alphaValue = new double[input.numQC];
			betaValue = new double[input.numQC];
			for (int i = 0; i < input.numBay; i++) {
				for (int j = 0; j < input.numQC; j++) {
					XhkValue[i][j] = upperPro.getValue(Xhk[i][j]);
				}
			}
			CmaxKsiUpperValue = upperPro.getValues(CmaxKsiUpper);
			alphaValue = upperPro.getValues(alpha);
			betaValue = upperPro.getValues(beta);
			System.out.println("Solution status = " + upperPro.getStatus());
			System.out.println("Solution value  = " + upperPro.getObjValue());
			input.writerOverall.write("Iteration " + Iteration + " upperObj:" + upperPro.getObjValue() + "\r\n");
			
			if (iterationDetail) {
				File file = new File("./result of benchmark/test_" + input.testID);
				if (!file.exists()) {// 如果文件夹不存在
					file.mkdir();// 创建文件夹
				}
				upperPro.writeSolution("./result of benchmark/test_" + input.testID + "/test_" + input.testID
						+ " iteration_" + Iteration + "upper.txt");
			}
			
			if (variableOutput) {
				for (int i = 0; i < input.numBay; i++) {
					for (int j = 0; j < input.numQC; j++) {
						System.out.println("X" + i + "," + j + ": Value = " + XhkValue[i][j]);
					}
				}
				for (int i = 0; i < num_ksi; i++) {
					System.out.println("CmaxKsi" + i + ": Value = " + CmaxKsiUpperValue[i]);
				}
				for (int i = 0; i < input.numQC; i++) {
					System.out.println("alpha" + i + ": Value = " + alphaValue[i]);
				}

				for (int i = 0; i < input.numQC; i++) {
					System.out.println("beta" + i + ": Value = " + betaValue[i]);
				}
			}	
		} else {
			System.out.println("上层问题未成功求解");
			System.exit(0);
		}
	}
	
	public void lowerModelLShape(Input input, boolean variableOutput, boolean iterationDetail) throws IloException {
//		int num_ksi = (int) Math.round(Math.pow(numScenarioByBay, input.numBay)); 
		int num_ksi = 1; 
		lowerProList = new IloCplex[num_ksi];
//		for (int i = 0; i < num_ksi; i++) {
//			lowerProList[i] = new IloCplex();
//		}
		jValue = input.numBay - input.safetyMargin - 1;
		probabilityKsi = new double[num_ksi];
		for (int i = 0; i < num_ksi; i++) {
			probabilityKsi[i] = 1 / num_ksi;
		}
		// 变量
		CmaxKsiLower = new IloNumVar[num_ksi];
		YKsiHK = new IloNumVar[num_ksi][input.numBay][input.numQC];
		QKsiJK = new IloNumVar[num_ksi][jValue][input.numQC - 1];
		RKsiK = new IloNumVar[num_ksi][input.numQC];
		YksiHKValue = new double[num_ksi][input.numBay][input.numQC];
		QKsiJKValue = new double[num_ksi][jValue][input.numQC - 1];
		RKsiKValue = new double[num_ksi][input.numQC];
		
		// 约束
		lower4_12 = new IloRange[num_ksi][jValue][input.numQC - 1];
		lower4_13 = new IloRange[num_ksi][input.numQC];
		
		// 矩阵
		int num_row = jValue * (input.numQC - 1) + input.numQC;
		int num_column = input.numBay * input.numQC + 2 * input.numQC;
		if (matrixT == null) {
			matrixT = new double[num_row][num_column];
		}
		if (vectorH == null) {
			vectorH = new double[num_row];
		}
		// 记录当前行
		int current_row = 0;
		// 对偶变量取值
		pai_value = new double[num_ksi][num_row];
		
		// 按照模型进行建模
		for (int i = 0; i < num_ksi; i++) {
			System.out.println("下层模型第" + i + "个scenario求解");
			lowerProList[i] = new IloCplex();
			// 设置变量
			CmaxKsiLower[i] = lowerProList[i].numVar(0, Double.MAX_VALUE, IloNumVarType.Float, "CmaxKsi_" + i);
			for (int j = 0; j < jValue; j++) {
				for (int j2 = 0; j2 < input.numQC - 1; j2++) {
					QKsiJK[i][j][j2] = lowerProList[i].numVar(0, Double.MAX_VALUE, IloNumVarType.Float, "QKsiJK_" + i + "," + j + "," + j2);
				}
			}
			for (int j = 0; j < input.numQC; j++) {
				RKsiK[i][j] = lowerProList[i].numVar(0, Double.MAX_VALUE, IloNumVarType.Float, "RKsiK_" + i + "," + j);
			}
			for (int j = 0; j < input.numBay; j++) {
				for (int j2 = 0; j2 < input.numQC; j2++) {
					YKsiHK[i][j][j2] = lowerProList[i].numVar(0, Double.MAX_VALUE, IloNumVarType.Float);
				}
			}
			
			// 目标函数
			IloNumExpr obj = lowerProList[i].numExpr();
			obj = lowerProList[i].sum(obj, CmaxKsiLower[i]);
			lowerProList[i].addMinimize(obj);
			
			// 设置约束
			//	4-12
			for (int j = 0; j < jValue; j++) {
				for (int j2 = 0; j2 < input.numQC - 1; j2++) {
					IloNumExpr expr = lowerProList[i].numExpr();
					double left_constant = 0, right_constant = 0;
					
					for (int k = 0; k <= j; k++) {
						left_constant += input.processTime[k] * XhkValue[k][j2];
						if (i == 0) {
							matrixT[current_row][k * input.numQC + j2] += input.processTime[k];
						}
						expr = lowerProList[i].sum(expr, YKsiHK[i][k][j2]);
					}
					expr = lowerProList[i].sum(expr, lowerProList[i].prod(QKsiJK[i][j][j2], -1));
					left_constant += (j + 1 - alphaValue[j2]) * input.traverseTime;
					if (i == 0) {
						matrixT[current_row][input.numBay * input.numQC + j2] += -input.traverseTime;
						vectorH[current_row] += -(j + 1) * input.traverseTime;
					}
					
					for (int k = input.safetyMargin + 1; k <= j + input.safetyMargin + 1 ; k++) {
						right_constant += input.processTime[k] * XhkValue[k][j2 + 1];
						if (i == 0) {
							matrixT[current_row][k * input.numQC + j2 + 1] = - input.processTime[k];
						}
						expr = lowerProList[i].sum(expr, lowerProList[i].prod(YKsiHK[i][k][j2 + 1], -1));
					}
					right_constant += (j + 1 + input.safetyMargin + 1 - alphaValue[j2 + 1]) * input.traverseTime;
					if (i == 0) {
						matrixT[current_row][input.numBay * input.numQC + j2 + 1] += input.traverseTime;
						vectorH[current_row] += (j + 1 + input.safetyMargin + 1) * input.traverseTime;
						current_row++;
					}
					lower4_12[i][j][j2] = lowerProList[i].addEq(expr, right_constant - left_constant, "Constraint4_12:" + i +"," + j + "," + j2);
				}
			}
			// 4-13
			for (int j = 0; j < input.numQC; j++) {
				double left_constant = 0;
				IloNumExpr expr = lowerProList[i].numExpr();
				for (int k = 0; k < input.numBay; k++) {
					left_constant += input.processTime[k] * XhkValue[k][j];
					if (i == 0) {
						matrixT[current_row][k * input.numQC + j] += input.processTime[k];
					}
					expr = lowerProList[i].sum(expr, YKsiHK[i][k][j]);
				}
				left_constant += (betaValue[j] - alphaValue[j]) * input.traverseTime;
				if (i == 0) {
					matrixT[current_row][input.numBay * input.numQC + j] += -input.traverseTime;
					matrixT[current_row][input.numBay * input.numQC + input.numQC + j] += input.traverseTime;
					current_row++;
				}
				expr = lowerProList[i].sum(expr, RKsiK[i][j]);
				expr = lowerProList[i].sum(expr, lowerProList[i].prod(-1, CmaxKsiLower[i]));
				lower4_13[i][j] = lowerProList[i].addLe(expr, -left_constant, "Contraint4_13" + i +"," + j);
			}
			
			// 求解问题
			if (lowerProList[i].solve()) {
				endTime = System.currentTimeMillis();
				System.out.println("Solution status = " + lowerProList[i].getStatus());
				System.out.println("Solution value  = " + lowerProList[i].getObjValue());
				
				if (variableOutput) {
					
					for (int j = 0; j < input.numBay; j++) {
						for (int j2 = 0; j2 < input.numQC; j2++) {
							YksiHKValue[i][j][j2] = lowerProList[i].getValue(YKsiHK[i][j][j2]);
							System.out.println("YKsiHK" + i + "," + j + "," + j2 + "Value: " + YksiHKValue[i][j][j2]);
						}
					}
					
					for (int j = 0; j < jValue; j++) {
						for (int j2 = 0; j2 < input.numQC - 1; j2++) {
							QKsiJKValue[i][j][j2] = lowerProList[i].getValue(QKsiJK[i][j][j2]);
							System.out.println("QKsiJK" + i + "," + j + "," + j2 + "Value: " + QKsiJKValue[i][j][j2]);
						}
					}
					
					for (int j = 0; j < input.numQC; j++) {
						RKsiKValue[i][j] = lowerProList[i].getValue(RKsiK[i][j]);
						System.out.println("RKsiK" + i + "," + j + "Value: " + RKsiKValue[i][j]);
					}
					
				}
				
				

//				if (iteration) {
//					input.writerOverall
//							.write("Iteration " + Iteration + " lowerObj(SD):" + lowerPro.getObjValue() + "\r\n");
//					input.writerOverall.write("\r\n");
//
//					File file = new File("./result of benchmark/test_" + input.testID);
//					if (!file.exists()) {// 如果文件夹不存在
//						file.mkdir();// 创建文件夹
//					}
//
//					if (iterationDetail) {
//						lowerPro.writeSolution("./result of benchmark/test_" + input.testID + "/test_" + input.testID
//								+ " iteration_" + Iteration + "lowerSD.txt");
//					}
//				}
				
				// 获取对偶变量取值
				// 约束4_12
				for (int j = 0; j < jValue; j++) {
					for (int j2 = 0; j2 < input.numQC - 1; j2++) {
						pai_value[i][j * (input.numQC - 1) + j2] = lowerProList[i].getDual(lower4_12[i][j][j2]);
						System.out.println("pai_" + i + '_' + (j * (input.numQC - 1) + j2) + ":" + pai_value[i][j * (input.numQC - 1) + j2]);
					}
				}
				// 约束4_13
				for (int j = 0; j < input.numQC; j++) {
					pai_value[i][jValue * (input.numQC - 1) + j] = lowerProList[i].getDual(lower4_13[i][j]);
					System.out.println("pai_" + i + '_' + (jValue * (input.numQC - 1) + j) + ":" + pai_value[i][jValue * (input.numQC - 1) + j]);
				}

			} else {
				System.out.println("下层问题求解失败");
				System.exit(0);
			}
		}
		
		
	}
}
