package ccg_unidirectional;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import ilog.concert.IloException;
import ilog.concert.IloNumExpr;
import ilog.concert.IloNumVar;
import ilog.concert.IloNumVarType;
import ilog.concert.IloObjective;
import ilog.concert.IloRange;
import ilog.cplex.IloCplex;

public class Model {
	IloCplex upperPro;
	IloCplex lowerPro;
	IloCplex determinsticAttackPro; // 用于确定性模型进行破坏的计算
	IloObjective lowerProObj;
	// 上层问题变量
	IloNumVar[][] Xhk;
	IloNumVar[] CmaxUpper;
	double[][] XhkValue;
	IloNumVar[] alpha;
	double[] alphaValue;
	IloNumVar[] beta;
	double[] betaValue;
	// 下层问题变量
	IloNumVar[][] Yhk;
	IloNumVar[] CmaxLower;
	IloNumVar[] Zh;
	IloNumVar[][] Niujk;
	IloNumVar[] Miuk;
	IloNumVar[][] Gammajk;
	IloNumVar[] Epsilonk;
	IloNumVar[][] Zetahk;
	IloNumVar[] Eta;
	double[] ZhValue;
	double[][] YhkValue;
	// 强对偶部分
	IloNumVar[][][] Whjk;
	IloNumVar[][] Ehk;
	// 上层问题约束
	IloRange[] upper1_2;
	IloRange[][] upper1_3;
	IloRange[][] upper1_4;
	IloRange[] upper1_8;
	IloRange[] upper1_9;
	IloRange[][] upper1_10;
	IloRange[] upper1_11;
	// 下层问题约束
	IloRange[][] lower2_2;
	IloRange[] lower2_3;
	IloRange[] lower2_6;
	IloRange[][] lower2_8To2_10;
	IloRange[] lower2_11;
	IloRange[][] lower2_14;
	IloRange[][] lower2_15;
	IloRange[] lower2_16;
	IloRange[] lower2_17;
	IloRange[][] lower2_18To2_21;
	IloRange[][] lower2_22;
	IloRange[] lower2_23;
	IloRange[] lower2_24;
	// 强对偶部分
	IloRange[][][] lower3_2;
	IloRange[][][] lower3_3;
	IloRange[][][] lower3_4;
	IloRange[][] lower3_6;
	IloRange[][] lower3_7;
	IloRange[][] lower3_8;

	long startTime; // 算法运行开始时间
	long endTime; // 算法运行结束时间
	private static final double M = 3000;
	private static final double ACCURANCY = 0.1; // C&CG算法终止的精度
	private static final double FLOATACCURACY = 0.001;
	double terminalSpan = 1800; // 程序最大运行时间，以s为单位
	int Iteration = 1;
	int IterationBenders = 0; // benders的迭代次数
	int IterationCcg = 0; // C&CG的迭代次数
	int jValue; // 贝位数减去安全距离再减1
	int changeThreshold = 4; // 同一Senario后将相应的Benders换成C&CG

	// BD与C&CG结合的部分
	ArrayList<UncertainSenario> uncertainSenario = new ArrayList<UncertainSenario>();

	// variableOutput表示是否在Console输出相关的变量取值
	// iterationDetail表示是否保存迭代过程中的模型求解信息
	public void upperModel(Input input, boolean variableOutput, boolean iterationDetail) throws IOException {
		try {
			startTime = System.currentTimeMillis();
			upperPro = new IloCplex();
			// 变量
			Xhk = new IloNumVar[input.numBay][input.numQC];
			CmaxUpper = new IloNumVar[1];
			alpha = new IloNumVar[input.numQC];
			beta = new IloNumVar[input.numQC];
			// 约束
			upper1_2 = new IloRange[input.numBay];
			upper1_3 = new IloRange[input.numBay][input.numQC];
			upper1_4 = new IloRange[input.numBay][input.numQC];
			upper1_8 = new IloRange[input.numQC];
			upper1_9 = new IloRange[input.numQC];

			// 设置变量
			for (int i = 0; i < input.numBay; i++) {
				for (int j = 0; j < input.numQC; j++) {
					Xhk[i][j] = upperPro.numVar(0, 1, IloNumVarType.Int, "X" + i + j);
				}
			}

			CmaxUpper[0] = upperPro.numVar(0, Double.MAX_VALUE, IloNumVarType.Float, "CmaxUpper");

			for (int i = 0; i < input.numQC; i++) {
				alpha[i] = upperPro.numVar(1, input.numBay, IloNumVarType.Float, "alpha" + i);
			}

			for (int i = 0; i < input.numQC; i++) {
				beta[i] = upperPro.numVar(1, input.numBay, IloNumVarType.Float, "beta" + i);
			}

			// 目标函数
			IloNumExpr upperObj = upperPro.prod(CmaxUpper[0], 1);
			upperPro.addMinimize(upperObj);

			// 上层问题约束
			// 约束1-2
			for (int i = 0; i < input.numBay; i++) {
				if (input.processTime[i] != 0) {
					IloNumExpr expr = upperPro.numExpr();
					for (int j = 0; j < input.numQC; j++) {
						expr = upperPro.sum(expr, upperPro.prod(Xhk[i][j], 1));
					}
					upper1_2[i] = upperPro.addEq(expr, 1, "constraint1_2 " + i);
				}
			}

			// 约束1-3
			for (int i = 0; i < input.numBay; i++) {
				for (int j = 0; j < input.numQC; j++) {
					IloNumExpr expr = upperPro.numExpr();
					expr = upperPro.sum(expr, upperPro.prod(Xhk[i][j], M - i - 1));
					expr = upperPro.sum(expr, upperPro.prod(alpha[j], 1));
					upper1_3[i][j] = upperPro.addLe(expr, M, "constraint1_3 " + i + "," + j);
				}
			}

			// 约束1-4
			for (int i = 0; i < input.numBay; i++) {
				for (int j = 0; j < input.numQC; j++) {
					IloNumExpr expr = upperPro.numExpr();
					expr = upperPro.sum(expr, upperPro.prod(Xhk[i][j], -M - i - 1));
					expr = upperPro.sum(expr, upperPro.prod(beta[j], 1));
					upper1_4[i][j] = upperPro.addGe(expr, -M, "constraint1_4 " + i + "," + j);
				}
			}

			// 约束1-8
			for (int i = 0; i < input.numQC - 1; i++) {
				IloNumExpr expr = upperPro.numExpr();
				expr = upperPro.sum(expr, upperPro.prod(alpha[i], 1));
				expr = upperPro.sum(expr, upperPro.prod(alpha[i + 1], -1));
				upper1_8[i] = upperPro.addLe(expr, -input.safetyMargin - 1, "constraint1_8 " + i);
			}

			// 约束1-9
			for (int i = 0; i < input.numQC - 1; i++) {
				IloNumExpr expr = upperPro.numExpr();
				expr = upperPro.sum(expr, upperPro.prod(beta[i], 1));
				expr = upperPro.sum(expr, upperPro.prod(beta[i + 1], -1));
				upper1_9[i] = upperPro.addLe(expr, -input.safetyMargin - 1, "constraint1_9 " + i);
			}

			if (upperPro.solve()) {
				System.out.println("求解成功了");
				upperPro.exportModel("MP.lp");
				XhkValue = new double[input.numBay][input.numQC];
				alphaValue = new double[input.numQC];
				betaValue = new double[input.numQC];
				for (int i = 0; i < input.numBay; i++) {
					for (int j = 0; j < input.numQC; j++) {
						XhkValue[i][j] = upperPro.getValue(Xhk[i][j]);
					}
				}
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
		} catch (IloException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	// variableOutput表示是否在Console输出相关的变量取值
	// iterationDetail表示是否保存迭代过程中的模型求解信息
	public void lowerModelKKT(Input input, boolean variableOutput, boolean iterationDetail) throws IOException {
		try {
			lowerPro = new IloCplex();
			jValue = input.numBay - input.safetyMargin - 1;
			// 变量
			Yhk = new IloNumVar[input.numBay][input.numQC];
			CmaxLower = new IloNumVar[1];
			Zh = new IloNumVar[input.numBay];
			Niujk = new IloNumVar[jValue][input.numQC];
			Miuk = new IloNumVar[input.numQC];
			Gammajk = new IloNumVar[jValue][input.numQC];
			Epsilonk = new IloNumVar[input.numQC];
			Zetahk = new IloNumVar[input.numBay][input.numQC];
			Eta = new IloNumVar[1];
			// 约束
			lower2_2 = new IloRange[jValue][input.numQC];
			lower2_3 = new IloRange[input.numQC];
			lower2_6 = new IloRange[1];
			lower2_8To2_10 = new IloRange[input.numBay][input.numQC];
			lower2_11 = new IloRange[1];
			lower2_14 = new IloRange[jValue][input.numQC - 1];
			lower2_15 = new IloRange[jValue][input.numQC - 1];
			lower2_16 = new IloRange[input.numQC];
			lower2_17 = new IloRange[input.numQC];
			lower2_18To2_21 = new IloRange[input.numBay][input.numQC];
			lower2_22 = new IloRange[input.numBay][input.numQC];
			lower2_23 = new IloRange[1];
			lower2_24 = new IloRange[1];

			// 设置变量
			for (int i = 0; i < input.numBay; i++) {
				for (int j = 0; j < input.numQC; j++) {
					Yhk[i][j] = lowerPro.numVar(0, Double.MAX_VALUE, IloNumVarType.Float, "Yhk_" + i + "," + j);
				}
			}

			CmaxLower[0] = lowerPro.numVar(0, Double.MAX_VALUE, IloNumVarType.Float, "CmaxLower");

			for (int i = 0; i < input.numBay; i++) {
				Zh[i] = lowerPro.numVar(0, 1, IloNumVarType.Float, "Zh_" + i);
			}

			for (int i = 0; i < jValue; i++) {
				for (int j = 0; j < input.numQC - 1; j++) {
					Niujk[i][j] = lowerPro.numVar(0, Double.MAX_VALUE, IloNumVarType.Float, "Niujk_" + i + "," + j);
				}
			}

			for (int i = 0; i < input.numQC; i++) {
				Miuk[i] = lowerPro.numVar(0, Double.MAX_VALUE, IloNumVarType.Float, "Miuk_" + i);
			}

			for (int i = 0; i < jValue; i++) {
				for (int j = 0; j < input.numQC; j++) {
					Gammajk[i][j] = lowerPro.numVar(0, 1, IloNumVarType.Int, "Gammajk_" + i + "," + j);
				}
			}

			for (int i = 0; i < input.numQC; i++) {
				Epsilonk[i] = lowerPro.numVar(0, 1, IloNumVarType.Int, "Epsilonk_" + i);
			}

			for (int i = 0; i < input.numBay; i++) {
				for (int j = 0; j < input.numQC; j++) {
					Zetahk[i][j] = lowerPro.numVar(0, 1, IloNumVarType.Int, "Zetahk_" + i + "," + j);
				}
			}
			Eta[0] = lowerPro.numVar(0, 1, IloNumVarType.Int, "Eta");

			// 设置约束
			// 2-2
			for (int i = 0; i < jValue; i++) {
				for (int j = 0; j < input.numQC - 1; j++) {
					IloNumExpr expr = lowerPro.numExpr();
					double rhs = 0;
					for (int k = 0; k <= i; k++) {
						rhs += -XhkValue[k][j] * input.processTime[k];
						expr = lowerPro.sum(expr,
								lowerPro.prod(Zh[k], input.uncertain * input.processTime[k] * XhkValue[k][j]));
						expr = lowerPro.sum(expr, lowerPro.prod(Yhk[k][j], 1));
					}
					for (int k = input.safetyMargin + 1; k < i + input.safetyMargin + 2; k++) {
						rhs += XhkValue[k][j + 1] * input.processTime[k];
						expr = lowerPro.sum(expr,
								lowerPro.prod(Zh[k], -input.uncertain * input.processTime[k] * XhkValue[k][j + 1]));
						expr = lowerPro.sum(expr, lowerPro.prod(Yhk[k][j + 1], -1));
					}
					rhs += (alphaValue[j] - alphaValue[j + 1] + input.safetyMargin + 1) * input.traverseTime;
					lower2_2[i][j] = lowerPro.addGe(expr, rhs, "Constraint2_2 " + i + "," + j);
				}
			}

			// 2-3
			for (int i = 0; i < input.numQC; i++) {
				IloNumExpr expr = lowerPro.numExpr();
				double rhs = 0;
				for (int j = 0; j < input.numBay; j++) {
					rhs += -input.processTime[j] * XhkValue[j][i];
					expr = lowerPro.sum(expr,
							lowerPro.prod(Zh[j], input.uncertain * input.processTime[j] * XhkValue[j][i]));
					expr = lowerPro.sum(expr, lowerPro.prod(Yhk[j][i], 1));
				}
				rhs += -(betaValue[i] - alphaValue[i]) * input.traverseTime;
				expr = lowerPro.sum(expr, lowerPro.prod(CmaxLower[0], -1));
				lower2_3[i] = lowerPro.addLe(expr, rhs, "Constraint2_3 " + i);
			}

			// 2-6
			IloNumExpr expr = lowerPro.numExpr();
			for (int i = 0; i < input.numBay; i++) {
				expr = lowerPro.sum(expr, lowerPro.prod(Zh[i], 1));
			}
			lower2_6[0] = lowerPro.addLe(expr, input.budget, "Constraint2_6");

			// 2-8
			for (int i = input.safetyMargin + 1; i < input.numBay - input.safetyMargin - 1; i++) {
				for (int j = 1; j < input.numQC - 1; j++) {
					expr = lowerPro.numExpr();
					for (int k = i; k < input.numBay - input.safetyMargin - 1; k++) {
						expr = lowerPro.sum(expr, lowerPro.prod(Niujk[k][j], 1));
					}
					for (int j2 = i - input.safetyMargin - 1; j2 < input.numBay - input.safetyMargin - 1; j2++) {
						expr = lowerPro.sum(expr, lowerPro.prod(Niujk[j2][j - 1], -1));
					}
					expr = lowerPro.sum(expr, lowerPro.prod(Miuk[j], -1));
					lower2_8To2_10[i][j] = lowerPro.addLe(expr, 0, "Constraint2_8 " + i + "," + j);
				}
			}

			// 2-9
			// part1
			for (int i = input.safetyMargin + 1; i < input.numBay - input.safetyMargin - 1; i++) {
				expr = lowerPro.numExpr();
				for (int j = i; j < input.numBay - input.safetyMargin - 1; j++) {
					expr = lowerPro.sum(expr, lowerPro.prod(Niujk[j][0], 1));
				}
				expr = lowerPro.sum(expr, lowerPro.prod(Miuk[0], -1));
				lower2_8To2_10[i][0] = lowerPro.addLe(expr, 0, "Constraint2_9 " + i + ",0");
			}

			// part2
			for (int i = 0; i < input.safetyMargin + 1; i++) {
				for (int j = 0; j < input.numQC - 1; j++) {
					expr = lowerPro.numExpr();
					for (int k = i; k < input.numBay - input.safetyMargin - 1; k++) {
						expr = lowerPro.sum(expr, lowerPro.prod(Niujk[k][j], 1));
					}
					expr = lowerPro.sum(expr, lowerPro.prod(Miuk[j], -1));
					lower2_8To2_10[i][j] = lowerPro.addLe(expr, 0, "Constraint2_9 " + i + "," + j);
				}
			}

			// 2-10
			// part1
			for (int i = input.numBay - input.safetyMargin - 1; i < input.numBay; i++) {
				for (int j = 1; j < input.numQC; j++) {
					expr = lowerPro.numExpr();
					for (int j2 = i - input.safetyMargin - 1; j2 < input.numBay - input.safetyMargin - 1; j2++) {
						expr = lowerPro.sum(expr, lowerPro.prod(Niujk[j2][j - 1], -1));
					}
					expr = lowerPro.sum(expr, lowerPro.prod(Miuk[j], -1));
					lower2_8To2_10[i][j] = lowerPro.addLe(expr, 0, "Constraint2_10 " + i + "," + j);
				}
			}

			// part2
			for (int i = input.safetyMargin + 1; i < input.numBay - input.safetyMargin - 1; i++) {
				expr = lowerPro.numExpr();
				for (int j = i - input.safetyMargin - 1; j < input.numBay - input.safetyMargin - 1; j++) {
					expr = lowerPro.sum(expr, lowerPro.prod(Niujk[j][input.numQC - 2], -1));
				}
				expr = lowerPro.sum(expr, lowerPro.prod(Miuk[input.numQC - 1], -1));
				int temp = input.numQC - 1;
				lower2_8To2_10[i][input.numQC - 1] = lowerPro.addLe(expr, 0, "Constraint2_10 " + i + "," + temp);
			}

			// 2-11
			expr = lowerPro.numExpr();
			for (int i = 0; i < input.numQC; i++) {
				expr = lowerPro.sum(expr, lowerPro.prod(Miuk[i], 1));
			}
			lower2_11[0] = lowerPro.addLe(expr, 1, "Constraint2_11 ");

			// 2-14
			for (int i = 0; i < jValue; i++) {
				for (int j = 0; j < input.numQC - 1; j++) {
					expr = lowerPro.numExpr();
					double rhs = 0;
					for (int k = 0; k <= i; k++) {
						rhs += -XhkValue[k][j] * input.processTime[k];
						expr = lowerPro.sum(expr,
								lowerPro.prod(Zh[k], input.uncertain * input.processTime[k] * XhkValue[k][j]));
						expr = lowerPro.sum(expr, lowerPro.prod(Yhk[k][j], 1));
					}
					for (int k = input.safetyMargin + 1; k <= i + input.safetyMargin + 1; k++) {
						rhs += XhkValue[k][j + 1] * input.processTime[k];
						expr = lowerPro.sum(expr,
								lowerPro.prod(Zh[k], -input.uncertain * input.processTime[k] * XhkValue[k][j + 1]));
						expr = lowerPro.sum(expr, lowerPro.prod(Yhk[k][j + 1], -1));
					}
					rhs += (alphaValue[j] - alphaValue[j + 1] + input.safetyMargin + 1) * input.traverseTime;
					expr = lowerPro.sum(expr, lowerPro.prod(Gammajk[i][j], -M));
					lower2_14[i][j] = lowerPro.addLe(expr, rhs, "Constraint2_14 " + i + "," + j);
				}
			}

			// 2-15
			for (int i = 0; i < jValue; i++) {
				for (int j = 0; j < input.numQC - 1; j++) {
					expr = lowerPro.numExpr();
					expr = lowerPro.sum(expr, lowerPro.prod(Niujk[i][j], 1));
					expr = lowerPro.sum(expr, lowerPro.prod(Gammajk[i][j], M));
					lower2_15[i][j] = lowerPro.addLe(expr, M, "Constraint2_15 " + i + "," + j);
				}
			}

			// 2-16
			for (int i = 0; i < input.numQC; i++) {
				expr = lowerPro.numExpr();
				double rhs = 0;
				for (int j = 0; j < input.numBay; j++) {
					rhs += input.processTime[j] * XhkValue[j][i];
					expr = lowerPro.sum(expr,
							lowerPro.prod(Zh[j], -input.uncertain * input.processTime[j] * XhkValue[j][i]));
					expr = lowerPro.sum(expr, lowerPro.prod(Yhk[j][i], -1));
				}
				rhs += (betaValue[i] - alphaValue[i]) * input.traverseTime;
				expr = lowerPro.sum(expr, lowerPro.prod(CmaxLower[0], 1));
				expr = lowerPro.sum(expr, lowerPro.prod(Epsilonk[i], -M));
				lower2_16[i] = lowerPro.addLe(expr, rhs, "Constraint2_16 " + i);
			}

			// 2-17
			for (int i = 0; i < input.numQC; i++) {
				expr = lowerPro.numExpr();
				expr = lowerPro.sum(expr, lowerPro.prod(Miuk[i], 1));
				expr = lowerPro.sum(expr, lowerPro.prod(Epsilonk[i], M));
				lower2_17[i] = lowerPro.addLe(expr, M, "Constraint2_17 " + i);
			}

			// 2-18
			for (int i = input.safetyMargin + 1; i < input.numBay - input.safetyMargin - 1; i++) {
				for (int j = 1; j < input.numQC - 1; j++) {
					expr = lowerPro.numExpr();
					for (int k = i; k < input.numBay - input.safetyMargin - 1; k++) {
						expr = lowerPro.sum(expr, lowerPro.prod(Niujk[k][j], -1));
					}
					for (int j2 = i - input.safetyMargin - 1; j2 < input.numBay - input.safetyMargin - 1; j2++) {
						expr = lowerPro.sum(expr, lowerPro.prod(Niujk[j2][j - 1], 1));
					}
					expr = lowerPro.sum(expr, lowerPro.prod(Miuk[j], 1));
					expr = lowerPro.sum(expr, lowerPro.prod(Zetahk[i][j], -M));
					lower2_18To2_21[i][j] = lowerPro.addLe(expr, 0, "Constraint2_18 " + i + "," + j);
				}
			}

			// 2-19
			// part1
			for (int i = input.safetyMargin + 1; i < input.numBay - input.safetyMargin - 1; i++) {
				expr = lowerPro.numExpr();
				for (int j = i; j < input.numBay - input.safetyMargin - 1; j++) {
					expr = lowerPro.sum(expr, lowerPro.prod(Niujk[j][0], -1));
				}
				expr = lowerPro.sum(expr, lowerPro.prod(Miuk[0], 1));
				expr = lowerPro.sum(expr, lowerPro.prod(Zetahk[i][0], -M));
				lower2_18To2_21[i][0] = lowerPro.addLe(expr, 0, "Constraint2_19 " + i + ",0");
			}

			// part2
			for (int i = 0; i < input.safetyMargin + 1; i++) {
				for (int j = 0; j < input.numQC - 1; j++) {
					expr = lowerPro.numExpr();
					for (int k = i; k < input.numBay - input.safetyMargin - 1; k++) {
						expr = lowerPro.sum(expr, lowerPro.prod(Niujk[k][j], -1));
					}
					expr = lowerPro.sum(expr, lowerPro.prod(Miuk[j], 1));
					expr = lowerPro.sum(expr, lowerPro.prod(Zetahk[i][j], -M));
					lower2_18To2_21[i][j] = lowerPro.addLe(expr, 0, "Constraint2_19 " + i + "," + j);
				}
			}

			// 2-20
			// part1
			for (int i = input.numBay - input.safetyMargin - 1; i < input.numBay; i++) {
				for (int j = 1; j < input.numQC; j++) {
					expr = lowerPro.numExpr();
					for (int j2 = i - input.safetyMargin - 1; j2 < input.numBay - input.safetyMargin - 1; j2++) {
						expr = lowerPro.sum(expr, lowerPro.prod(Niujk[j2][j - 1], 1));
					}
					expr = lowerPro.sum(expr, lowerPro.prod(Miuk[j], 1));
					expr = lowerPro.sum(expr, lowerPro.prod(Zetahk[i][j], -M));
					lower2_18To2_21[i][j] = lowerPro.addLe(expr, 0, "Constraint2_20 " + i + "," + j);
				}
			}

			// part2
			for (int i = input.safetyMargin + 1; i < input.numBay - input.safetyMargin - 1; i++) {
				expr = lowerPro.numExpr();
				for (int j = i - input.safetyMargin - 1; j < input.numBay - input.safetyMargin - 1; j++) {
					expr = lowerPro.sum(expr, lowerPro.prod(Niujk[j][input.numQC - 2], 1));
				}
				expr = lowerPro.sum(expr, lowerPro.prod(Miuk[input.numQC - 1], 1));
				expr = lowerPro.sum(expr, lowerPro.prod(Zetahk[i][input.numQC - 1], -M));
				int temp = input.numQC - 1;
				lower2_18To2_21[i][input.numQC - 1] = lowerPro.addLe(expr, 0, "Constraint2_20 " + i + "," + temp);
			}

			// 2-21
			// part1
			for (int i = 0; i < input.safetyMargin + 1; i++) {
				expr = lowerPro.numExpr();
				expr = lowerPro.sum(expr, lowerPro.prod(Miuk[input.numQC - 1], 1));
				expr = lowerPro.sum(expr, lowerPro.prod(Zetahk[i][input.numQC - 1], -M));
				lower2_18To2_21[i][input.numQC - 1] = lowerPro.addLe(expr, 0, "Constraint2_21 " + i);
			}

			// part2
			for (int i = input.numBay - input.safetyMargin - 1; i < input.numBay; i++) {
				expr = lowerPro.numExpr();
				expr = lowerPro.sum(expr, lowerPro.prod(Miuk[0], 1));
				expr = lowerPro.sum(expr, lowerPro.prod(Zetahk[i][0], -M));
				lower2_18To2_21[i][0] = lowerPro.addLe(expr, 0, "Constraint2_21 " + i);
			}

			// 2-22
			for (int i = 0; i < input.numBay; i++) {
				for (int j = 0; j < input.numQC; j++) {
					expr = lowerPro.numExpr();
					expr = lowerPro.sum(expr, lowerPro.prod(Yhk[i][j], 1));
					expr = lowerPro.sum(expr, lowerPro.prod(Zetahk[i][j], M));
					lower2_22[i][j] = lowerPro.addLe(expr, M, "Constraint2_22 " + i + "," + j);
				}
			}

			// 2-23
			expr = lowerPro.numExpr();
			for (int i = 0; i < input.numQC; i++) {
				expr = lowerPro.sum(expr, lowerPro.prod(Miuk[i], -1));
			}
			expr = lowerPro.sum(expr, lowerPro.prod(Eta[0], -M));
			lower2_23[0] = lowerPro.addLe(expr, -1, "Constraint2_23");

			// 2-24
			expr = lowerPro.numExpr();
			expr = lowerPro.sum(expr, lowerPro.prod(CmaxLower[0], 1));
			expr = lowerPro.sum(expr, lowerPro.prod(Eta[0], M));
			lower2_24[0] = lowerPro.addLe(expr, M, "Constraint2_24");

			// 设置目标函数
			// 原问题目标函数
			IloNumExpr obj = lowerPro.numExpr();
			obj = lowerPro.sum(obj, CmaxLower[0]);
			lowerProObj = lowerPro.addMaximize(obj);
//			lowerProObj = lowerPro.addMinimize(obj);
////		 对偶问题目标函数值
//			IloNumExpr dualObj = lowerPro.numExpr();
//			for (int i = 0; i < input.numBay - input.safetyMargin - 1; i++) {
//				for (int j = 0; j < input.numQC - 1; j++) {
//					double coefficient = 0;
//					for (int k = 0; k < i + 1; k++) {
//						coefficient += -input.processTime[k] * XhkValue[k][j];
//					}
//					for (int k = input.safetyMargin + 1; k < i + input.safetyMargin + 1; k++) {
//						coefficient += input.processTime[k] * XhkValue[k][j + 1];
//					}
//					coefficient += (input.safetyMargin + 1 + alphaValue[j] - alphaValue[j + 1]) * input.traverseTime;
//					dualObj = lowerPro.sum(dualObj, lowerPro.prod(Niujk[i][j], coefficient));
//				}
//			}
//			for (int i = 0; i < input.numQC; i++) {
//				double coefficient = 0;
//				for (int j = 0; j < input.numBay; j++) {
//					coefficient += input.processTime[j] * XhkValue[j][i];
//				}
//				coefficient += (betaValue[i] -alphaValue[i]) * input.traverseTime;
//				dualObj = lowerPro.sum(dualObj, lowerPro.prod(Miuk[i], coefficient));
//			}
//			lowerProObj = lowerPro.addMaximize(dualObj);
			lowerPro.exportModel("SP.lp");

			if (lowerPro.solve()) {
				System.out.println("Solution status = " + lowerPro.getStatus());
				System.out.println("Solution value  = " + lowerPro.getObjValue());

				if (variableOutput) {
					for (int i = 0; i < input.numBay; i++) {
						for (int j = 0; j < input.numQC; j++) {
							System.out.println("Y" + i + "," + j + ": Value = " + lowerPro.getValue(Yhk[i][j]));
						}
					}
					for (int i = 0; i < input.numBay; i++) {
						System.out.println("Z" + i + ": Value = " + lowerPro.getValue(Zh[i]));
					}
				}

				input.writerOverall
						.write("Iteration " + Iteration + " lowerObj(KKT):" + lowerPro.getObjValue() + "\r\n");
				input.writerOverall.write("\r\n");

				if (iterationDetail) {
					lowerPro.writeSolution("./result of benchmark/test_" + input.testID + "/test_" + input.testID
							+ " iteration_" + Iteration + "lowerKKT.txt");
				}

			} else {
				System.out.println("下层问题求解失败");
				System.exit(0);
			}
		} catch (IloException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

	// iteration为true时表示进行的是迭代运算，为false时表示进行的是deterministic时的attack
	// variableOutput表示是否在Console输出相关的变量取值
	// iterationDetail表示是否保存迭代过程中的模型求解信息
	public void lowerModelSD(Input input, boolean iteration, boolean variableOutput, boolean iterationDetail)
			throws IOException {
		try {
			lowerPro = new IloCplex();
			jValue = input.numBay - input.safetyMargin - 1;
			// 变量
			CmaxLower = new IloNumVar[1];
			Zh = new IloNumVar[input.numBay];
			Niujk = new IloNumVar[jValue][input.numQC - 1];
			Miuk = new IloNumVar[input.numQC];
			Whjk = new IloNumVar[input.numBay][jValue][input.numQC - 1];
			Ehk = new IloNumVar[input.numBay][input.numQC];
			// 约束
			lower2_6 = new IloRange[1];
			lower2_8To2_10 = new IloRange[input.numBay][input.numQC];
			lower2_11 = new IloRange[1];
			lower3_2 = new IloRange[input.numBay][input.numBay - input.safetyMargin - 1][input.numQC - 1];
			lower3_3 = new IloRange[input.numBay][input.numBay - input.safetyMargin - 1][input.numQC - 1];
			lower3_4 = new IloRange[input.numBay][input.numBay - input.safetyMargin - 1][input.numQC - 1];
			lower3_6 = new IloRange[input.numBay][input.numQC];
			lower3_7 = new IloRange[input.numBay][input.numQC];
			lower3_8 = new IloRange[input.numBay][input.numQC];
			// 设置变量
			CmaxLower[0] = lowerPro.numVar(0, Double.MAX_VALUE, IloNumVarType.Float, "CmaxLower");
			for (int i = 0; i < input.numBay; i++) {
				Zh[i] = lowerPro.numVar(0, 1, IloNumVarType.Int, "Zh_" + i);
			}
			for (int i = 0; i < jValue; i++) {
				for (int j = 0; j < input.numQC - 1; j++) {
					Niujk[i][j] = lowerPro.numVar(0, Double.MAX_VALUE, IloNumVarType.Float, "Niujk_" + i + "," + j);
				}
			}
			for (int i = 0; i < input.numQC; i++) {
				Miuk[i] = lowerPro.numVar(0, Double.MAX_VALUE, IloNumVarType.Float, "Miuk_" + i);
			}
			for (int i = 0; i < input.numBay; i++) {
				for (int j = 0; j < jValue; j++) {
					for (int j2 = 0; j2 < input.numQC - 1; j2++) {
						Whjk[i][j][j2] = lowerPro.numVar(0, Double.MAX_VALUE, IloNumVarType.Float,
								"Whjk_" + i + "," + j + "," + j2);
					}
				}
			}
			for (int i = 0; i < input.numBay; i++) {
				for (int j = 0; j < input.numQC; j++) {
					Ehk[i][j] = lowerPro.numVar(0, Double.MAX_VALUE, IloNumVarType.Float, "Ehk_" + i + "," + j);
				}
			}

			// 设置约束
			// 2-6
			IloNumExpr expr = lowerPro.numExpr();
			for (int i = 0; i < input.numBay; i++) {
				expr = lowerPro.sum(expr, lowerPro.prod(Zh[i], 1));
			}
			lower2_6[0] = lowerPro.addLe(expr, (int) input.budget, "Constraint2_6");
			// 2-8
			for (int i = input.safetyMargin + 1; i < input.numBay - input.safetyMargin - 1; i++) {
				for (int j = 1; j < input.numQC - 1; j++) {
					expr = lowerPro.numExpr();
					for (int k = i; k < input.numBay - input.safetyMargin - 1; k++) {
						expr = lowerPro.sum(expr, lowerPro.prod(Niujk[k][j], 1));
					}
					for (int j2 = i - input.safetyMargin - 1; j2 < input.numBay - input.safetyMargin - 1; j2++) {
						expr = lowerPro.sum(expr, lowerPro.prod(Niujk[j2][j - 1], -1));
					}
					expr = lowerPro.sum(expr, lowerPro.prod(Miuk[j], -1));
					lower2_8To2_10[i][j] = lowerPro.addLe(expr, 0, "Constraint2_8 " + i + "," + j);
				}
			}
			// 2-9
			// part1
			for (int i = input.safetyMargin + 1; i < input.numBay - input.safetyMargin - 1; i++) {
				expr = lowerPro.numExpr();
				for (int j = i; j < input.numBay - input.safetyMargin - 1; j++) {
					expr = lowerPro.sum(expr, lowerPro.prod(Niujk[j][0], 1));
				}
				expr = lowerPro.sum(expr, lowerPro.prod(Miuk[0], -1));
				lower2_8To2_10[i][0] = lowerPro.addLe(expr, 0, "Constraint2_9 " + i + ",0");
			}
			// part2
			for (int i = 0; i < input.safetyMargin + 1; i++) {
				for (int j = 0; j < input.numQC - 1; j++) {
					expr = lowerPro.numExpr();
					for (int k = i; k < input.numBay - input.safetyMargin - 1; k++) {
						expr = lowerPro.sum(expr, lowerPro.prod(Niujk[k][j], 1));
					}
					expr = lowerPro.sum(expr, lowerPro.prod(Miuk[j], -1));
					lower2_8To2_10[i][j] = lowerPro.addLe(expr, 0, "Constraint2_9 " + i + "," + j);
				}
			}
			// 2-10
			// part1
			for (int i = input.numBay - input.safetyMargin - 1; i < input.numBay; i++) {
				for (int j = 1; j < input.numQC; j++) {
					expr = lowerPro.numExpr();
					for (int j2 = i - input.safetyMargin - 1; j2 < input.numBay - input.safetyMargin - 1; j2++) {
						expr = lowerPro.sum(expr, lowerPro.prod(Niujk[j2][j - 1], -1));
					}
					expr = lowerPro.sum(expr, lowerPro.prod(Miuk[j], -1));
					lower2_8To2_10[i][j] = lowerPro.addLe(expr, 0, "Constraint2_10 " + i + "," + j);
				}
			}
			// part2
			for (int i = input.safetyMargin + 1; i < input.numBay - input.safetyMargin - 1; i++) {
				expr = lowerPro.numExpr();
				for (int j = i - input.safetyMargin - 1; j < input.numBay - input.safetyMargin - 1; j++) {
					expr = lowerPro.sum(expr, lowerPro.prod(Niujk[j][input.numQC - 2], -1));
				}
				expr = lowerPro.sum(expr, lowerPro.prod(Miuk[input.numQC - 1], -1));
				int temp = input.numQC - 1;
				lower2_8To2_10[i][input.numQC - 1] = lowerPro.addLe(expr, 0, "Constraint2_10 " + i + "," + temp);
			}
			// 2-11
			expr = lowerPro.numExpr();
			for (int i = 0; i < input.numQC; i++) {
				expr = lowerPro.sum(expr, lowerPro.prod(Miuk[i], 1));
			}
			lower2_11[0] = lowerPro.addLe(expr, 1, "Constraint2_11 ");
			// 3-2
			for (int i = 0; i < input.numBay; i++) {
				for (int j = 0; j < jValue; j++) {
					for (int j2 = 0; j2 < input.numQC - 1; j2++) {
						expr = lowerPro.numExpr();
						expr = lowerPro.sum(expr, lowerPro.prod(Whjk[i][j][j2], 1));
						expr = lowerPro.sum(expr, lowerPro.prod(Niujk[j][j2], -1));
						lower3_2[i][j][j2] = lowerPro.addLe(expr, 0, "Constraint3_2 " + i + "," + j + "," + j2);
					}
				}
			}
			// 3-3
			for (int i = 0; i < input.numBay; i++) {
				for (int j = 0; j < jValue; j++) {
					for (int j2 = 0; j2 < input.numQC - 1; j2++) {
						expr = lowerPro.numExpr();
						expr = lowerPro.sum(expr, lowerPro.prod(Whjk[i][j][j2], 1));
						expr = lowerPro.sum(expr, lowerPro.prod(Zh[i], -M));
						lower3_3[i][j][j2] = lowerPro.addLe(expr, 0, "Constraint3_3 " + i + "," + j + "," + j2);
					}
				}
			}
			// 3-4
			for (int i = 0; i < input.numBay; i++) {
				for (int j = 0; j < jValue; j++) {
					for (int j2 = 0; j2 < input.numQC - 1; j2++) {
						expr = lowerPro.numExpr();
						expr = lowerPro.sum(expr, lowerPro.prod(Whjk[i][j][j2], 1));
						expr = lowerPro.sum(expr, lowerPro.prod(Niujk[j][j2], -1));
						expr = lowerPro.sum(expr, lowerPro.prod(Zh[i], -M));
						lower3_4[i][j][j2] = lowerPro.addGe(expr, -M, "Constraint3_4 " + i + "," + j + "," + j2);
					}
				}
			}
			// 3-6
			for (int i = 0; i < input.numBay; i++) {
				for (int j = 0; j < input.numQC; j++) {
					expr = lowerPro.numExpr();
					expr = lowerPro.sum(expr, lowerPro.prod(Ehk[i][j], 1));
					expr = lowerPro.sum(expr, lowerPro.prod(Miuk[j], -1));
					lower3_6[i][j] = lowerPro.addLe(expr, 0, "Constraint3_6 " + i + "," + j);
				}
			}
			// 3-7
			for (int i = 0; i < input.numBay; i++) {
				for (int j = 0; j < input.numQC; j++) {
					expr = lowerPro.numExpr();
					expr = lowerPro.sum(expr, lowerPro.prod(Ehk[i][j], 1));
					expr = lowerPro.sum(expr, lowerPro.prod(Zh[i], -M));
					lower3_7[i][j] = lowerPro.addLe(expr, 0, "Constraint3_7 " + i + "," + j);
				}
			}
			// 3-8
			for (int i = 0; i < input.numBay; i++) {
				for (int j = 0; j < input.numQC; j++) {
					expr = lowerPro.numExpr();
					expr = lowerPro.sum(expr, lowerPro.prod(Ehk[i][j], 1));
					expr = lowerPro.sum(expr, lowerPro.prod(Miuk[j], -1));
					expr = lowerPro.sum(expr, lowerPro.prod(Zh[i], -M));
					lower3_8[i][j] = lowerPro.addGe(expr, -M, "Constraint3_8 " + i + "," + j);
				}
			}
			// 设置目标函数
			IloNumExpr dualObj = lowerPro.numExpr();
			for (int i = 0; i < input.numBay - input.safetyMargin - 1; i++) {
				for (int j = 0; j < input.numQC - 1; j++) {
					for (int j2 = 0; j2 <= i; j2++) {
						dualObj = lowerPro.sum(dualObj,
								lowerPro.prod(Niujk[i][j], -input.processTime[j2] * XhkValue[j2][j]));
						dualObj = lowerPro.sum(dualObj, lowerPro.prod(Whjk[j2][i][j],
								-input.processTime[j2] * XhkValue[j2][j] * input.uncertain));
					}
					for (int j2 = input.safetyMargin + 1; j2 <= i + input.safetyMargin + 1; j2++) {
						dualObj = lowerPro.sum(dualObj,
								lowerPro.prod(Niujk[i][j], input.processTime[j2] * XhkValue[j2][j + 1]));
						dualObj = lowerPro.sum(dualObj, lowerPro.prod(Whjk[j2][i][j],
								input.processTime[j2] * XhkValue[j2][j + 1] * input.uncertain));
					}
					dualObj = lowerPro.sum(dualObj, lowerPro.prod(Niujk[i][j],
							(input.safetyMargin + 1 + alphaValue[j] - alphaValue[j + 1]) * input.traverseTime));
				}
			}

			for (int i = 0; i < input.numQC; i++) {
				for (int j = 0; j < input.numBay; j++) {
					dualObj = lowerPro.sum(dualObj, lowerPro.prod(Miuk[i], input.processTime[j] * XhkValue[j][i]));
					dualObj = lowerPro.sum(dualObj,
							lowerPro.prod(Ehk[j][i], input.processTime[j] * XhkValue[j][i] * input.uncertain));
				}
				dualObj = lowerPro.sum(dualObj,
						lowerPro.prod(Miuk[i], (betaValue[i] - alphaValue[i]) * input.traverseTime));
			}
			lowerProObj = lowerPro.addMaximize(dualObj);
//			lowerPro.exportModel("SP.lp");
			if (lowerPro.solve()) {
				System.out.println("Solution status = " + lowerPro.getStatus());
				System.out.println("Solution value  = " + lowerPro.getObjValue());
				ZhValue = lowerPro.getValues(Zh);

				if (variableOutput) {
					for (int i = 0; i < input.numBay; i++) {
						System.out.println("Z" + i + ": Value = " + ZhValue[i]);
					}
				}

				if (iteration) {
					input.writerOverall
							.write("Iteration " + Iteration + " lowerObj(SD):" + lowerPro.getObjValue() + "\r\n");
					input.writerOverall.write("\r\n");

					File file = new File("./result of benchmark/test_" + input.testID);
					if (!file.exists()) {// 如果文件夹不存在
						file.mkdir();// 创建文件夹
					}

					if (iterationDetail) {
						lowerPro.writeSolution("./result of benchmark/test_" + input.testID + "/test_" + input.testID
								+ " iteration_" + Iteration + "lowerSD.txt");
					}
				}
			} else {
				System.out.println("下层问题求解失败");
				System.exit(0);
			}
		} catch (IloException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

	// method选择KKT或者Strong
	// Duality，variableOutput为True时在console输出变量取值，iterationDetail为True时储存每一步迭代过程中结果的详细信息
	public void ccgIteration(Input input, String method, boolean variableOutput, boolean iterationDetail)
			throws IOException {
		try {
			double lowerBound = upperPro.getObjValue();
			double upperBound = lowerPro.getObjValue();
			double bestupperBound = upperBound; // 最好的上界

			while (upperBound - lowerBound >= ACCURANCY && (endTime - startTime) / 1000 < terminalSpan) {
				Iteration += 1;
				System.out.println("！！！第" + Iteration + "次加入cut！！！");

				// 创建变量
				IloNumVar[][] Yhkl = new IloNumVar[input.numBay][input.numQC];
				// 创建约束
				IloRange[][] upper1_10 = new IloRange[jValue][input.numQC - 1];
				IloRange[] upper1_11 = new IloRange[input.numQC];
				IloRange[] upper1_14 = new IloRange[1];

				// 设置变量
				for (int i = 0; i < input.numBay; i++) {
					for (int j = 0; j < input.numQC; j++) {
						Yhkl[i][j] = upperPro.numVar(0, Double.MAX_VALUE, IloNumVarType.Float, "Yhkl " + i + "," + j);
					}
				}

				// 设置约束
				// 1-10
				double[] ZhValue = new double[input.numBay];
				ZhValue = lowerPro.getValues(Zh);
				for (int i = 0; i < jValue; i++) {
					for (int j = 0; j < input.numQC - 1; j++) {
						IloNumExpr expr = upperPro.numExpr();
						for (int k = 0; k < i + 1; k++) {
							expr = upperPro.sum(expr, upperPro.prod(Xhk[k][j],
									input.processTime[k] + input.uncertain * input.processTime[k] * ZhValue[k]));
							expr = upperPro.sum(expr, upperPro.prod(Yhkl[k][j], 1));
						}
						expr = upperPro.sum(expr, upperPro.prod(alpha[j], -input.traverseTime));
						for (int k = input.safetyMargin + 1; k < i + input.safetyMargin + 2; k++) {
							expr = upperPro.sum(expr, upperPro.prod(Xhk[k][j + 1],
									-(input.processTime[k] + input.uncertain * input.processTime[k] * ZhValue[k])));
							expr = upperPro.sum(expr, upperPro.prod(Yhkl[k][j + 1], -1));
						}
						expr = upperPro.sum(expr, upperPro.prod(alpha[j + 1], input.traverseTime));
						upper1_10[i][j] = upperPro.addGe(expr, (input.safetyMargin + 1) * input.traverseTime,
								"Constraint1_10 " + i + "," + j);
					}
				}

				// 1-11
				for (int i = 0; i < input.numQC; i++) {
					IloNumExpr expr = upperPro.numExpr();
					for (int j = 0; j < input.numBay; j++) {
						expr = upperPro.sum(expr, upperPro.prod(Xhk[j][i],
								input.processTime[j] + input.uncertain * input.processTime[j] * ZhValue[j]));
						expr = upperPro.sum(expr, upperPro.prod(Yhkl[j][i], 1));
					}
					expr = upperPro.sum(expr, upperPro.prod(beta[i], input.traverseTime));
					expr = upperPro.sum(expr, upperPro.prod(alpha[i], -input.traverseTime));
					expr = upperPro.sum(expr, upperPro.prod(CmaxUpper[0], -1));
					upper1_11[i] = upperPro.addLe(expr, 0, "Constraint1_11 " + i);
				}

				// 将1-14以Callback的形式加入到模型中
				SpeciafiedValueAbortCallback checkIntegerSolution = new SpeciafiedValueAbortCallback(upperPro,
						lowerBound);
				upperPro.clearCallbacks();
				upperPro.use(checkIntegerSolution);

				// 根据设置的算法终止时间计算本次迭代的最长时间
				double leftTime = terminalSpan - (System.currentTimeMillis() - startTime) / 1000;
				upperPro.setParam(IloCplex.DoubleParam.TiLim, leftTime);

				// 对新的上层问题进行求解
				if (upperPro.solve()) {
					input.writerOverall
							.write("Iteration " + Iteration + " upperObj:" + upperPro.getObjValue() + "\r\n");

					if (iterationDetail) {
						upperPro.writeSolution("./result of benchmark/test_" + input.testID + "/test_" + input.testID
								+ " iteration_" + Iteration + "upper.txt");
					}

					System.out.println("上层问题第" + Iteration + "次求解成功了");
					lowerBound = upperPro.getValue(CmaxUpper[0]);
					System.out.println("下界 = " + lowerBound);

					// 更新X值
					for (int i = 0; i < XhkValue.length; i++) {
						for (int j = 0; j < XhkValue[0].length; j++) {
							XhkValue[i][j] = upperPro.getValue(Xhk[i][j]);
						}
					}

					// 跟新alpha和beta
					alphaValue = upperPro.getValues(alpha);
					betaValue = upperPro.getValues(beta);

					// 记录X的取值，也就是泊位和岸桥的分配关系，还有alpha和beta的取值，需要先对内容进行清空
					input.resetWriterAssignment();
					assignmentRecord(input);

					if (variableOutput) {
						for (int i = 0; i < input.numBay; i++) {
							for (int j = 0; j < input.numQC; j++) {
								System.out.println("X" + i + "," + j + ": Value = " + XhkValue[i][j]);
							}
						}

						for (int i = 0; i < input.numQC; i++) {
							System.out.println("alpha " + i + ": Value = " + alphaValue[i]);
						}

						for (int i = 0; i < input.numQC; i++) {
							System.out.println("beta " + i + ": Value = " + betaValue[i]);
						}

						for (int i = 0; i < input.numBay; i++) {
							for (int j = 0; j < input.numQC; j++) {
								System.out.println("Y" + i + "," + j + ": Value = " + upperPro.getValue(Yhkl[i][j]));
							}
						}
					}

				} else {
					System.out.println("上层问题未成功求解");
					System.exit(0);
				}

				if (method.equals("KKT")) {
					// 对子问题的部分约束进行更新
					// 2-2
					for (int i = 0; i < jValue; i++) {
						for (int j = 0; j < input.numQC - 1; j++) {
							IloNumExpr expr = lowerPro.numExpr();
							double rhs = 0;
							for (int k = 0; k <= i; k++) {
								rhs += -XhkValue[k][j] * input.processTime[k];
								expr = lowerPro.sum(expr,
										lowerPro.prod(Zh[k], input.uncertain * input.processTime[k] * XhkValue[k][j]));
								expr = lowerPro.sum(expr, lowerPro.prod(Yhk[k][j], 1));
							}

							for (int k = input.safetyMargin + 1; k < i + input.safetyMargin + 2; k++) {
								rhs += XhkValue[k][j + 1] * input.processTime[k];
								expr = lowerPro.sum(expr, lowerPro.prod(Zh[k],
										-input.uncertain * input.processTime[k] * XhkValue[k][j + 1]));
								expr = lowerPro.sum(expr, lowerPro.prod(Yhk[k][j + 1], -1));
							}
							rhs += (alphaValue[j] - alphaValue[j + 1] + input.safetyMargin + 1) * input.traverseTime;
							lowerPro.remove(lower2_2[i]);
							lower2_2[i][j] = lowerPro.addGe(expr, rhs, "Constraint2_2 " + i + "," + j);
						}
					}

					// 2-3
					for (int i = 0; i < input.numQC; i++) {
						IloNumExpr expr = lowerPro.numExpr();
						double rhs = 0;
						for (int j = 0; j < input.numBay; j++) {
							rhs += -input.processTime[j] * XhkValue[j][i];
							expr = lowerPro.sum(expr,
									lowerPro.prod(Zh[j], input.uncertain * input.processTime[j] * XhkValue[j][i]));
							expr = lowerPro.sum(expr, lowerPro.prod(Yhk[j][i], 1));
						}
						rhs += -(betaValue[i] - alphaValue[i]) * input.traverseTime;
						expr = lowerPro.sum(expr, lowerPro.prod(CmaxLower[0], -1));
						lowerPro.remove(lower2_3[i]);
						lower2_3[i] = lowerPro.addLe(expr, rhs, "Constraint2_3 " + i);
					}

					// 2-14
					for (int i = 0; i < jValue; i++) {
						for (int j = 0; j < input.numQC - 1; j++) {
							IloNumExpr expr = lowerPro.numExpr();
							double rhs = 0;
							for (int k = 0; k <= i; k++) {
								rhs += -XhkValue[k][j] * input.processTime[k];
								expr = lowerPro.sum(expr,
										lowerPro.prod(Zh[k], input.uncertain * input.processTime[k] * XhkValue[k][j]));
								expr = lowerPro.sum(expr, lowerPro.prod(Yhk[k][j], 1));
							}

							for (int k = input.safetyMargin + 1; k <= i + input.safetyMargin + 1; k++) {
								rhs += XhkValue[k][j + 1] * input.processTime[k];
								expr = lowerPro.sum(expr, lowerPro.prod(Zh[k],
										-input.uncertain * input.processTime[k] * XhkValue[k][j + 1]));
								expr = lowerPro.sum(expr, lowerPro.prod(Yhk[k][j + 1], -1));
							}
							rhs += (alphaValue[j] - alphaValue[j + 1] + input.safetyMargin + 1) * input.traverseTime;
							expr = lowerPro.sum(expr, lowerPro.prod(Gammajk[i][j], -M));
							lowerPro.remove(lower2_14[i][j]);
							lower2_14[i][j] = lowerPro.addLe(expr, rhs, "Constraint2_14 " + i + "," + j);
						}
					}

					// 2-16
					for (int i = 0; i < input.numQC; i++) {
						IloNumExpr expr = lowerPro.numExpr();
						double rhs = 0;
						for (int j = 0; j < input.numBay; j++) {
							rhs += input.processTime[j] * XhkValue[j][i];
							expr = lowerPro.sum(expr,
									lowerPro.prod(Zh[j], -input.uncertain * input.processTime[j] * XhkValue[j][i]));
							expr = lowerPro.sum(expr, lowerPro.prod(Yhk[j][i], -1));
						}
						rhs += (betaValue[i] - alphaValue[i]) * input.traverseTime;
						expr = lowerPro.s(um(expr, lowerPro.prod(CmaxLower[0], 1));
						expr = lowerPro.sum(expr, lowerPro.prod(Epsilonk[i], -M));
						lowerPro.remove(lower2_16[i]);
						lower2_16[i] = lowerPro.addLe(expr, rhs, "Constraint2_16 " + i);
					}
				} else if (method.equals("SD")) {

					// 对子问题目标函数值进行更新
					IloNumExpr dualObj = lowerPro.numExpr();
					for (int i = 0; i < input.numBay - input.safetyMargin - 1; i++) {
						for (int j = 0; j < input.numQC - 1; j++) {
							for (int j2 = 0; j2 <= i; j2++) {
								dualObj = lowerPro.sum(dualObj,
										lowerPro.prod(Niujk[i][j], -input.processTime[j2] * XhkValue[j2][j]));
								dualObj = lowerPro.sum(dualObj, lowerPro.prod(Whjk[j2][i][j],
										-input.processTime[j2] * XhkValue[j2][j] * input.uncertain));
							}

							for (int j2 = input.safetyMargin + 1; j2 <= i + input.safetyMargin + 1; j2++) {
								dualObj = lowerPro.sum(dualObj,
										lowerPro.prod(Niujk[i][j], input.processTime[j2] * XhkValue[j2][j + 1]));
								dualObj = lowerPro.sum(dualObj, lowerPro.prod(Whjk[j2][i][j],
										input.processTime[j2] * XhkValue[j2][j + 1] * input.uncertain));
							}
							dualObj = lowerPro.sum(dualObj, lowerPro.prod(Niujk[i][j],
									(input.safetyMargin + 1 + alphaValue[j] - alphaValue[j + 1]) * input.traverseTime));
						}
					}

					for (int i = 0; i < input.numQC; i++) {
						for (int j = 0; j < input.numBay; j++) {
							dualObj = lowerPro.sum(dualObj,
									lowerPro.prod(Miuk[i], input.processTime[j] * XhkValue[j][i]));
							dualObj = lowerPro.sum(dualObj,
									lowerPro.prod(Ehk[j][i], input.processTime[j] * XhkValue[j][i] * input.uncertain));
						}
						dualObj = lowerPro.sum(dualObj,
								lowerPro.prod(Miuk[i], (betaValue[i] - alphaValue[i]) * input.traverseTime));
					}
					lowerProObj.setExpr(dualObj);
					lowerPro.exportModel("SP(SD).lp");
				} else {
					System.out.println("Method参数错误！");
				}

				// 对新的下层问题进行求解
				if (lowerPro.solve()) {
					endTime = System.currentTimeMillis();
					System.out.println("下层问题第" + Iteration + "次求解成功了");
					if (method.equals("KKT")) {
						input.writerOverall
								.write("Iteration " + Iteration + " lowerObj(KKT):" + lowerPro.getObjValue() + "\r\n");
						input.writerOverall.write("\r\n");

						if (iterationDetail) {
							lowerPro.writeSolution("./result of benchmark/test_" + input.testID + "/test_"
									+ input.testID + " iteration_" + Iteration + "lowerKKT.txt");
						}

						upperBound = lowerPro.getValue(CmaxLower[0]);

						if (variableOutput) {
							for (int i = 0; i < input.numBay; i++) {
								for (int j = 0; j < input.numQC; j++) {
									System.out.println("Y" + i + "," + j + ": Value = " + lowerPro.getValue(Yhk[i][j]));
								}
							}
						}

					} else if (method.equals("SD")) {
						input.writerOverall
								.write("Iteration " + Iteration + " lowerObj(SD):" + lowerPro.getObjValue() + "\r\n");
						input.writerOverall.write("\r\n");

						if (iterationDetail) {
							lowerPro.writeSolution("./result of benchmark/test_" + input.testID + "/test_"
									+ input.testID + " iteration_" + Iteration + "lowerSD.txt");
						}

						upperBound = lowerPro.getObjValue();
						System.out.println("对偶变量的值待输出");
					} else {
						System.out.println("Method参数错误！");
					}

					System.out.println("上界 = " + upperBound);

					if (upperBound <= bestupperBound) {
						bestupperBound = upperBound;
						ZhValue = lowerPro.getValues(Zh);
						// 先清空之前的内容
						input.resetWriterSenario();
						senarioRecord(input);
					}

					if (variableOutput) {
						for (int i = 0; i < input.numBay; i++) {
							System.out.println("Z" + i + ": Value = " + lowerPro.getValue(Zh[i]));
						}
					}

				} else {
					System.out.println("下层问题未成功求解");
					System.exit(0);
				}
			}
			long totalComputingTime = endTime - startTime;
			input.writerOverall.write("Total computing time:" + totalComputingTime + "\r\n");
			input.writerOverall.write("Iteration:" + Iteration);
			input.writerOverall.close();
			input.writerAssignment.close();
			input.writerSenario.close();
		} catch (IloException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	public void BDIteration(Input input, boolean variableOutput, boolean iterationDetail) throws IOException {
		try {
			double lowerBound = upperPro.getObjValue();
			double upperBound = lowerPro.getObjValue();
			double bestupperBound = upperBound; // 最好的上界
			// 获取一些变量值
			double[][] NiujklValue = new double[jValue][input.numQC - 1];
			double[][][] WhjklValue = new double[input.numBay][jValue][input.numQC - 1];
			double[] MiuklValue = new double[input.numQC];
			double[][] EhklValue = new double[input.numBay][input.numQC];
			IloRange[] upper1_15 = new IloRange[1];

			while (upperBound - lowerBound >= ACCURANCY && (endTime - startTime) / 1000 < terminalSpan) {
				Iteration += 1;
				System.out.println("！！！第" + Iteration + "次加入cut！！！");

				// 获取变量取值
				for (int i = 0; i < jValue; i++) {
					NiujklValue[i] = lowerPro.getValues(Niujk[i]);
				}

				for (int i = 0; i < input.numBay; i++) {
					for (int j = 0; j < jValue; j++) {
						WhjklValue[i][j] = lowerPro.getValues(Whjk[i][j]);
					}
				}

				MiuklValue = lowerPro.getValues(Miuk);

				for (int i = 0; i < input.numBay; i++) {
					EhklValue[i] = lowerPro.getValues(Ehk[i]);
				}

				// 设置约束
				IloNumExpr expr = upperPro.numExpr();
				expr = upperPro.sum(expr, upperPro.prod(-1, CmaxUpper[0]));

				double constantTerm = 0;
				for (int i = 0; i < jValue; i++) {
					for (int j = 0; j < input.numQC - 1; j++) {
						for (int j2 = 0; j2 <= i; j2++) {
							double coefficient = input.processTime[j2] * NiujklValue[i][j];
							coefficient += input.uncertain * input.processTime[j2] * WhjklValue[j2][i][j];
							coefficient = -coefficient;
							expr = upperPro.sum(expr, upperPro.prod(coefficient, Xhk[j2][j]));
						}

						for (int j2 = input.safetyMargin + 1; j2 <= i + input.safetyMargin + 1; j2++) {
							double coefficient = input.processTime[j2] * NiujklValue[i][j];
							coefficient += input.uncertain * input.processTime[j2] * WhjklValue[j2][i][j];
							expr = upperPro.sum(expr, upperPro.prod(coefficient, Xhk[j2][j + 1]));
						}

						expr = upperPro.sum(expr, upperPro.prod(input.traverseTime * NiujklValue[i][j], alpha[j]));
						expr = upperPro.sum(expr, upperPro.prod(-input.traverseTime * NiujklValue[i][j], alpha[j + 1]));
						constantTerm += -input.traverseTime * NiujklValue[i][j] * (input.safetyMargin + 1);
					}
				}

				for (int i = 0; i < input.numQC; i++) {
					for (int j = 0; j < input.numBay; j++) {
						double coefficient = input.processTime[j] * MiuklValue[i];
						coefficient += input.uncertain * input.processTime[j] * EhklValue[j][i];
						expr = upperPro.sum(expr, upperPro.prod(coefficient, Xhk[j][i]));
					}

					expr = upperPro.sum(expr, upperPro.prod(input.traverseTime * MiuklValue[i], beta[i]));
					expr = upperPro.sum(expr, upperPro.prod(-input.traverseTime * MiuklValue[i], alpha[i]));
				}
				upper1_15[0] = upperPro.addLe(expr, constantTerm, "Constraint1_15 ");

				SpeciafiedValueAbortCallback checkIntegerSolution = new SpeciafiedValueAbortCallback(upperPro,
						lowerBound);
				upperPro.clearCallbacks();
				upperPro.use(checkIntegerSolution);

				// 根据设置的算法终止时间计算本次迭代的最长时间
				double leftTime = terminalSpan - (System.currentTimeMillis() - startTime) / 1000;
				upperPro.setParam(IloCplex.DoubleParam.TiLim, leftTime);

				// 对新的上层问题进行求解
				if (upperPro.solve()) {
					input.writerOverall
							.write("Iteration " + Iteration + " upperObj:" + upperPro.getObjValue() + "\r\n");

					if (iterationDetail) {
						upperPro.writeSolution("./result of benchmark/test_" + input.testID + "/test_" + input.testID
								+ " iteration_" + Iteration + "upper.txt");
					}

					System.out.println("上层问题第" + Iteration + "次求解成功了");
					lowerBound = upperPro.getValue(CmaxUpper[0]);
					System.out.println("下界 = " + lowerBound);

					// 更新X值
					for (int i = 0; i < XhkValue.length; i++) {
						for (int j = 0; j < XhkValue[0].length; j++) {
							XhkValue[i][j] = upperPro.getValue(Xhk[i][j]);
						}
					}

					// 跟新alpha和beta
					alphaValue = upperPro.getValues(alpha);
					betaValue = upperPro.getValues(beta);

					if (variableOutput) {
						for (int i = 0; i < input.numBay; i++) {
							for (int j = 0; j < input.numQC; j++) {
								System.out.println("X" + i + "," + j + ": Value = " + XhkValue[i][j]);
							}
						}

						for (int i = 0; i < input.numQC; i++) {
							System.out.println("alpha " + i + ": Value = " + alphaValue[i]);
						}

						for (int i = 0; i < input.numQC; i++) {
							System.out.println("beta " + i + ": Value = " + betaValue[i]);
						}
					}
				} else {
					System.out.println("上层问题未成功求解");
					System.exit(0);
				}

				// 对子问题目标函数值进行更新
				IloNumExpr dualObj = lowerPro.numExpr();
				for (int i = 0; i < input.numBay - input.safetyMargin - 1; i++) {
					for (int j = 0; j < input.numQC - 1; j++) {
						for (int j2 = 0; j2 <= i; j2++) {
							dualObj = lowerPro.sum(dualObj,
									lowerPro.prod(Niujk[i][j], -input.processTime[j2] * XhkValue[j2][j]));
							dualObj = lowerPro.sum(dualObj, lowerPro.prod(Whjk[j2][i][j],
									-input.processTime[j2] * XhkValue[j2][j] * input.uncertain));
						}

						for (int j2 = input.safetyMargin + 1; j2 <= i + input.safetyMargin + 1; j2++) {
							dualObj = lowerPro.sum(dualObj,
									lowerPro.prod(Niujk[i][j], input.processTime[j2] * XhkValue[j2][j + 1]));
							dualObj = lowerPro.sum(dualObj, lowerPro.prod(Whjk[j2][i][j],
									input.processTime[j2] * XhkValue[j2][j + 1] * input.uncertain));
						}

						dualObj = lowerPro.sum(dualObj, lowerPro.prod(Niujk[i][j],
								(input.safetyMargin + 1 + alphaValue[j] - alphaValue[j + 1]) * input.traverseTime));
					}
				}

				for (int i = 0; i < input.numQC; i++) {
					for (int j = 0; j < input.numBay; j++) {
						dualObj = lowerPro.sum(dualObj, lowerPro.prod(Miuk[i], input.processTime[j] * XhkValue[j][i]));
						dualObj = lowerPro.sum(dualObj,
								lowerPro.prod(Ehk[j][i], input.processTime[j] * XhkValue[j][i] * input.uncertain));
					}
					dualObj = lowerPro.sum(dualObj,
							lowerPro.prod(Miuk[i], (betaValue[i] - alphaValue[i]) * input.traverseTime));
				}
				lowerProObj.setExpr(dualObj);
//				lowerPro.exportModel("SP(SD).lp");

				// 对新的下层问题进行求解
				if (lowerPro.solve()) {
					endTime = System.currentTimeMillis();
					System.out.println("下层问题第" + Iteration + "次求解成功了");

					input.writerOverall
							.write("Iteration " + Iteration + " lowerObj(SD):" + lowerPro.getObjValue() + "\r\n");
					input.writerOverall.write("\r\n");
//					lowerPro.writeSolution("./result of benchmark/test_" + input.testID + "/test_" + input.testID + " iteration_" + Iteration + "lowerSD.txt");
					upperBound = lowerPro.getObjValue();
					System.out.println("上界 = " + upperBound);

					if (upperBound <= bestupperBound) {
						bestupperBound = upperBound;
						ZhValue = lowerPro.getValues(Zh);
						
						// 记录X的取值，也就是泊位和岸桥的分配关系，还有alpha和beta的取值，需要先对内容进行清空
						input.resetWriterAssignment();
						assignmentRecord(input);
						
						// 先清空之前的内容,记录发生uncertainty的贝位索引
						input.resetWriterSenario();
						senarioRecord(input);
					}

					if (variableOutput) {
						for (int i = 0; i < input.numBay; i++) {
							System.out.println("Z" + i + ": Value = " + lowerPro.getValue(Zh[i]));
						}
					}
				} else {
					System.out.println("下层问题未成功求解");
					System.exit(0);
				}
			}
			long totalComputingTime = endTime - startTime;
			input.writerOverall.write("Total computing time:" + totalComputingTime + "\r\n");
			input.writerOverall.write("Iteration:" + Iteration);
			input.writerOverall.close();
			input.writerAssignment.close();
			input.writerSenario.close();
		} catch (IloException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	// variableOutput表示是否在Console输出相关的变量取值
	// iterationDetail表示是否保存迭代过程中的模型求解信息
	public void CCGBDIteration(Input input, boolean variableOutput, boolean iterationDetail) throws IOException {
		try {
			double lowerBound = upperPro.getObjValue();
			double upperBound = lowerPro.getObjValue();
			double bestupperBound = upperBound;

			double[] ZhValue = new double[input.numBay];

			while (upperBound - lowerBound >= ACCURANCY && (endTime - startTime) / 1000 < terminalSpan) {
				Iteration += 1;
				System.out.println("！！！第" + Iteration + "次加入cut！！！");

				ZhValue = lowerPro.getValues(Zh);

				// 确定Zh变量的取值
				int[] ZhSet = new int[input.budget];
				int index = 0;

				for (int i = 0; i < ZhValue.length; i++) {
					if (Math.abs(ZhValue[i] - 1) <= FLOATACCURACY) {
						ZhSet[index] = i;
						index += 1;
					}
				}

				UncertainSenario senario = null;
				for (UncertainSenario temp : uncertainSenario) {
					if (StaticMethod.intListEqual(ZhSet, temp.ZhSet)) {
						temp.frequency += 1;
						senario = temp;
						break;
					}
				}

				if (senario == null) {
					senario = new UncertainSenario();
					senario.ZhSet = ZhSet;
					senario.frequency += 1;
					uncertainSenario.add(senario);
				}

				if (senario.frequency >= changeThreshold) {
					String txtPath = "./result of benchmark/uncertaintyRecord";
					String txtName = "uncertainRecordInstance" + input.testID + ".txt";
					File file = new File(txtPath);
					if (!file.exists()) {// 如果文件夹不存在
						file.mkdir();// 创建文件夹
					}

//					//如果文件已经存在，则需要创建新的空文件
//					File file1 = new File(txtPath + "/" + txtName);
//					if (file1.exists()) {
//						StaticMethod.txtFileCreatandWrite(txtPath, txtName);
//					}

					StringBuilder uncertaintyRecord = new StringBuilder();
					for (int i = 0; i < senario.ZhSet.length; i++) {
						uncertaintyRecord.append(senario.ZhSet[i]);
						uncertaintyRecord.append("\t");
					}

					StaticMethod.txtFileCreatandWrite(txtPath, txtName, uncertaintyRecord + "\r\n");
					// 需要先将之前benders加入的同样的uncertainty senario的cut去除
					for (IloRange removeCut : senario.upper1_15Iteration) {
						upperPro.remove(removeCut);
					}

					uncertainSenario.remove(senario);

					// C&CG
					IterationCcg += 1;
					System.out.println("！！！第" + IterationCcg + "次加入C&CG cut！！！");

					// 创建变量
					IloNumVar[][] Yhkl = new IloNumVar[input.numBay][input.numQC];
					// 创建约束
					IloRange[][] upper1_10 = new IloRange[jValue][input.numQC - 1];
					IloRange[] upper1_11 = new IloRange[input.numQC];

					// 设置变量
					for (int i = 0; i < input.numBay; i++) {
						for (int j = 0; j < input.numQC; j++) {
							Yhkl[i][j] = upperPro.numVar(0, Double.MAX_VALUE, IloNumVarType.Float,
									"Yhkl " + i + "," + j);
						}
					}

					// 设置约束
					// 1-10
					for (int i = 0; i < jValue; i++) {
						for (int j = 0; j < input.numQC - 1; j++) {
							IloNumExpr expr = upperPro.numExpr();
							for (int k = 0; k < i + 1; k++) {
								expr = upperPro.sum(expr, upperPro.prod(Xhk[k][j],
										input.processTime[k] + input.uncertain * input.processTime[k] * ZhValue[k]));
								expr = upperPro.sum(expr, upperPro.prod(Yhkl[k][j], 1));
							}
							expr = upperPro.sum(expr, upperPro.prod(alpha[j], -input.traverseTime));

							for (int k = input.safetyMargin + 1; k < i + input.safetyMargin + 2; k++) {
								expr = upperPro.sum(expr, upperPro.prod(Xhk[k][j + 1],
										-(input.processTime[k] + input.uncertain * input.processTime[k] * ZhValue[k])));
								expr = upperPro.sum(expr, upperPro.prod(Yhkl[k][j + 1], -1));
							}
							expr = upperPro.sum(expr, upperPro.prod(alpha[j + 1], input.traverseTime));
							upper1_10[i][j] = upperPro.addGe(expr, (input.safetyMargin + 1) * input.traverseTime,
									"Constraint1_10 " + i + "," + j);
						}
					}

					// 1-11
					for (int i = 0; i < input.numQC; i++) {
						IloNumExpr expr = upperPro.numExpr();
						for (int j = 0; j < input.numBay; j++) {
							expr = upperPro.sum(expr, upperPro.prod(Xhk[j][i],
									input.processTime[j] + input.uncertain * input.processTime[j] * ZhValue[j]));
							expr = upperPro.sum(expr, upperPro.prod(Yhkl[j][i], 1));
						}
						expr = upperPro.sum(expr, upperPro.prod(beta[i], input.traverseTime));
						expr = upperPro.sum(expr, upperPro.prod(alpha[i], -input.traverseTime));
						expr = upperPro.sum(expr, upperPro.prod(CmaxUpper[0], -1));
						upper1_11[i] = upperPro.addLe(expr, 0, "Constraint1_11 " + i);
					}
				} else {
					// Benders
					// 用于保留变量取值
					double[][] NiujklValue = new double[jValue][input.numQC - 1];
					double[][][] WhjklValue = new double[input.numBay][jValue][input.numQC - 1];
					double[] MiuklValue = new double[input.numQC];
					double[][] EhklValue = new double[input.numBay][input.numQC];
					IloRange[] upper1_15 = new IloRange[1];

					// 获取变量取值
					IterationBenders += 1;
					System.out.println("！！！第" + IterationBenders + "次加入Benders cut！！！");

					for (int i = 0; i < jValue; i++) {
						NiujklValue[i] = lowerPro.getValues(Niujk[i]);
					}

					for (int i = 0; i < input.numBay; i++) {
						for (int j = 0; j < jValue; j++) {
							WhjklValue[i][j] = lowerPro.getValues(Whjk[i][j]);
						}
					}

					MiuklValue = lowerPro.getValues(Miuk);

					for (int i = 0; i < input.numBay; i++) {
						EhklValue[i] = lowerPro.getValues(Ehk[i]);
					}

					// 设置约束
					IloNumExpr expr = upperPro.numExpr();
					expr = upperPro.sum(expr, upperPro.prod(-1, CmaxUpper[0]));

					double constantTerm = 0;
					for (int i = 0; i < jValue; i++) {
						for (int j = 0; j < input.numQC - 1; j++) {
							for (int j2 = 0; j2 <= i; j2++) {
								double coefficient = input.processTime[j2] * NiujklValue[i][j];
								coefficient += input.uncertain * input.processTime[j2] * WhjklValue[j2][i][j];
								coefficient = -coefficient;
								expr = upperPro.sum(expr, upperPro.prod(coefficient, Xhk[j2][j]));
							}

							for (int j2 = input.safetyMargin + 1; j2 <= i + input.safetyMargin + 1; j2++) {
								double coefficient = input.processTime[j2] * NiujklValue[i][j];
								coefficient += input.uncertain * input.processTime[j2] * WhjklValue[j2][i][j];
								expr = upperPro.sum(expr, upperPro.prod(coefficient, Xhk[j2][j + 1]));
							}

							expr = upperPro.sum(expr, upperPro.prod(input.traverseTime * NiujklValue[i][j], alpha[j]));
							expr = upperPro.sum(expr,
									upperPro.prod(-input.traverseTime * NiujklValue[i][j], alpha[j + 1]));
							constantTerm += -input.traverseTime * NiujklValue[i][j] * (input.safetyMargin + 1);
						}
					}

					for (int i = 0; i < input.numQC; i++) {
						for (int j = 0; j < input.numBay; j++) {
							double coefficient = input.processTime[j] * MiuklValue[i];
							coefficient += input.uncertain * input.processTime[j] * EhklValue[j][i];
							expr = upperPro.sum(expr, upperPro.prod(coefficient, Xhk[j][i]));
						}

						expr = upperPro.sum(expr, upperPro.prod(input.traverseTime * MiuklValue[i], beta[i]));
						expr = upperPro.sum(expr, upperPro.prod(-input.traverseTime * MiuklValue[i], alpha[i]));
					}
					upper1_15[0] = upperPro.addLe(expr, constantTerm, "Constraint1_15 ");
					senario.upper1_15Iteration.add(upper1_15[0]);
				}

				SpeciafiedValueAbortCallback checkIntegerSolution = new SpeciafiedValueAbortCallback(upperPro,
						lowerBound);
				upperPro.clearCallbacks();
				upperPro.use(checkIntegerSolution);

				// 根据设置的算法终止时间计算本次迭代的最长时间
				double leftTime = terminalSpan - (System.currentTimeMillis() - startTime) / 1000;
				upperPro.setParam(IloCplex.DoubleParam.TiLim, leftTime);
				// 对新的上层问题进行求解
				if (upperPro.solve()) {
					input.writerOverall
							.write("Iteration " + Iteration + " upperObj:" + upperPro.getObjValue() + "\r\n");
					if (iterationDetail) {
						upperPro.writeSolution("./result of benchmark/test_" + input.testID + "/test_" + input.testID
								+ " iteration_" + Iteration + "upper.txt");
					}

					System.out.println("上层问题第" + Iteration + "次求解成功了");
					lowerBound = upperPro.getValue(CmaxUpper[0]);
					System.out.println("下界 = " + lowerBound);
					// 更新X值
					for (int i = 0; i < XhkValue.length; i++) {
						for (int j = 0; j < XhkValue[0].length; j++) {
							XhkValue[i][j] = upperPro.getValue(Xhk[i][j]);
						}
					}

					// 跟新alpha和beta
					alphaValue = upperPro.getValues(alpha);
					betaValue = upperPro.getValues(beta);

					// 记录X的取值，也就是泊位和岸桥的分配关系，还有alpha和beta的取值，需要先对内容进行清空
					input.resetWriterAssignment();
					assignmentRecord(input);

					if (variableOutput) {
						for (int i = 0; i < input.numBay; i++) {
							for (int j = 0; j < input.numQC; j++) {
								System.out.println("X" + i + "," + j + ": Value = " + XhkValue[i][j]);
							}
						}

						for (int i = 0; i < input.numQC; i++) {
							System.out.println("alpha " + i + ": Value = " + alphaValue[i]);
						}

						for (int i = 0; i < input.numQC; i++) {
							System.out.println("beta " + i + ": Value = " + betaValue[i]);
						}
					}
				} else {
					System.out.println("上层问题未成功求解");
					System.exit(0);
				}

				// 对子问题目标函数值进行更新
				IloNumExpr dualObj = lowerPro.numExpr();
				for (int i = 0; i < input.numBay - input.safetyMargin - 1; i++) {
					for (int j = 0; j < input.numQC - 1; j++) {
						for (int j2 = 0; j2 <= i; j2++) {
							dualObj = lowerPro.sum(dualObj,
									lowerPro.prod(Niujk[i][j], -input.processTime[j2] * XhkValue[j2][j]));
							dualObj = lowerPro.sum(dualObj, lowerPro.prod(Whjk[j2][i][j],
									-input.processTime[j2] * XhkValue[j2][j] * input.uncertain));
						}

						for (int j2 = input.safetyMargin + 1; j2 <= i + input.safetyMargin + 1; j2++) {
							dualObj = lowerPro.sum(dualObj,
									lowerPro.prod(Niujk[i][j], input.processTime[j2] * XhkValue[j2][j + 1]));
							dualObj = lowerPro.sum(dualObj, lowerPro.prod(Whjk[j2][i][j],
									input.processTime[j2] * XhkValue[j2][j + 1] * input.uncertain));
						}
						dualObj = lowerPro.sum(dualObj, lowerPro.prod(Niujk[i][j],
								(input.safetyMargin + 1 + alphaValue[j] - alphaValue[j + 1]) * input.traverseTime));
					}
				}

				for (int i = 0; i < input.numQC; i++) {
					for (int j = 0; j < input.numBay; j++) {
						dualObj = lowerPro.sum(dualObj, lowerPro.prod(Miuk[i], input.processTime[j] * XhkValue[j][i]));
						dualObj = lowerPro.sum(dualObj,
								lowerPro.prod(Ehk[j][i], input.processTime[j] * XhkValue[j][i] * input.uncertain));
					}
					dualObj = lowerPro.sum(dualObj,
							lowerPro.prod(Miuk[i], (betaValue[i] - alphaValue[i]) * input.traverseTime));
				}
				lowerProObj.setExpr(dualObj);
//				lowerPro.exportModel("SP(SD).lp");

				// 对新的下层问题进行求解
				if (lowerPro.solve()) {
					endTime = System.currentTimeMillis();
					System.out.println("下层问题第" + Iteration + "次求解成功了");

					input.writerOverall
							.write("Iteration " + Iteration + " lowerObj(SD):" + lowerPro.getObjValue() + "\r\n");
					input.writerOverall.write("\r\n");
					if (iterationDetail) {
						lowerPro.writeSolution("./result of benchmark/test_" + input.testID + "/test_" + input.testID
								+ " iteration_" + Iteration + "lowerSD.txt");
					}
					upperBound = lowerPro.getObjValue();
					System.out.println("上界 = " + upperBound);

					if (upperBound <= bestupperBound) {
						bestupperBound = upperBound;
						ZhValue = lowerPro.getValues(Zh);
						// 先清空之前的内容
						input.resetWriterSenario();
						senarioRecord(input);
					}

					if (variableOutput) {
						for (int i = 0; i < input.numBay; i++) {
							System.out.println("Z" + i + ": Value = " + lowerPro.getValue(Zh[i]));
						}
					}
				} else {
					System.out.println("下层问题未成功求解");
					System.exit(0);
				}
			}
			long totalComputingTime = endTime - startTime;
			input.writerOverall.write("Total computing time:" + totalComputingTime);
			input.writerOverall.write("\r\n");
			input.writerOverall.write("Iteration:" + Iteration);
			input.writerOverall.write("\r\n");
			input.writerOverall.write("BendersCut:" + IterationBenders);
			input.writerOverall.write("\r\n");
			input.writerOverall.write("C&CGCut:" + IterationCcg);
			input.writerOverall.close();
			input.writerAssignment.close();
			input.writerSenario.close();
		} catch (IloException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

	// variableOutput表示是否在Console输出相关的变量取值
	public void deterministicModel(Input input, boolean variableOutput, boolean resultSave) throws IOException {
		try {
			startTime = System.currentTimeMillis();
			upperPro = new IloCplex();
			int jValue = input.numBay - input.safetyMargin - 1;
//			upperPro.setOut(null); // 关闭输出日志
			// 变量
			Xhk = new IloNumVar[input.numBay][input.numQC];
			CmaxUpper = new IloNumVar[1];
			alpha = new IloNumVar[input.numQC];
			beta = new IloNumVar[input.numQC];
			Yhk = new IloNumVar[input.numBay][input.numQC];
			// 约束
			upper1_2 = new IloRange[input.numBay];
			upper1_3 = new IloRange[input.numBay][input.numQC];
			upper1_4 = new IloRange[input.numBay][input.numQC];
			upper1_8 = new IloRange[input.numQC];
			upper1_9 = new IloRange[input.numQC];
			upper1_10 = new IloRange[jValue][input.numQC - 1];
			upper1_11 = new IloRange[input.numQC];

			// 设置变量
			for (int i = 0; i < input.numBay; i++) {
				for (int j = 0; j < input.numQC; j++) {
					Xhk[i][j] = upperPro.numVar(0, 1, IloNumVarType.Int, "X" + i + "," + j);
				}
			}

			CmaxUpper[0] = upperPro.numVar(0, Double.MAX_VALUE, IloNumVarType.Float, "CmaxUpper");

			for (int i = 0; i < input.numQC; i++) {
				alpha[i] = upperPro.numVar(1, input.numBay, IloNumVarType.Float, "alpha" + i);
			}

			for (int i = 0; i < input.numQC; i++) {
				beta[i] = upperPro.numVar(1, input.numBay, IloNumVarType.Float, "beta" + i);
			}

			for (int i = 0; i < input.numBay; i++) {
				for (int j = 0; j < input.numQC; j++) {
					Yhk[i][j] = upperPro.numVar(0, Double.MAX_VALUE, IloNumVarType.Float, "Yhk_" + i + "," + j);
				}
			}

			// 目标函数
			IloNumExpr upperObj = upperPro.prod(CmaxUpper[0], 1);
			upperPro.addMinimize(upperObj);

			// 上层问题约束
			// 约束1-2
			for (int i = 0; i < input.numBay; i++) {
				if (input.processTime[i] != 0) {
					IloNumExpr expr = upperPro.numExpr();
					for (int j = 0; j < input.numQC; j++) {
						expr = upperPro.sum(expr, upperPro.prod(Xhk[i][j], 1));
					}
					upper1_2[i] = upperPro.addEq(expr, 1, "constraint1_2 " + i);
				}
			}

			// 约束1-3
			for (int i = 0; i < input.numBay; i++) {
				for (int j = 0; j < input.numQC; j++) {
					IloNumExpr expr = upperPro.numExpr();
					expr = upperPro.sum(expr, upperPro.prod(Xhk[i][j], M - i - 1));
					expr = upperPro.sum(expr, upperPro.prod(alpha[j], 1));
					upper1_3[i][j] = upperPro.addLe(expr, M, "constraint1_3 " + i + "," + j);
				}
			}

			// 约束1-4
			for (int i = 0; i < input.numBay; i++) {
				for (int j = 0; j < input.numQC; j++) {
					IloNumExpr expr = upperPro.numExpr();
					expr = upperPro.sum(expr, upperPro.prod(Xhk[i][j], -M - i - 1));
					expr = upperPro.sum(expr, upperPro.prod(beta[j], 1));
					upper1_4[i][j] = upperPro.addGe(expr, -M, "constraint1_4 " + i + "," + j);
				}
			}

			// 约束1-8
			for (int i = 0; i < input.numQC - 1; i++) {
				IloNumExpr expr = upperPro.numExpr();
				expr = upperPro.sum(expr, upperPro.prod(alpha[i], 1));
				expr = upperPro.sum(expr, upperPro.prod(alpha[i + 1], -1));
				upper1_8[i] = upperPro.addLe(expr, -input.safetyMargin - 1, "constraint1_8 " + i);
			}

			// 约束1-9
			for (int i = 0; i < input.numQC - 1; i++) {
				IloNumExpr expr = upperPro.numExpr();
				expr = upperPro.sum(expr, upperPro.prod(beta[i], 1));
				expr = upperPro.sum(expr, upperPro.prod(beta[i + 1], -1));
				upper1_9[i] = upperPro.addLe(expr, -input.safetyMargin - 1, "constraint1_9 " + i);
			}

			// 约束1-10
			for (int i = 0; i < jValue; i++) {
				for (int j = 0; j < input.numQC - 1; j++) {
					IloNumExpr expr = upperPro.numExpr();
					for (int k = 0; k <= i; k++) {
						expr = upperPro.sum(expr, upperPro.prod(Xhk[k][j], input.processTime[k]));
						expr = upperPro.sum(expr, upperPro.prod(Yhk[k][j], 1));
					}
					for (int k = input.safetyMargin + 1; k < i + input.safetyMargin + 2; k++) {
						expr = upperPro.sum(expr, upperPro.prod(Xhk[k][j + 1], -input.processTime[k]));
						expr = upperPro.sum(expr, upperPro.prod(Yhk[k][j + 1], -1));
					}
					expr = upperPro.sum(expr, upperPro.prod(alpha[j], -input.traverseTime));
					expr = upperPro.sum(expr, upperPro.prod(alpha[j + 1], input.traverseTime));
					upperPro.addGe(expr, (input.safetyMargin + 1) * input.traverseTime,
							"Constraint1_10 " + i + "," + j);
				}
			}

			// 约束1-11
			for (int i = 0; i < input.numQC; i++) {
				IloNumExpr expr = upperPro.numExpr();
				for (int j = 0; j < input.numBay; j++) {
					expr = upperPro.sum(expr, upperPro.prod(Xhk[j][i], input.processTime[j]));
					expr = upperPro.sum(expr, upperPro.prod(Yhk[j][i], 1));
				}
				expr = upperPro.sum(expr, upperPro.prod(beta[i], input.traverseTime));
				expr = upperPro.sum(expr, upperPro.prod(alpha[i], -input.traverseTime));
				expr = upperPro.sum(expr, upperPro.prod(CmaxUpper[0], -1));
				upperPro.addLe(expr, 0, "Constraint1_11 " + i);
			}

//			upperPro.exportModel("determine.lp");
			if (upperPro.solve()) {
				endTime = System.currentTimeMillis();
				double totalTime = endTime - startTime;
				File file = new File("./determinstic");

				if (!file.exists()) {// 如果文件夹不存在
					file.mkdir();// 创建文件夹
				}

				if (resultSave) {
					upperPro.writeSolution("./determinstic/test_" + input.testID + ".txt");
					FileWriter writer = new FileWriter("./determinstic/test_" + input.testID + ".txt", true);
					writer.write("运行的时间为：" + totalTime);
					writer.close();
				}

				StaticMethod.txtFileCreatandWrite("./determinstic", "deterministicOverall.txt",
						input.testID + "\t" + upperPro.getObjValue() + "\t" + totalTime + "\n");
				System.out.println("求解成功了");

				XhkValue = new double[input.numBay][input.numQC];
				alphaValue = new double[input.numQC];
				betaValue = new double[input.numQC];

				System.out.println("Solution status = " + upperPro.getStatus());
				System.out.println("Solution value  = " + upperPro.getObjValue());

				for (int i = 0; i < input.numBay; i++) {
					for (int j = 0; j < input.numQC; j++) {
						XhkValue[i][j] = upperPro.getValue(Xhk[i][j]);
					}
				}
				alphaValue = upperPro.getValues(alpha);
				betaValue = upperPro.getValues(beta);

				if (variableOutput) {
					for (int i = 0; i < input.numBay; i++) {
						for (int j = 0; j < input.numQC; j++) {
							System.out.println("X" + i + "," + j + ": Value = " + XhkValue[i][j]);
						}
					}

					for (int i = 0; i < input.numQC; i++) {
						System.out.println("alpha" + i + ": Value = " + alphaValue[i]);
					}

					for (int i = 0; i < input.numQC; i++) {
						System.out.println("beta" + i + ": Value = " + betaValue[i]);
					}

					for (int i = 0; i < input.numBay; i++) {
						for (int j = 0; j < input.numQC; j++) {
							System.out.println("Y" + i + "," + j + ": Value = " + upperPro.getValue(Yhk[i][j]));
						}
					}
				}
			} else {
				System.out.println("确定性问题未成功求解");
				System.exit(0);
			}
		} catch (IloException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	// 必须在deterministicModel() 运行后再运行
	public void deterministicAttackModel(Input input, boolean resultSave) {
		try {
			lowerModelSD(input, false, false, false);

			if (resultSave) {
				lowerPro.writeSolution("./determinstic/attackTest_" + input.testID + ".txt");
			}

			String tempString = "";
			for (int i = 0; i < ZhValue.length; i++) {
				if (Math.abs(ZhValue[i] - 1) <= FLOATACCURACY) {
					tempString += "\t" + i;
				}
			}
			tempString += "\n";

			StaticMethod.txtFileCreatandWrite("./determinstic", "deterministicAttackOverall.txt",
					input.testID + "\t" + String.format("%.2f", lowerPro.getObjValue()) + tempString);
		} catch (IOException | IloException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	// 当以确定性模型作为初始解时，求完之后需要剔除相关约束
	public void deterministicInitial(Input input, boolean variableOutput, boolean resultSave) throws IOException {
		deterministicModel(input, variableOutput, resultSave);

		try {
			input.writerOverall.write("Iteration " + Iteration + " upperObj:" + upperPro.getObjValue() + "\r\n");
		} catch (IloException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}

		int jValue = input.numBay - input.safetyMargin - 1;
		try {
			// 约束1-10
			for (int i = 0; i < jValue; i++) {
				for (int j = 0; j < input.numQC - 1; j++) {
					upperPro.remove(upper1_10[i][j]);
				}
			}

			// 约束1-11
			for (int i = 0; i < input.numQC; i++) {
				upperPro.remove(upper1_11[i]);
			}
		} catch (IloException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

	// 固定代表贝位和岸桥指派关系的X变量，固定代表uncertainty的Z变量后，求解的问题
	// 相当于求解固定uncertainty后最小化的下层问题
	public void assignmentUncertaintyFix(Input input) throws IOException {
		try {
			lowerPro = new IloCplex();
			jValue = input.numBay - input.safetyMargin - 1;
			// 变量
			Yhk = new IloNumVar[input.numBay][input.numQC];
			CmaxLower = new IloNumVar[1];
			// 约束
			lower2_2 = new IloRange[jValue][input.numQC];
			lower2_3 = new IloRange[input.numQC];

			// 设置变量
			for (int i = 0; i < input.numBay; i++) {
				for (int j = 0; j < input.numQC; j++) {
					Yhk[i][j] = lowerPro.numVar(0, Double.MAX_VALUE, IloNumVarType.Float, "Yhk_" + i + "," + j);
				}
			}

			CmaxLower[0] = lowerPro.numVar(0, Double.MAX_VALUE, IloNumVarType.Float, "CmaxLower");

			
			// 读取X和Z变量的取值
			Schedule scheduleData = new Schedule(input.numBay, input.numQC);
			String assignmentPath = "./result of benchmark/assignment/best assignment of test_" + input.testID + ".txt";
			String uncertaintyPath = "./result of benchmark/senario/uncertainty senario of test_" + input.testID + ".txt";
			scheduleData.dataInput(assignmentPath, uncertaintyPath);
			
			// 设置约束
			// 2-2
			for (int i = 0; i < jValue; i++) {
				for (int j = 0; j < input.numQC - 1; j++) {
					IloNumExpr expr = lowerPro.numExpr();
					double rhs = 0;
					
					for (int k = 0; k <= i; k++) {
						rhs += -scheduleData.XhkValue[k][j] * input.processTime[k] * (1 + input.uncertain * scheduleData.ZhValue[k]);
						expr = lowerPro.sum(expr, lowerPro.prod(Yhk[k][j], 1));
					}
					
					for (int k = input.safetyMargin + 1; k < i + input.safetyMargin + 2; k++) {
						rhs += scheduleData.XhkValue[k][j + 1] * input.processTime[k] * (1 + input.uncertain * scheduleData.ZhValue[k]);
						expr = lowerPro.sum(expr, lowerPro.prod(Yhk[k][j + 1], -1));
					}
					
					rhs += (scheduleData.alphaValue[j] - scheduleData.alphaValue[j + 1] + input.safetyMargin + 1) * input.traverseTime;
					lower2_2[i][j] = lowerPro.addGe(expr, rhs, "Constraint2_2 " + i + "," + j);
				}
			}

			// 2-3
			for (int i = 0; i < input.numQC; i++) {
				IloNumExpr expr = lowerPro.numExpr();
				double rhs = 0;
				for (int j = 0; j < input.numBay; j++) {
					rhs += -scheduleData.XhkValue[j][i] * input.processTime[j] * (1 + input.uncertain * scheduleData.ZhValue[j]);
					expr = lowerPro.sum(expr, lowerPro.prod(Yhk[j][i], 1));
				}
				rhs += -(scheduleData.betaValue[i] - scheduleData.alphaValue[i]) * input.traverseTime;
				expr = lowerPro.sum(expr, lowerPro.prod(CmaxLower[0], -1));
				lower2_3[i] = lowerPro.addLe(expr, rhs, "Constraint2_3 " + i);
			}
			
			// 设置目标函数
			// 原问题目标函数
			IloNumExpr obj = lowerPro.numExpr();
			obj = lowerPro.sum(obj, CmaxLower[0]);
			lowerProObj = lowerPro.addMinimize(obj);
			
			if (lowerPro.solve()) {
				System.out.println("Solution status = " + lowerPro.getStatus());
				System.out.println("Solution value  = " + lowerPro.getObjValue());
				
			}
		} catch (IloException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}

	// 用于记录上层问题输出的指派方案，也就是变量X、alpha和beta的取值
	public void assignmentRecord(Input input) {
		try {
			for (int i = 0; i < XhkValue.length; i++) {
				for (int j = 0; j < XhkValue[0].length; j++) {
					input.writerAssignment.write("X" + i + "," + j + ":" + (int) Math.round(XhkValue[i][j]) + "\n");
				}
			}

			for (int i = 0; i < alphaValue.length; i++) {
				input.writerAssignment.write("alpha_" + i + ":" + (int) Math.round(alphaValue[i]) + "\n");
			}

			for (int i = 0; i < betaValue.length; i++) {
				input.writerAssignment.write("beta_" + i + ":" + (int) Math.round(betaValue[i]) + "\n");
			}

		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	// 记录下层问题Attack发生的贝位索引和
	public void senarioRecord(Input input) {
		try {
			for (int i = 0; i < ZhValue.length; i++) {
				input.writerSenario.write("Z" + i + ":" + Math.round(ZhValue[i]) + "\n");
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
}
