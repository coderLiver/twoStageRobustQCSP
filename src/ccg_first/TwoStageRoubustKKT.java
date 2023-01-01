package ccg_first;

import ilog.concert.IloNumExpr;
import ilog.concert.IloNumVar;
import ilog.concert.IloNumVarType;
import ilog.concert.IloObjective;
import ilog.concert.IloRange;
import ilog.cplex.IloCplex;
import ilog.concert.IloException;
import ilog.concert.IloIntVar;

public class TwoStageRoubustKKT {
	IloCplex upperPro;
	IloCplex lowerPro;
	IloObjective lowerProObj;
	// 上层问题变量
	IloNumVar[][][] Xhki;
	IloNumVar[] Cmax1;
	double[][][] XhkiValue;
	// 下层问题变量
	IloNumVar[] beta_h;
	IloNumVar[] lamda_h;
	IloNumVar[][][] gamma_hki;
	IloNumVar[] Ch;
	IloNumVar[] zh;
	IloNumVar[] miu_h;
	IloNumVar[][][] niu_hki;
	IloNumVar[] ksi_h;
	IloNumVar[] oh;
	IloNumVar[] pai;
	IloNumVar[] Cmax;
	// 上层问题约束取值范围
	IloRange[] upper1_2;
	IloRange[][] upper1_3;
	IloRange[][] upper1_4;
	// 下层问题约束取值范围
	IloRange[] lower2_2;
	IloRange[][][] lower2_3;
	IloRange[] lower2_4;
	IloRange[] lower2_5;
	IloRange[] lower2_6;
	IloRange[] lower2_7;
	IloRange[] lower2_8;
	IloRange[][][] lower2_9;
	IloRange[][][] lower2_10;
	IloRange[] lower2_11;
	IloRange[] lower2_12;
	IloRange[] lower2_13;
	IloRange[] lower2_14;
	IloRange[] lower2_15;
	IloRange[] lower2_16;
	IloRange[] lower2_17;

	private static final double M = 1000;
	private static final double M1 = 2000;
	private static final double ACCURANCY = 1;
	int Iteration = 0;

	public void upperModel(int numBay, int numQC) {
		try {
			upperPro = new IloCplex();
			// 变量
			Xhki = new IloNumVar[numBay][numQC][numBay];
			Cmax1 = new IloNumVar[1];
			// 约束
			upper1_2 = new IloRange[numBay];
			upper1_3 = new IloRange[numQC][numBay];
			upper1_4 = new IloRange[numQC][numBay];
			// 设置变量
			for (int i = 0; i < numBay; i++) {
				for (int j = 0; j < numQC; j++) {
					for (int j2 = 0; j2 < numBay; j2++) {
						Xhki[i][j][j2] = upperPro.numVar(0, 1, IloNumVarType.Int, "X" + i + j + j2);
					}
				}
			}

			Cmax1[0] = upperPro.numVar(0, Double.MAX_VALUE, IloNumVarType.Float, "Cmax");

			// 目标函数
			IloNumExpr masterObj = upperPro.prod(Cmax1[0], 1);
			upperPro.addMinimize(masterObj);

			// 上层问题约束
			// 约束1-2
			for (int i = 0; i < numBay; i++) {
				IloNumExpr expr = upperPro.numExpr();
				for (int j = 0; j < numQC; j++) {
					for (int j2 = 0; j2 < numBay; j2++) {
						expr = upperPro.sum(expr, upperPro.prod(Xhki[i][j][j2], 1));
					}
				}
				upper1_2[i] = upperPro.addEq(expr, 1, "constraint1_2 " + i);
			}

			// 约束1-3
			for (int i = 0; i < numQC; i++) {
				for (int j = 0; j < numBay; j++) {
					IloNumExpr expr = upperPro.numExpr();
					for (int k = 0; k < numBay; k++) {
						expr = upperPro.sum(expr, upperPro.prod(Xhki[k][i][j], 1));
					}
					upper1_3[i][j] = upperPro.addLe(expr, 1, "constraint1_3 " + i + j);
				}
			}

			// 约束1-4
			for (int i = 0; i < numQC; i++) {
				for (int j = 1; j < numBay; j++) { // 注意这条约束i是从2开始的
					IloNumExpr expr = upperPro.numExpr();
					for (int k = 0; k < numBay; k++) {
						expr = upperPro.sum(expr, upperPro.prod(Xhki[k][i][j], 1));
						expr = upperPro.sum(expr, upperPro.prod(Xhki[k][i][j - 1], -1));
					}
					upper1_4[i][j] = upperPro.addLe(expr, 0, "constraint1_4 " + i + j);
				}
			}

			if (upperPro.solve()) {
				System.out.println("求解成功了");
				XhkiValue = new double[numBay][numQC][numBay];
				for (int i = 0; i < numBay; i++) {
					for (int j = 0; j < numQC; j++) {
						XhkiValue[i][j] = upperPro.getValues(Xhki[i][j]);
					}
				}

				System.out.println("Solution status = " + upperPro.getStatus());
				System.out.println("Solution value  = " + upperPro.getObjValue());

				for (int i = 0; i < numBay; i++) {
					for (int j = 0; j < numQC; j++) {
						for (int j2 = 0; j2 < numBay; j2++) {
							System.out.println("X" + i + j + j2 + ": Value = " + XhkiValue[i][j][j2]);
						}
					}
				}

			} else {
				System.out.println("上层问题未成功求解");
			}
		} catch (IloException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	public void lowerModel(Input input) {
		try {
			int numBay = input.numBay;
			int numQC = input.numQC;
			lowerPro = new IloCplex();
			// 变量
			beta_h = new IloNumVar[numBay];
			lamda_h = new IloNumVar[numBay];
			gamma_hki = new IloNumVar[numBay][numQC][numBay];
			Ch = new IloNumVar[numBay];
			zh = new IloNumVar[numBay];
			miu_h = new IloNumVar[numBay];
			niu_hki = new IloNumVar[numBay][numQC][numBay];
			ksi_h = new IloNumVar[numBay];
			oh = new IloNumVar[numBay];
			pai = new IloNumVar[1];
			Cmax = new IloNumVar[1];
			// 约束
			lower2_2 = new IloRange[numBay];
			lower2_3 = new IloRange[numBay][numQC][numBay];
			lower2_4 = new IloRange[numBay];
			lower2_5 = new IloRange[numBay];
			lower2_6 = new IloRange[1];
			lower2_7 = new IloRange[numBay];
			lower2_8 = new IloRange[numBay];
			lower2_9 = new IloRange[numBay][numQC][numBay];
			lower2_10 = new IloRange[numBay][numQC][numBay];
			lower2_11 = new IloRange[numBay];
			lower2_12 = new IloRange[numBay];
			lower2_13 = new IloRange[numBay];
			lower2_14 = new IloRange[numBay];
			lower2_15 = new IloRange[1];
			lower2_16 = new IloRange[1];
			lower2_17 = new IloRange[1];
			// 设置变量
			for (int i = 0; i < numBay; i++) {
				beta_h[i] = lowerPro.numVar(0, Double.MAX_VALUE, IloNumVarType.Float, "beta_" + i);
			}
			for (int i = 0; i < numBay; i++) {
				lamda_h[i] = lowerPro.intVar(0, Integer.MAX_VALUE, "lamda_" + i);
//				lamda_h[i] = lowerPro.numVar(0, Double.MAX_VALUE, IloNumVarType.Float, "lamda_" + i);
			}
			for (int i = 0; i < numBay; i++) {
				for (int j = 0; j < numQC; j++) {
					for (int j2 = 1; j2 < numBay; j2++) { // 实际索引从1开始
						gamma_hki[i][j][j2] = lowerPro.numVar(0, Double.MAX_VALUE, IloNumVarType.Float,
								"gamma_" + i + j + j2);
					}
				}
			}
			for (int i = 0; i < numBay; i++) {
				Ch[i] = lowerPro.numVar(0, Double.MAX_VALUE, IloNumVarType.Float, "C_" + i);
			}
			for (int i = 0; i < numBay; i++) {
				zh[i] = lowerPro.numVar(0, 1, IloNumVarType.Float, "z_" + i);
			}
			for (int i = 0; i < numBay; i++) {
				miu_h[i] = lowerPro.intVar(0, 1, "miu_" + i);
//				miu_h[i] = lowerPro.numVar(0, 1, IloNumVarType.Int, "miu_" + i);
			}
			for (int i = 0; i < numBay; i++) {
				for (int j = 0; j < numQC; j++) {
					for (int j2 = 0; j2 < numBay; j2++) {
						niu_hki[i][j][j2] = lowerPro.intVar(0, 1, "niu_" + i + j + j2);
//						niu_hki[i][j][j2] = lowerPro.numVar(0, 1, IloNumVarType.Int, "niu_" + i + j + j2);
					}
				}
			}
			for (int i = 0; i < numBay; i++) {
				ksi_h[i] = lowerPro.intVar(0, 1, "ksi_" + i);
//				ksi_h[i] = lowerPro.numVar(0, 1, IloNumVarType.Int, "ksi_" + i);
			}
			for (int i = 0; i < numBay; i++) {
				oh[i] = lowerPro.intVar(0, 1, "o" + i);
//				oh[i] = lowerPro.numVar(0, 1, IloNumVarType.Int, "o" + i);
			}
			pai[0] = lowerPro.numVar(0, 1, IloNumVarType.Int, "pai");
			Cmax[0] = lowerPro.numVar(Double.MIN_VALUE, Double.MAX_VALUE, IloNumVarType.Float, "Cmax");

			// 设置约束
			// 2-2
			for (int i = 0; i < numBay; i++) {
				IloNumExpr expr = lowerPro.numExpr();
				double p1 = input.processTime[i]; // 稳定处理时间
				double p2 = p1 * input.uncertain; // 波动处理时间
				for (int j = 0; j < numQC; j++) {
					double t = Math.abs((input.initLocation[j] - 1 - i) * input.traverseTime); // 与实际索引差1
					expr = lowerPro.sum(expr, (p1 + t) * XhkiValue[i][j][0]);
					expr = lowerPro.sum(expr, lowerPro.prod(XhkiValue[i][j][0] * p2, zh[i]));
				}
				expr = lowerPro.sum(expr, lowerPro.prod(-1, Ch[i]));
				lower2_2[i] = lowerPro.addLe(expr, 0, "Constraint2_2 " + i);
			}
			// 2-3
			for (int i = 0; i < numBay; i++) {
				for (int j = 0; j < numQC; j++) {
					for (int j2 = 1; j2 < numBay; j2++) {
						IloNumExpr expr = lowerPro.numExpr();
//						IloNumExpr expr1 = lowerPro.numExpr();
						double p1 = input.processTime[i];
						double p2 = p1 * input.uncertain;
						for (int k = 0; k < numBay; k++) {
							double t = Math.abs(k - i) * input.traverseTime;
							expr = lowerPro.sum(expr, lowerPro.prod(XhkiValue[k][j][j2 - 1], Ch[k]));
							expr = lowerPro.sum(expr, t * XhkiValue[k][j][j2 - 1]);
						}
						expr = lowerPro.sum(expr, p1);
						expr = lowerPro.sum(expr, lowerPro.prod(p2, zh[i]));
						expr = lowerPro.sum(expr, lowerPro.prod(-1, Ch[i]));
						expr = lowerPro.sum(expr, -M * (1 - XhkiValue[i][j][j2]));
						lower2_3[i][j][j2] = lowerPro.addLe(expr, 0, "Constraint2_3 " + i + j + j2);
					}
				}
			}
			// 2-4
			IloNumExpr expr = lowerPro.numExpr();
			for (int i = 0; i < numBay; i++) {
				expr = lowerPro.numExpr();
				expr = lowerPro.sum(expr, lowerPro.prod(1, Ch[i]));
				expr = lowerPro.sum(expr, lowerPro.prod(-1, Cmax[0]));
				lower2_4[i] = lowerPro.addLe(expr, 0, "Constraint2_4 " + i);
			}
			// 2-5
			for (int i = 0; i < numBay; i++) {
				expr = lowerPro.numExpr();
				expr = lowerPro.sum(expr, lowerPro.prod(-1, beta_h[i]));
				expr = lowerPro.sum(expr, lowerPro.prod(1, lamda_h[i]));
				for (int j = 0; j < numQC; j++) {
					for (int j2 = 1; j2 < numBay; j2++) { // 索引实际从2开始的
						expr = lowerPro.sum(expr, lowerPro.prod(1, gamma_hki[i][j][j2]));
					}
				}
				for (int j = 0; j < numQC; j++) {
					for (int j2 = 1; j2 < numBay; j2++) { // 索引实际从2开始的
						for (int k = 0; k < numBay; k++) {
							expr = lowerPro.sum(expr, lowerPro.prod(-XhkiValue[i][j][j2 - 1], gamma_hki[k][j][j2]));
						}
					}
				}
				lower2_5[i] = lowerPro.addLe(expr, 0, "Constraint2_5 " + i);
			}
			// 2-6
			expr = lowerPro.numExpr();
			for (int i = 0; i < Ch.length; i++) {
				expr = lowerPro.sum(expr, lowerPro.prod(1, beta_h[i]));
			}
			lower2_6[0] = lowerPro.addLe(expr, 1, "Constraint2_6");
			// 2-7
			for (int i = 0; i < numBay; i++) {
				expr = lowerPro.numExpr();
				double p1 = input.processTime[i]; // 稳定处理时间
				double p2 = p1 * input.uncertain; // 波动处理时间
				for (int j = 0; j < numQC; j++) {
					double t = Math.abs(input.initLocation[j] - 1 - i) * input.traverseTime;
					expr = lowerPro.sum(expr, (p1 + t) * XhkiValue[i][j][0]);
					expr = lowerPro.sum(expr, lowerPro.prod(XhkiValue[i][j][0] * p2, zh[i]));
				}
				expr = lowerPro.sum(expr, lowerPro.prod(-1, Ch[i]));
				expr = lowerPro.sum(expr, lowerPro.prod(M, miu_h[i]));
				lower2_7[i] = lowerPro.addGe(expr, 0, "Constraint2_7" + i);
			}
			// 2-8
			for (int i = 0; i < numBay; i++) {
				expr = lowerPro.numExpr();
				expr = lowerPro.sum(expr, lowerPro.prod(1, lamda_h[i]));
				expr = lowerPro.sum(expr, lowerPro.prod(M, miu_h[i]));
				lower2_8[i] = lowerPro.addLe(expr, M, "Constraint2_8" + i);
			}
			// 2-9
			for (int i = 0; i < numBay; i++) {
				for (int j = 0; j < numQC; j++) {
					for (int j2 = 1; j2 < numBay; j2++) {
						expr = lowerPro.numExpr();
//						IloNumExpr expr1 = lowerPro.numExpr();
						double p1 = input.processTime[i];
						double p2 = p1 * input.uncertain;
						for (int k = 0; k < numBay; k++) {
							double t = Math.abs(k - i) * input.traverseTime;
							expr = lowerPro.sum(expr, lowerPro.prod(XhkiValue[k][j][j2 - 1], Ch[k]));
							expr = lowerPro.sum(expr, t * XhkiValue[k][j][j2 - 1]);
						}
						expr = lowerPro.sum(expr, p1);
						expr = lowerPro.sum(expr, lowerPro.prod(p2, zh[i]));
						expr = lowerPro.sum(expr, lowerPro.prod(-1, Ch[i]));
						expr = lowerPro.sum(expr, -M * (1 - XhkiValue[i][j][j2]));
						expr = lowerPro.sum(expr, lowerPro.prod(M1, niu_hki[i][j][j2]));
						lower2_9[i][j][j2] = lowerPro.addGe(expr, 0, "Constraint2_9 " + i + j + j2);
					}
				}
			}
			// 2-10
			for (int i = 0; i < numBay; i++) {
				for (int j = 0; j < numQC; j++) {
					for (int j2 = 1; j2 < numBay; j2++) { // 实际索引从2开始
						expr = lowerPro.numExpr();
						expr = lowerPro.sum(expr, gamma_hki[i][j][j2]);
						expr = lowerPro.sum(expr, lowerPro.prod(M, niu_hki[i][j][j2]));
						lower2_10[i][j][j2] = lowerPro.addLe(expr, M, "Constraint2_10 " + i + j + j2);
					}
				}
			}
			// 2-11
			for (int i = 0; i < numBay; i++) {
				expr = lowerPro.numExpr();
				expr = lowerPro.sum(expr, lowerPro.prod(-1, Ch[i]));
				expr = lowerPro.sum(expr, Cmax[0]);
				expr = lowerPro.sum(expr, lowerPro.prod(-M, ksi_h[i]));
				lower2_11[i] = lowerPro.addLe(expr, 0, "Constraint2_11 " + i);
			}
			// 2-12
			for (int i = 0; i < numBay; i++) {
				expr = lowerPro.numExpr();
				expr = lowerPro.sum(expr, beta_h[i]);
				expr = lowerPro.sum(expr, lowerPro.prod(M, ksi_h[i]));
				lower2_12[i] = lowerPro.addLe(expr, M, "Constraint2_12 " + i);
			}
			// 2-13
			for (int i = 0; i < numBay; i++) {
				expr = lowerPro.numExpr();
				expr = lowerPro.sum(expr, lowerPro.prod(-1, beta_h[i]));
				expr = lowerPro.sum(expr, lowerPro.prod(1, lamda_h[i]));
				for (int j = 0; j < numQC; j++) {
					for (int j2 = 1; j2 < numBay; j2++) { // 索引实际从2开始的
						expr = lowerPro.sum(expr, lowerPro.prod(1, gamma_hki[i][j][j2]));
					}
				}
				for (int j = 0; j < numQC; j++) {
					for (int j2 = 1; j2 < numBay; j2++) { // 索引实际从2开始的
						for (int k = 0; k < numBay; k++) {
							expr = lowerPro.sum(expr, lowerPro.prod(-XhkiValue[i][j][j2 - 1], gamma_hki[k][j][j2]));
						}
					}
				}
				expr = lowerPro.sum(expr, lowerPro.prod(M, oh[i]));
				lower2_13[i] = lowerPro.addGe(expr, 0, "Constraint2_13 " + i);
			}
			// 2-14
			for (int i = 0; i < numBay; i++) {
				expr = lowerPro.numExpr();
				expr = lowerPro.sum(expr, Ch[i]);
				expr = lowerPro.sum(expr, lowerPro.prod(M, oh[i]));
				lower2_14[i] = lowerPro.addLe(expr, M, "Constraint2_14 " + i);
			}
			// 2-15
			expr = lowerPro.numExpr();
			for (int i = 0; i < Ch.length; i++) {
				expr = lowerPro.sum(expr, lowerPro.prod(1, beta_h[i]));
			}
			expr = lowerPro.sum(expr, lowerPro.prod(M, pai[0]));
			lower2_15[0] = lowerPro.addGe(expr, 1, "Constraint2_15");
			// 2-16
			expr = lowerPro.numExpr();
			expr = lowerPro.sum(expr, Cmax[0]);
			expr = lowerPro.sum(expr, lowerPro.prod(M, pai[0]));
			lower2_16[0] = lowerPro.addLe(expr, M, "Constraint2_16");
			// 2-17
			expr = lowerPro.numExpr();
			for (int i = 0; i < numBay; i++) {
				expr = lowerPro.sum(expr, zh[i]);
			}
			lower2_17[0] = lowerPro.addLe(expr, input.budget, "Constraint2_17");

			// 设置目标函数
			// 原问题目标函数
			IloNumExpr obj = lowerPro.numExpr();
			obj = lowerPro.sum(obj, Cmax[0]);
			lowerProObj = lowerPro.addMaximize(obj);
//			// 对偶问题目标函数
//			IloNumExpr obj = lowerPro.numExpr();
//			for (int i = 0; i < numBay; i++) {
//				for (int j = 0; j < numQC; j++) {
//					double p1 = input.processTime[i]; // 稳定处理时间
//					double p2 = p1 * input.uncertain; // 波动处理时间
//					double t = Math.abs(input.initLocation[j] - 1 - i) * input.traverseTime;
//					obj = lowerPro.sum(obj, lowerPro.prod((p1 + t) * XhkiValue[i][j][0], lamda_h[i]));
//					for (int k = 1; k < numBay; k++) {
//						double constant = p1;
//						for (int l = 0; l < numBay; l++) {
//							double t1 = Math.abs(l - i) * input.traverseTime;
//							constant += t1 * XhkiValue[l][j][k - 1];
//						}
//						constant -= M * (1 - XhkiValue[i][j][k]);
//						obj = lowerPro.sum(obj, lowerPro.prod(constant, gamma_hki[i][j][k]));
//					}
//					
//				}
//			}
//			lowerPro.addMaximize(obj);

			if (lowerPro.solve()) {
				double[] beta_hValue = lowerPro.getValues(beta_h);
				double[] lamda_hValue = lowerPro.getValues(lamda_h);
				double[][][] gamma_hkiValue = new double[numBay][numQC][numBay];
				for (int i = 0; i < gamma_hkiValue.length; i++) {
					for (int j = 0; j < gamma_hkiValue[0].length; j++) {
						for (int j2 = 1; j2 < gamma_hkiValue[0][0].length; j2++) { // 实际索引从2开始
							gamma_hkiValue[i][j][j2] = lowerPro.getValue(gamma_hki[i][j][j2]);
						}
					}
				}
				double[] ChValue = lowerPro.getValues(Ch);
				double[] zhValue = lowerPro.getValues(zh);
				double[] miu_hValue = lowerPro.getValues(miu_h);
				double[][][] niu_hkiValue = new double[numBay][numQC][numBay];
				for (int i = 0; i < niu_hkiValue.length; i++) {
					for (int j = 0; j < niu_hkiValue[0].length; j++) {
						for (int j2 = 1; j2 < niu_hkiValue[0][0].length; j2++) { // 实际索引从2开始
							niu_hkiValue[i][j][j2] = lowerPro.getValue(niu_hki[i][j][j2]);
						}
					}
				}
				double[] ksi_hValue = lowerPro.getValues(ksi_h);
				double[] ohValue = lowerPro.getValues(oh);
				double[] paiValue = lowerPro.getValues(pai);
				double[] CmaxValue = lowerPro.getValues(Cmax);

				System.out.println("Solution status = " + lowerPro.getStatus());
				System.out.println("Solution value  = " + lowerPro.getObjValue());

				for (int i = 0; i < beta_hValue.length; i++) {
					System.out.println("beta_" + i + ": Value = " + beta_hValue[i]);
				}
				for (int i = 0; i < lamda_hValue.length; i++) {
					System.out.println("lamda_" + i + ": Value = " + lamda_hValue[i]);
				}
				for (int i = 0; i < numBay; i++) {
					for (int j = 0; j < numQC; j++) {
						for (int j2 = 1; j2 < numBay; j2++) { // 索引从2开始
							System.out.println("gamma_" + i + j + j2 + ": Value = " + gamma_hkiValue[i][j][j2]);
						}
					}
				}
				for (int i = 0; i < ChValue.length; i++) {
					System.out.println("C" + i + ": Value = " + ChValue[i]);
				}
				for (int i = 0; i < zhValue.length; i++) {
					System.out.println("z" + i + ": Value = " + zhValue[i]);
				}
				for (int i = 0; i < miu_hValue.length; i++) {
					System.out.println("miu_" + i + ": Value = " + miu_hValue[i]);
				}
				for (int i = 0; i < numBay; i++) {
					for (int j = 0; j < numQC; j++) {
						for (int j2 = 1; j2 < numBay; j2++) { // 索引从2开始
							System.out.println("niu_" + i + j + j2 + ": Value = " + niu_hkiValue[i][j][j2]);
						}
					}
				}
				for (int i = 0; i < ksi_hValue.length; i++) {
					System.out.println("ksi_" + i + ": Value = " + ksi_hValue[i]);
				}
				for (int i = 0; i < ohValue.length; i++) {
					System.out.println("o" + i + ": Value = " + ohValue[i]);
				}
				System.out.println("beta: Value = " + paiValue[0]);
				System.out.println("Cmax: Value = " + CmaxValue[0]);
			} else {
				System.out.println("下层问题求解失败");
			}
		} catch (IloException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	public void ccgIteration(Input input) {
		try {
			double lowerBound = upperPro.getObjValue();
			double upperBound = lowerPro.getObjValue();

			while (upperBound - lowerBound >= ACCURANCY) {
				Iteration += 1;
				int numBay = input.numBay;
				int numQC = input.numQC;
				System.out.println("！！！第" + Iteration + "次加入cut！！！");

				// 创建变量
				IloNumVar[] Chl = new IloNumVar[numBay];
				IloNumVar[][][] Dhkil = new IloNumVar[numBay][numQC][numBay - 1];
				// 创建约束
				IloRange[] upper1_5 = new IloRange[numBay]; // 约束1-5
				IloRange[][][] upper1_6 = new IloRange[numBay][numQC][numBay]; // 约束1-6
				IloRange[][][] upper1_7 = new IloRange[numBay][numQC][numBay - 1]; // 约束1-7
				IloRange[][][] upper1_8 = new IloRange[numBay][numQC][numBay - 1]; // 约束1-8
				IloRange[][][] upper1_9 = new IloRange[numBay][numQC][numBay - 1]; // 约束1-9
				IloRange[] upper1_10 = new IloRange[numBay]; // 约束1-10

				// 设置变量
				for (int i = 0; i < numBay; i++) {
					Chl[i] = upperPro.numVar(0, Double.MAX_VALUE, IloNumVarType.Float, "Ch" + Iteration + "_" + i);
				}
				for (int i = 0; i < numBay; i++) {
					for (int j = 0; j < numQC; j++) {
						for (int j2 = 0; j2 < numBay - 1; j2++) {
							Dhkil[i][j][j2] = upperPro.numVar(0, Double.MAX_VALUE, IloNumVarType.Float,
									"Dhki" + Iteration + "_" + i);
						}
					}
				}

				// 设置约束
				// 1-5
				double[] zhValue = lowerPro.getValues(zh);
				for (int i = 0; i < numBay; i++) {
					IloNumExpr expr = upperPro.numExpr();
					for (int j = 0; j < numQC; j++) {
						double p1 = input.processTime[i];
						double p2 = p1 * input.uncertain;
						double t = Math.abs(input.initLocation[j] - 1 - i) * input.traverseTime;
						expr = upperPro.sum(expr, upperPro.prod(p1 + p2 * zhValue[i] + t, Xhki[i][j][0]));
					}
					expr = upperPro.sum(expr, upperPro.prod(-1, Chl[i]));
					upper1_5[i] = upperPro.addLe(expr, 0, "Constraint1_5 " + Iteration + i);
				}
				// 1-6
				for (int i = 0; i < numBay; i++) {
					for (int j = 0; j < numQC; j++) {
						for (int j2 = 1; j2 < numBay; j2++) {
							IloNumExpr expr = upperPro.numExpr();
							for (int k = 0; k < numBay; k++) {
								expr = upperPro.sum(expr, Dhkil[k][j][j2 - 1]);
								double t = Math.abs(i - k) * input.traverseTime;
								expr = upperPro.sum(expr, upperPro.prod(t, Xhki[k][j][j2 - 1]));
							}
							double p1 = input.processTime[i];
							double p2 = p1 * input.uncertain;
							expr = upperPro.sum(expr, p1 + p2 * zhValue[i]);
							expr = upperPro.sum(expr, upperPro.prod(-1, Chl[i]));
							expr = upperPro.sum(expr, upperPro.prod(M, Xhki[i][j][j2]));
							upper1_6[i][j][j2] = upperPro.addLe(expr, M, "Constraint1_6 " + Iteration + i + j + j2);
						}
					}
				}
				// 1-7
				for (int i = 0; i < numBay; i++) {
					for (int j = 0; j < numQC; j++) {
						for (int j2 = 0; j2 < numBay - 1; j2++) {
							IloNumExpr expr = upperPro.numExpr();
							expr = upperPro.sum(expr, Dhkil[i][j][j2]);
							expr = upperPro.sum(expr, upperPro.prod(-1, Chl[i]));
							upper1_7[i][j][j2] = upperPro.addLe(expr, 0, "Constraint1_7 " + Iteration + i + j + j2);
						}
					}
				}
				// 1-8
				for (int i = 0; i < numBay; i++) {
					for (int j = 0; j < numQC; j++) {
						for (int j2 = 0; j2 < numBay - 1; j2++) {
							IloNumExpr expr = upperPro.numExpr();
							expr = upperPro.sum(expr, Dhkil[i][j][j2]);
							expr = upperPro.sum(expr, upperPro.prod(-M, Xhki[i][j][j2]));
							upper1_8[i][j][j2] = upperPro.addLe(expr, 0, "Constraint1_8 " + Iteration + i + j + j2);
						}
					}
				}
				// 1-9
				for (int i = 0; i < numBay; i++) {
					for (int j = 0; j < numQC; j++) {
						for (int j2 = 0; j2 < numBay - 1; j2++) {
							IloNumExpr expr = upperPro.numExpr();
							expr = upperPro.sum(expr, Chl[i]);
							expr = upperPro.sum(expr, upperPro.prod(M, Xhki[i][j][j2]));
							expr = upperPro.sum(expr, upperPro.prod(-1, Dhkil[i][j][j2]));
							upper1_9[i][j][j2] = upperPro.addLe(expr, M, "Constraint1_9 " + Iteration + i + j + j2);
						}
					}
				}
				// 1-10
				for (int i = 0; i < numBay; i++) {
					IloNumExpr expr = upperPro.numExpr();
					expr = upperPro.sum(expr, Chl[i]);
					expr = upperPro.sum(expr, upperPro.prod(-1, Cmax1[0]));
					upper1_10[i] = upperPro.addLe(expr, 0, "Constraint1_10 " + Iteration + i);
				}

				// 对新的上层问题进行求解
				if (upperPro.solve()) {
					System.out.println("上层问题第" + Iteration + "次求解成功了");
					lowerBound = upperPro.getValue(Cmax1[0]);
					System.out.println("下界 = " + lowerBound);
					// 更新X值
					for (int i = 0; i < XhkiValue.length; i++) {
						for (int j = 0; j < XhkiValue[0].length; j++) {
							XhkiValue[i][j] = upperPro.getValues(Xhki[i][j]);
						}
					}
					for (int i = 0; i < numBay; i++) {
						for (int j = 0; j < numQC; j++) {
							for (int j2 = 0; j2 < numBay; j2++) {
								System.out.println("X" + i + j + j2 + ": Value = " + XhkiValue[i][j][j2]);
							}
						}
					}
				} else {
					System.out.println("上层问题未成功求解");
				}

				// 对子问题的部分约束进行更新
				// 2-2
				for (int i = 0; i < numBay; i++) {
					IloNumExpr expr = lowerPro.numExpr();
					double p1 = input.processTime[i]; // 稳定处理时间
					double p2 = p1 * input.uncertain; // 波动处理时间
					for (int j = 0; j < numQC; j++) {
						double t = Math.abs((input.initLocation[j] - 1 - i) * input.traverseTime); // 与实际索引差1
						expr = lowerPro.sum(expr, (p1 + t) * XhkiValue[i][j][0]);
						expr = lowerPro.sum(expr, lowerPro.prod(XhkiValue[i][j][0] * p2, zh[i]));
					}
					expr = lowerPro.sum(expr, lowerPro.prod(-1, Ch[i]));
					lowerPro.remove(lower2_2[i]);
					lower2_2[i] = lowerPro.addLe(expr, 0, "Constraint2_2 " + i);
				}
				// 2-3
				for (int i = 0; i < numBay; i++) {
					for (int j = 0; j < numQC; j++) {
						for (int j2 = 1; j2 < numBay; j2++) {
							IloNumExpr expr = lowerPro.numExpr();
//							IloNumExpr expr1 = lowerPro.numExpr();
							double p1 = input.processTime[i];
							double p2 = p1 * input.uncertain;
							for (int k = 0; k < numBay; k++) {
								double t = Math.abs(k - i) * input.traverseTime;
								expr = lowerPro.sum(expr, lowerPro.prod(XhkiValue[k][j][j2 - 1], Ch[k]));
								expr = lowerPro.sum(expr, t * XhkiValue[k][j][j2 - 1]);
							}
							expr = lowerPro.sum(expr, p1);
							expr = lowerPro.sum(expr, lowerPro.prod(p2, zh[i]));
							expr = lowerPro.sum(expr, lowerPro.prod(-1, Ch[i]));
							expr = lowerPro.sum(expr, -M * (1 - XhkiValue[i][j][j2]));
							lowerPro.remove(lower2_3[i][j][j2]);
							lower2_3[i][j][j2] = lowerPro.addLe(expr, 0, "Constraint2_3 " + i + j + j2);
						}
					}
				}
				// 2-5
				IloNumExpr expr = lowerPro.numExpr();
				for (int i = 0; i < numBay; i++) {
					expr = lowerPro.numExpr();
					expr = lowerPro.sum(expr, lowerPro.prod(-1, beta_h[i]));
					expr = lowerPro.sum(expr, lowerPro.prod(1, lamda_h[i]));
					for (int j = 0; j < numQC; j++) {
						for (int j2 = 1; j2 < numBay; j2++) { // 索引实际从2开始的
							expr = lowerPro.sum(expr, lowerPro.prod(1, gamma_hki[i][j][j2]));
						}
					}
					for (int j = 0; j < numQC; j++) {
						for (int j2 = 1; j2 < numBay; j2++) { // 索引实际从2开始的
							for (int k = 0; k < numBay; k++) {
								expr = lowerPro.sum(expr, lowerPro.prod(-XhkiValue[i][j][j2 - 1], gamma_hki[k][j][j2]));
							}
						}
					}
					lowerPro.remove(lower2_5[i]);
					lower2_5[i] = lowerPro.addLe(expr, 0, "Constraint2_5 " + i);
				}
				// 2-7
				for (int i = 0; i < numBay; i++) {
					expr = lowerPro.numExpr();
					double p1 = input.processTime[i]; // 稳定处理时间
					double p2 = p1 * input.uncertain; // 波动处理时间
					for (int j = 0; j < numQC; j++) {
						double t = Math.abs(input.initLocation[j] - 1 - i) * input.traverseTime;
						expr = lowerPro.sum(expr, (p1 + t) * XhkiValue[i][j][0]);
						expr = lowerPro.sum(expr, lowerPro.prod(XhkiValue[i][j][0] * p2, zh[i]));
					}
					expr = lowerPro.sum(expr, lowerPro.prod(-1, Ch[i]));
					expr = lowerPro.sum(expr, lowerPro.prod(M, miu_h[i]));
					lowerPro.remove(lower2_7[i]);
					lower2_7[i] = lowerPro.addGe(expr, 0, "Constraint2_7" + i);
				}
				// 2-9
				for (int i = 0; i < numBay; i++) {
					for (int j = 0; j < numQC; j++) {
						for (int j2 = 1; j2 < numBay; j2++) {
							expr = lowerPro.numExpr();
//							IloNumExpr expr1 = lowerPro.numExpr();
							double p1 = input.processTime[i];
							double p2 = p1 * input.uncertain;
							for (int k = 0; k < numBay; k++) {
								double t = Math.abs(k - i) * input.traverseTime;
								expr = lowerPro.sum(expr, lowerPro.prod(XhkiValue[k][j][j2 - 1], Ch[k]));
								expr = lowerPro.sum(expr, t * XhkiValue[k][j][j2 - 1]);
							}
							expr = lowerPro.sum(expr, p1);
							expr = lowerPro.sum(expr, lowerPro.prod(p2, zh[i]));
							expr = lowerPro.sum(expr, lowerPro.prod(-1, Ch[i]));
							expr = lowerPro.sum(expr, -M * (1 - XhkiValue[i][j][j2]));
							expr = lowerPro.sum(expr, lowerPro.prod(M1, niu_hki[i][j][j2]));
							lowerPro.remove(lower2_9[i][j][j2]);
							lower2_9[i][j][j2] = lowerPro.addGe(expr, 0, "Constraint2_9 " + i + j + j2);
						}
					}
				}
				// 2-13
				for (int i = 0; i < numBay; i++) {
					expr = lowerPro.numExpr();
					expr = lowerPro.sum(expr, lowerPro.prod(-1, beta_h[i]));
					expr = lowerPro.sum(expr, lowerPro.prod(1, lamda_h[i]));
					for (int j = 0; j < numQC; j++) {
						for (int j2 = 1; j2 < numBay; j2++) { // 索引实际从2开始的
							expr = lowerPro.sum(expr, lowerPro.prod(1, gamma_hki[i][j][j2]));
						}
					}
					for (int j = 0; j < numQC; j++) {
						for (int j2 = 1; j2 < numBay; j2++) { // 索引实际从2开始的
							for (int k = 0; k < numBay; k++) {
								expr = lowerPro.sum(expr, lowerPro.prod(-XhkiValue[i][j][j2 - 1], gamma_hki[k][j][j2]));
							}
						}
					}
					expr = lowerPro.sum(expr, lowerPro.prod(M, oh[i]));
					lowerPro.remove(lower2_13[i]);
					lower2_13[i] = lowerPro.addGe(expr, 0, "Constraint2_13 " + i);
				}

				// 对新的下层问题进行求解
				if (lowerPro.solve()) {
					System.out.println("下层问题第" + Iteration + "次求解成功了");
					upperBound = lowerPro.getValue(Cmax[0]);
					System.out.println("上界 = " + upperBound);
				} else {
					System.out.println("下层问题未成功求解");
				}
//				lowerModel(numBay, numQC, input);
//				upperBound = lowerPro.getObjValue();
			}
		} catch (IloException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

	public void crossForbidCut(Input input) {

		try {
			int numBay = input.numBay;
			int numQC = input.numQC;
			IloCplex finalPro = new IloCplex();
			IloIntVar[][] Yhh = new IloIntVar[numBay][numBay];
			IloNumVar[] Ch1 = new IloNumVar[numBay];
			IloNumVar[] Cmax2 = new IloNumVar[1];
			IloRange[] final3_2 = new IloRange[numBay];
			IloRange[][][] final3_3 = new IloRange[numBay][numQC][numBay];
			IloRange[] final3_4 = new IloRange[numBay];
			IloRange[][] final3_5 = new IloRange[numBay][numBay]; // 对应Y变量的第一条约束
			IloRange[][] final3_6 = new IloRange[numBay][numBay]; // 对应Y变量的第二条约束
			IloRange[][] final3_7 = new IloRange[numBay][numBay]; // 对应Y变量的第三条约束
			// 设置变量
			for (int i = 0; i < numBay; i++) {
				for (int j = 0; j < numBay; j++) {
					Yhh[i][j] = finalPro.intVar(0, 1, "Y" + i + j);
				}
			}
			for (int i = 0; i < numBay; i++) {
				Ch1[i] = finalPro.numVar(0, Double.MAX_VALUE, IloNumVarType.Float, "C_" + i);
			}
			Cmax2[0] = finalPro.numVar(0, Double.MAX_VALUE, IloNumVarType.Float, "Cmax");

			// 设置约束
			// 获取z值
			double[] zhValue = lowerPro.getValues(zh);
			// 3-2
			for (int i = 0; i < numBay; i++) {
				IloNumExpr expr = finalPro.numExpr();
				double p1 = input.processTime[i]; // 稳定处理时间
				double p2 = p1 * input.uncertain; // 波动处理时间
				for (int j = 0; j < numQC; j++) {
					double t = Math.abs((input.initLocation[j] - 1 - i) * input.traverseTime); // 与实际索引差1
					expr = finalPro.sum(expr, (p1 + t) * XhkiValue[i][j][0] + XhkiValue[i][j][0] * p2 * zhValue[j]);
				}
				expr = finalPro.sum(expr, finalPro.prod(-1, Ch1[i]));
				final3_2[i] = finalPro.addLe(expr, 0, "Constraint3_2 " + i);
			}
			// 3-3
			for (int i = 0; i < numBay; i++) {
				for (int j = 0; j < numQC; j++) {
					for (int j2 = 1; j2 < numBay; j2++) {
						IloNumExpr expr = finalPro.numExpr();
//									IloNumExpr expr1 = lowerPro.numExpr();
						double p1 = input.processTime[i];
						double p2 = p1 * input.uncertain;
						for (int k = 0; k < numBay; k++) {
							double t = Math.abs(k - i) * input.traverseTime;
							expr = finalPro.sum(expr, finalPro.prod(XhkiValue[k][j][j2 - 1], Ch1[k]));
							expr = finalPro.sum(expr, t * XhkiValue[k][j][j2 - 1]);
						}
						expr = finalPro.sum(expr, p1 + p2 * zhValue[i]);
						expr = finalPro.sum(expr, finalPro.prod(-1, Ch1[i]));
						expr = finalPro.sum(expr, -M * (1 - XhkiValue[i][j][j2]));
						final3_3[i][j][j2] = finalPro.addLe(expr, 0, "Constraint3_3 " + i + j + j2);
					}
				}
			}
			// 3-4
			IloNumExpr expr = finalPro.numExpr();
			for (int i = 0; i < numBay; i++) {
				expr = finalPro.numExpr();
				expr = finalPro.sum(expr, finalPro.prod(1, Ch1[i]));
				expr = finalPro.sum(expr, finalPro.prod(-1, Cmax2[0]));
				final3_4[i] = finalPro.addLe(expr, 0, "Constraint3_4 " + i);
			}
			// 3-5
			for (int i = 0; i < numBay; i++) {
				for (int j = 0; j < numBay; j++) {
					expr = finalPro.numExpr();
					expr = finalPro.sum(expr, Ch1[i]);
					expr = finalPro.sum(expr, finalPro.prod(-1, Ch1[j]));
					double p1 = input.processTime[j];
					double p2 = p1 * input.uncertain;
					double t = Math.abs(i - j) * input.traverseTime;
					expr = finalPro.sum(expr, p1 + t + p2 * zhValue[j]);
					expr = finalPro.sum(expr, finalPro.prod(M, Yhh[i][j]));
					final3_5[i][j] = finalPro.addGe(expr, 0, "Constraint3_5 " + i + j);
				}
			}
			// 3-6
			for (int i = 0; i < numBay; i++) {
				for (int j = 0; j < numBay; j++) {
					expr = finalPro.numExpr();
					expr = finalPro.sum(expr, Ch1[i]);
					expr = finalPro.sum(expr, finalPro.prod(-1, Ch1[j]));
					double p1 = input.processTime[j];
					double p2 = p1 * input.uncertain;
					double t = Math.abs(i - j) * input.traverseTime;
					expr = finalPro.sum(expr, p1 + t + p2 * zhValue[j] - M);
					expr = finalPro.sum(expr, finalPro.prod(M, Yhh[i][j]));
					final3_6[i][j] = finalPro.addLe(expr, 0, "Constraint3_6 " + i + j);
				}
			}
			// 3-7
			for (int i = 0; i < numBay; i++) {
				for (int j = i + 1; j < numBay; j++) {
					expr = finalPro.numExpr();
					for (int k = 0; k < numQC; k++) {
						for (int k2 = 0; k2 < numBay; k2++) {
							expr = finalPro.sum(expr, k * XhkiValue[i][k][k2] - k * XhkiValue[j][k][k2]);
						}
					}
					expr = finalPro.sum(expr, 1);
					expr = finalPro.sum(expr, finalPro.prod(-M, Yhh[i][j]));
					expr = finalPro.sum(expr, finalPro.prod(-M, Yhh[j][i]));
					final3_7[i][j] = finalPro.addLe(expr, 0, "Constraint3_7 " + i + j);
				}
			}
			
			finalPro.addMinimize(Cmax2[0]);

			// 对新的下层问题进行求解
			if (finalPro.solve()) {
				System.out.println("Solution status = " + finalPro.getStatus());
				System.out.println("Solution value  = " + finalPro.getObjValue());
				double[] ChValue = finalPro.getValues(Ch1);
				double[][] YhhValue = new double[numBay][numBay];
				for (int i = 0; i < YhhValue.length; i++) {
					YhhValue[i] = finalPro.getValues(Yhh[i]);
				}
				for (int i = 0; i < numBay; i++) {
					System.out.println("C" + i + ": Value = " + ChValue[i]);
				}
				for (int i = 0; i < numBay; i++) {
					for (int j = 0; j < numBay; j++) {
						System.out.println("Y" + i + j + ": Value = " + YhhValue[i][j]);
					}
				}
				for (int i = 0; i < numBay; i++) {
					System.out.println("z" + i + ": Value = " + zhValue[i]);
				}
			} else {
				System.out.println("最终问题未成功求解");
			}
		} catch (IloException e) {
			// TODO: handle exception
		}

	}
}
