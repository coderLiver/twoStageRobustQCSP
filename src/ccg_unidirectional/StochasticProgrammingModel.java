package ccg_unidirectional;

import java.io.File;
import java.io.IOException;
import ilog.concert.IloException;
import ilog.concert.IloNumExpr;
import ilog.concert.IloNumVar;
import ilog.concert.IloNumVarType;
import ilog.concert.IloObjective;
import ilog.concert.IloRange;
import ilog.cplex.IloCplex;
import ccg_unidirectional.Model;

public class StochasticProgrammingModel {
	IloCplex upperPro;
	IloCplex[] lowerProList;

	// �ϲ��������
//	IloNumVar[] CmaxKsiUpper;
//	double[] CmaxKsiUpperValue;
	IloNumVar[][] Xhk;
	double[][] XhkValue;
	IloNumVar[] alpha;
	double[] alphaValue;
	IloNumVar[] beta;
	double[] betaValue;
	IloNumVar[] theta;
	double thetaValue;
	// �²��������
	IloNumVar[] CmaxKsiLower;
	double[] CmaxKsiLowerValue;
	IloNumVar[][][] YKsiHK;
	double[][][] YksiHKValue;
	IloNumVar[][][] QKsiJK;
	double[][][] QKsiJKValue;
	IloNumVar[][] RKsiK;
	double[][] RKsiKValue;
	// �²������ż����
	double[][] paiValue;
	// �²�����Ŀ�꺯��
	IloObjective[] lowerProObjList;

	// �ϲ�����Լ��
	IloRange[] upper4_2;
	IloRange[][] upper4_3;
	IloRange[][] upper4_4;
	IloRange[] upper4_8;
	IloRange[] upper4_9;
	// �²�����Լ��
	IloRange[][][] lower4_12;
	IloRange[][] lower4_13;
	// �²������������
	double[][][] matrixT;
	double[][] vectorH;

	long startTime; // �㷨���п�ʼʱ��
	long endTime; // �㷨���н���ʱ��
	private static final double M = 3000;
	private static final double ACCURANCY = 0.1; // C&CG�㷨��ֹ�ľ���
	private static final double FLOATACCURACY = 0.001;
	double terminalSpan = 1800; // �����������ʱ�䣬��sΪ��λ
	int Iteration = 1;
	int jValue; // ��λ����ȥ��ȫ�����ټ�1
	int numScenarioByBay = 2; // ÿһ����λ��ҵʱ�����������
	double[] probabilityKsi; // ��ͬscenario��Ӧ�ĸ���

	public void upperModel(Input input, boolean variableOutput, boolean iterationDetail)
			throws IloException, IOException {
		System.out.println("�ϲ�ģ�Ϳ�ʼ���");
		startTime = System.currentTimeMillis();
		upperPro = new IloCplex();
		upperPro.setOut(null);
		// ����
		Xhk = new IloNumVar[input.numBay][input.numQC];
		alpha = new IloNumVar[input.numQC];
		beta = new IloNumVar[input.numQC];
		theta = new IloNumVar[1];
		// �²������ż����

		// Լ��
		upper4_2 = new IloRange[input.numBay];
		upper4_3 = new IloRange[input.numBay][input.numQC];
		upper4_4 = new IloRange[input.numBay][input.numQC];
		upper4_8 = new IloRange[input.numQC];
		upper4_9 = new IloRange[input.numQC];

		// ���ñ���
		for (int i = 0; i < input.numBay; i++) {
			for (int j = 0; j < input.numQC; j++) {
				Xhk[i][j] = upperPro.numVar(0, 1, IloNumVarType.Int, "X" + i + j);
			}
		}


		for (int i = 0; i < input.numQC; i++) {
			alpha[i] = upperPro.numVar(1, input.numBay, IloNumVarType.Float, "alpha" + i);
		}

		for (int i = 0; i < input.numQC; i++) {
			beta[i] = upperPro.numVar(1, input.numBay, IloNumVarType.Float, "beta" + i);
		}

		theta[0] = upperPro.numVar(Double.MIN_VALUE, Double.MAX_VALUE, IloNumVarType.Float, "Theta");
		// ����Ŀ�꺯��
		IloNumExpr upperObj = upperPro.numExpr();
		upperObj = upperPro.sum(upperObj, upperPro.prod(theta[0], 1));
		upperPro.addMinimize(upperObj);

		// �ϲ�����Լ��
		// Լ��4-2
		for (int i = 0; i < input.numBay; i++) {
			if (input.processTime[i] != 0) {
				IloNumExpr expr = upperPro.numExpr();
				for (int j = 0; j < input.numQC; j++) {
					expr = upperPro.sum(expr, upperPro.prod(Xhk[i][j], 1));
				}
				upper4_2[i] = upperPro.addEq(expr, 1, "constraint4_2 " + i);
			}
		}

		// Լ��4-3
		for (int i = 0; i < input.numBay; i++) {
			for (int j = 0; j < input.numQC; j++) {
				IloNumExpr expr = upperPro.numExpr();
				expr = upperPro.sum(expr, upperPro.prod(Xhk[i][j], M - i - 1));
				expr = upperPro.sum(expr, upperPro.prod(alpha[j], 1));
				upper4_3[i][j] = upperPro.addLe(expr, M, "constraint4_3 " + i + "," + j);
			}
		}

		// Լ��4-4
		for (int i = 0; i < input.numBay; i++) {
			for (int j = 0; j < input.numQC; j++) {
				IloNumExpr expr = upperPro.numExpr();
				expr = upperPro.sum(expr, upperPro.prod(Xhk[i][j], -M - i - 1));
				expr = upperPro.sum(expr, upperPro.prod(beta[j], 1));
				upper4_4[i][j] = upperPro.addGe(expr, -M, "constraint4_4 " + i + "," + j);
			}
		}

		// Լ��4-8
		for (int i = 0; i < input.numQC - 1; i++) {
			IloNumExpr expr = upperPro.numExpr();
			expr = upperPro.sum(expr, upperPro.prod(alpha[i], 1));
			expr = upperPro.sum(expr, upperPro.prod(alpha[i + 1], -1));
			upper4_8[i] = upperPro.addLe(expr, -input.safetyMargin - 1, "constraint4_8 " + i);
		}

		// Լ��4-9
		for (int i = 0; i < input.numQC - 1; i++) {
			IloNumExpr expr = upperPro.numExpr();
			expr = upperPro.sum(expr, upperPro.prod(beta[i], 1));
			expr = upperPro.sum(expr, upperPro.prod(beta[i + 1], -1));
			upper4_9[i] = upperPro.addLe(expr, -input.safetyMargin - 1, "constraint4_9 " + i);
		}
		
		// �������õ��㷨��ֹʱ����㱾�ε������ʱ��
		double leftTime = terminalSpan - (System.currentTimeMillis() - startTime) / 1000;
		upperPro.setParam(IloCplex.DoubleParam.TiLim, leftTime);
		
		// �������
		if (upperPro.solve()) {
			System.out.println("�������ɹ�");
			XhkValue = new double[input.numBay][input.numQC];
			alphaValue = new double[input.numQC];
			betaValue = new double[input.numQC];
			for (int i = 0; i < input.numBay; i++) {
				for (int j = 0; j < input.numQC; j++) {
					XhkValue[i][j] = upperPro.getValue(Xhk[i][j]);
				}
			}
//			CmaxKsiUpperValue = upperPro.getValues(CmaxKsiUpper);
			alphaValue = upperPro.getValues(alpha);
			betaValue = upperPro.getValues(beta);
			thetaValue = upperPro.getValue(theta[0]);
			System.out.println("Solution status = " + upperPro.getStatus());
			System.out.println("Solution value  = " + upperPro.getObjValue());
//			input.writerOverallStochastic.write("Iteration " + Iteration + " upperObj:" + upperPro.getObjValue() + "\r\n");

			if (iterationDetail) {
				File file = new File("./result of benchmark/stochastic_test_" + input.testID);
				if (!file.exists()) {// ����ļ��в�����
					file.mkdir();// �����ļ���
				}
				upperPro.writeSolution("./result of benchmark/stochastic_test_" + input.testID + "/test_" + input.testID
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
			System.out.println("�ϲ�����δ�ɹ����");
			System.exit(0);
		}
	}

	public void lowerModelLShape(Input input, boolean variableOutput, boolean iterationDetail) throws IloException {
		int num_ksi = input.num_ksi;
		lowerProList = new IloCplex[num_ksi];
		jValue = input.numBay - input.safetyMargin - 1;
		probabilityKsi = new double[num_ksi];
		for (int i = 0; i < num_ksi; i++) {
			probabilityKsi[i] = 1 / (double) num_ksi;
		}
		// ����
		CmaxKsiLower = new IloNumVar[num_ksi];
		CmaxKsiLowerValue = new double[num_ksi];
		YKsiHK = new IloNumVar[num_ksi][input.numBay][input.numQC];
		QKsiJK = new IloNumVar[num_ksi][jValue][input.numQC - 1];
		RKsiK = new IloNumVar[num_ksi][input.numQC];
		YksiHKValue = new double[num_ksi][input.numBay][input.numQC];
		QKsiJKValue = new double[num_ksi][jValue][input.numQC - 1];
		RKsiKValue = new double[num_ksi][input.numQC];

		// Լ��
		lower4_12 = new IloRange[num_ksi][jValue][input.numQC - 1];
		lower4_13 = new IloRange[num_ksi][input.numQC];

		// ����
		int num_row = jValue * (input.numQC - 1) + input.numQC;
		int num_column = input.numBay * input.numQC + 2 * input.numQC;
		if (matrixT == null) {
			matrixT = new double[num_ksi][num_row][num_column];
		}
		if (vectorH == null) {
			vectorH = new double[num_ksi][num_row];
		}
		
		// ��ż����ȡֵ
		paiValue = new double[num_ksi][num_row];

		// ����ģ�ͽ��н�ģ
		for (int i = 0; i < num_ksi; i++) {
			// ��¼��ǰ��
			int current_row = 0;
			System.out.println("�²�ģ�͵�" + i + "��scenario���");
			lowerProList[i] = new IloCplex();
			lowerProList[i].setOut(null);
			// ���ñ���
			CmaxKsiLower[i] = lowerProList[i].numVar(0, Double.MAX_VALUE, IloNumVarType.Float, "CmaxKsi_" + i);
			for (int j = 0; j < jValue; j++) {
				for (int j2 = 0; j2 < input.numQC - 1; j2++) {
					QKsiJK[i][j][j2] = lowerProList[i].numVar(0, Double.MAX_VALUE, IloNumVarType.Float,
							"QKsiJK_" + i + "," + j + "," + j2);
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

			// Ŀ�꺯��
			IloNumExpr obj = lowerProList[i].numExpr();
			obj = lowerProList[i].sum(obj, CmaxKsiLower[i]);
			lowerProList[i].addMinimize(obj);

			// ����Լ��
			// 4-12
			for (int j = 0; j < jValue; j++) {
				for (int j2 = 0; j2 < input.numQC - 1; j2++) {
					IloNumExpr expr = lowerProList[i].numExpr();
					double left_constant = 0, right_constant = 0;

					for (int k = 0; k <= j; k++) {
						left_constant += input.stochasticProcessTime[i][k] * XhkValue[k][j2];
						matrixT[i][current_row][k * input.numQC + j2] += input.stochasticProcessTime[i][k];
						expr = lowerProList[i].sum(expr, YKsiHK[i][k][j2]);
					}
					expr = lowerProList[i].sum(expr, lowerProList[i].prod(QKsiJK[i][j][j2], -1));
					left_constant += (j + 1 - alphaValue[j2]) * input.traverseTime;
					matrixT[i][current_row][input.numBay * input.numQC + j2] += -input.traverseTime;
					vectorH[i][current_row] += -(j + 1) * input.traverseTime;

					for (int k = input.safetyMargin + 1; k <= j + input.safetyMargin + 1; k++) {
						right_constant += input.stochasticProcessTime[i][k] * XhkValue[k][j2 + 1];
						matrixT[i][current_row][k * input.numQC + j2 + 1] = -input.stochasticProcessTime[i][k];
						expr = lowerProList[i].sum(expr, lowerProList[i].prod(YKsiHK[i][k][j2 + 1], -1));
					}
					
					right_constant += (j + 1 + input.safetyMargin + 1 - alphaValue[j2 + 1]) * input.traverseTime;
					matrixT[i][current_row][input.numBay * input.numQC + j2 + 1] += input.traverseTime;
					vectorH[i][current_row] += (j + 1 + input.safetyMargin + 1) * input.traverseTime;
					current_row++;
					lower4_12[i][j][j2] = lowerProList[i].addEq(expr, right_constant - left_constant,
							"Constraint4_12:" + i + "," + j + "," + j2);
				}
			}
			
			// 4-13
			for (int j = 0; j < input.numQC; j++) {
				double left_constant = 0;
				IloNumExpr expr = lowerProList[i].numExpr();
				for (int k = 0; k < input.numBay; k++) {
					left_constant += input.stochasticProcessTime[i][k] * XhkValue[k][j];
					matrixT[i][current_row][k * input.numQC + j] += input.stochasticProcessTime[i][k];
					expr = lowerProList[i].sum(expr, YKsiHK[i][k][j]);
				}
				left_constant += (betaValue[j] - alphaValue[j]) * input.traverseTime;
				matrixT[i][current_row][input.numBay * input.numQC + j] += -input.traverseTime;
				matrixT[i][current_row][input.numBay * input.numQC + input.numQC + j] += input.traverseTime;
				current_row++;
				expr = lowerProList[i].sum(expr, RKsiK[i][j]);
				expr = lowerProList[i].sum(expr, lowerProList[i].prod(-1, CmaxKsiLower[i]));
				lower4_13[i][j] = lowerProList[i].addLe(expr, -left_constant, "Contraint4_13" + i + "," + j);
			}

			// �������
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

				// ��ȡ��ż����ȡֵ
				// Լ��4_12
				int temp_index = 0;
				for (int j = 0; j < jValue; j++) {
					for (int j2 = 0; j2 < input.numQC - 1; j2++) {
						paiValue[i][temp_index] = lowerProList[i].getDual(lower4_12[i][j][j2]);
						temp_index++;
//						System.out.println("pai_" + i + '_' + (j * (input.numQC - 1) + j2) + ":"
//								+ paiValue[i][j * (input.numQC - 1) + j2]);
					}
				}
				// Լ��4_13
				for (int j = 0; j < input.numQC; j++) {
					paiValue[i][temp_index] = lowerProList[i].getDual(lower4_13[i][j]);
//					System.out.println("pai_" + i + '_' + (jValue * (input.numQC - 1) + j) + ":"
//							+ paiValue[i][jValue * (input.numQC - 1) + j]);
				}		
				
			} else {
				System.out.println("�²��������ʧ��");
				System.exit(0);
			}
		}
		// ��ȡCmax����ֵ
		for (int i = 0; i < num_ksi; i++) {
			CmaxKsiLowerValue[i] = lowerProList[i].getValue(CmaxKsiLower[i]);
		}
		
	}

	public void lowerModelIterate(Input input, boolean variableOutput, boolean iterationDetail) throws IloException, IOException {
		int num_ksi = input.num_ksi;
		System.out.println("�²������" + Iteration + "�ε���");
		// ��������½�ָL-shape�е�w��theta�������²���������½�
		double lowerBound = 0;
		double upperBound = lowerBound + 1;
		int num_row = jValue * (input.numQC - 1) + input.numQC;
		int num_column = input.numBay * input.numQC + 2 * input.numQC;
		double[] EVector = new double[num_column];
		double eScalar = 0;

		for (int i = 0; i < num_ksi; i++) {
			// ����E_{s+1}
			for (int j = 0; j < num_column; j++) {
				for (int j2 = 0; j2 < num_row; j2++) {
					EVector[j] += probabilityKsi[i] * paiValue[i][j2] * matrixT[i][j2][j];
				}
			}

			// ����e_{s+1}
			for (int j2 = 0; j2 < num_row; j2++) {
				eScalar += probabilityKsi[i] * paiValue[i][j2] * vectorH[i][j2];
			}
		}

		upperBound = eScalar;
		lowerBound = thetaValue;
		int temp_index1 = 0;
		for (int i = 0; i < input.numBay; i++) {
			for (int j = 0; j < input.numQC; j++) {
				upperBound -= EVector[temp_index1] * XhkValue[i][j];
				temp_index1++;
			}
		}
		for (int i = 0; i < input.numQC; i++) {
			upperBound -= EVector[temp_index1] * alphaValue[i];
			temp_index1++;
		}
		for (int i = 0; i < input.numQC; i++) {
			upperBound -= EVector[temp_index1] * betaValue[i];
			temp_index1++;
		}
		
		double total_obj = 0;
		
		for (int i = 0; i < num_ksi; i++) {
			total_obj += probabilityKsi[i] * CmaxKsiLowerValue[i];
		}
		input.writerOverallStochastic.write("Iteration " + Iteration + "upperBound:" + upperBound + "\r\n");
		input.writerOverallStochastic.write("Iteration " + Iteration + "lowerBound:" + lowerBound + "\r\n");
		input.writerOverallStochastic.write("\r\n");
		System.out.println("�Ͻ�Ϊ" + upperBound);
		System.out.println("�½�Ϊ" + lowerBound);

		while (upperBound - lowerBound >= ACCURANCY && (endTime - startTime) / 1000 < terminalSpan) {
			Iteration += 1;
			System.out.println("��������" + Iteration + "�μ���cut������");

			// ����Լ��
			IloRange[] lower4_10 = new IloRange[1];

			// ����Լ��
			// 4-10
			IloNumExpr temp_expr = upperPro.numExpr();
			int temp_index = 0;
			for (int i = 0; i < input.numBay; i++) {
				for (int j = 0; j < input.numQC; j++) {
					temp_expr = upperPro.sum(temp_expr, upperPro.prod(Xhk[i][j], EVector[temp_index]));
					temp_index++;
				}
			}

			for (int i = 0; i < input.numQC; i++) {
				temp_expr = upperPro.sum(temp_expr, upperPro.prod(alpha[i], EVector[temp_index]));
				temp_index++;
			}

			for (int i = 0; i < input.numQC; i++) {
				temp_expr = upperPro.sum(temp_expr, upperPro.prod(beta[i], EVector[temp_index]));
				temp_index++;
			}
			
			temp_expr = upperPro.sum(temp_expr, upperPro.prod(theta[0], 1));
			lower4_10[0] = upperPro.addGe(temp_expr, eScalar, "Constraint4_10");
			
			// �������õ��㷨��ֹʱ����㱾�ε������ʱ��
			double leftTime = terminalSpan - (System.currentTimeMillis() - startTime) / 1000;
			upperPro.setParam(IloCplex.DoubleParam.TiLim, leftTime);
			
			// �������
			if (upperPro.solve()) {
				System.out.println("�ϲ��������ɹ�");
				for (int i = 0; i < input.numBay; i++) {
					for (int j = 0; j < input.numQC; j++) {
						XhkValue[i][j] = upperPro.getValue(Xhk[i][j]);
					}
				}
//				CmaxKsiUpperValue = upperPro.getValues(CmaxKsiUpper);
				alphaValue = upperPro.getValues(alpha);
				betaValue = upperPro.getValues(beta);
				thetaValue = upperPro.getValue(theta[0]);
				System.out.println("Solution status = " + upperPro.getStatus());
				System.out.println("Solution value  = " + upperPro.getObjValue());

				if (iterationDetail) {
					File file = new File("./result of benchmark/test_" + input.testID);
					if (!file.exists()) {// ����ļ��в�����
						file.mkdir();// �����ļ���
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
				System.out.println("�ϲ�����δ�ɹ����");
				System.exit(0);
			}
			
			// ���������
			for (int i = 0; i < num_ksi; i++) {
				// ����Լ��
				// 4-12
				for (int j = 0; j < jValue; j++) {
					for (int j2 = 0; j2 < input.numQC - 1; j2++) {
						IloNumExpr expr = lowerProList[i].numExpr();
						double left_constant = 0, right_constant = 0;

						for (int k = 0; k <= j; k++) {
							left_constant += input.stochasticProcessTime[i][k] * XhkValue[k][j2];
							expr = lowerProList[i].sum(expr, YKsiHK[i][k][j2]);
						}
						expr = lowerProList[i].sum(expr, lowerProList[i].prod(QKsiJK[i][j][j2], -1));
						left_constant += (j + 1 - alphaValue[j2]) * input.traverseTime;
						
						
						for (int k = input.safetyMargin + 1; k <= j + input.safetyMargin + 1; k++) {
							right_constant += input.stochasticProcessTime[i][k] * XhkValue[k][j2 + 1];
							expr = lowerProList[i].sum(expr, lowerProList[i].prod(YKsiHK[i][k][j2 + 1], -1));
						}
						right_constant += (j + 1 + input.safetyMargin + 1 - alphaValue[j2 + 1]) * input.traverseTime;
						lowerProList[i].remove(lower4_12[i][j][j2]);
						lower4_12[i][j][j2] = lowerProList[i].addEq(expr, right_constant - left_constant,
								"Constraint4_12:" + i + "," + j + "," + j2);
					}
				}
				// 4-13
				for (int j = 0; j < input.numQC; j++) {
					double left_constant = 0;
					IloNumExpr expr = lowerProList[i].numExpr();
					for (int k = 0; k < input.numBay; k++) {
						left_constant += input.stochasticProcessTime[i][k] * XhkValue[k][j];
						expr = lowerProList[i].sum(expr, YKsiHK[i][k][j]);
					}
					left_constant += (betaValue[j] - alphaValue[j]) * input.traverseTime;
					expr = lowerProList[i].sum(expr, RKsiK[i][j]);
					expr = lowerProList[i].sum(expr, lowerProList[i].prod(-1, CmaxKsiLower[i]));
					lowerProList[i].remove(lower4_13[i][j]);
					lower4_13[i][j] = lowerProList[i].addLe(expr, -left_constant, "Contraint4_13" + i + "," + j);
				}

				// �������
				if (lowerProList[i].solve()) {
					endTime = System.currentTimeMillis();

					if (variableOutput) {

						for (int j = 0; j < input.numBay; j++) {
							for (int j2 = 0; j2 < input.numQC; j2++) {
								YksiHKValue[i][j][j2] = lowerProList[i].getValue(YKsiHK[i][j][j2]);
//								System.out.println("YKsiHK" + i + "," + j + "," + j2 + "Value: " + YksiHKValue[i][j][j2]);
							}
						}

						for (int j = 0; j < jValue; j++) {
							for (int j2 = 0; j2 < input.numQC - 1; j2++) {
								QKsiJKValue[i][j][j2] = lowerProList[i].getValue(QKsiJK[i][j][j2]);
//								System.out.println("QKsiJK" + i + "," + j + "," + j2 + "Value: " + QKsiJKValue[i][j][j2]);
							}
						}

						for (int j = 0; j < input.numQC; j++) {
							RKsiKValue[i][j] = lowerProList[i].getValue(RKsiK[i][j]);
//							System.out.println("RKsiK" + i + "," + j + "Value: " + RKsiKValue[i][j]);
						}

					}

//					if (iteration) {
//						input.writerOverall
//								.write("Iteration " + Iteration + " lowerObj(SD):" + lowerPro.getObjValue() + "\r\n");
//						input.writerOverall.write("\r\n");
	//
//						File file = new File("./result of benchmark/test_" + input.testID);
//						if (!file.exists()) {// ����ļ��в�����
//							file.mkdir();// �����ļ���
//						}
	//
//						if (iterationDetail) {
//							lowerPro.writeSolution("./result of benchmark/test_" + input.testID + "/test_" + input.testID
//									+ " iteration_" + Iteration + "lowerSD.txt");
//						}
//					}

					// ��ȡ��ż����ȡֵ
					// Լ��4_12
					for (int j = 0; j < jValue; j++) {
						for (int j2 = 0; j2 < input.numQC - 1; j2++) {
							paiValue[i][j * (input.numQC - 1) + j2] = lowerProList[i].getDual(lower4_12[i][j][j2]);
//							System.out.println("pai_" + i + '_' + (j * (input.numQC - 1) + j2) + ":"
//									+ paiValue[i][j * (input.numQC - 1) + j2]);
						}
					}
					// Լ��4_13
					for (int j = 0; j < input.numQC; j++) {
						paiValue[i][jValue * (input.numQC - 1) + j] = lowerProList[i].getDual(lower4_13[i][j]);
//						System.out.println("pai_" + i + '_' + (jValue * (input.numQC - 1) + j) + ":"
//								+ paiValue[i][jValue * (input.numQC - 1) + j]);
					}
					// ��ȡĿ�꺯��ֵ
					CmaxKsiLowerValue[i] = lowerProList[i].getValue(CmaxKsiLower[i]);
				} else {
					System.out.println("�²��������ʧ��");
					System.exit(0);
				}
			}
			
			eScalar = 0;
			EVector = new double[num_column];
			for (int i = 0; i < num_ksi; i++) {
				// ����E_{s+1}
				for (int j = 0; j < num_column; j++) {
					for (int j2 = 0; j2 < num_row; j2++) {
						EVector[j] += probabilityKsi[i] * paiValue[i][j2] * matrixT[i][j2][j];
					}
				}

				// ����e_{s+1}
				for (int j2 = 0; j2 < num_row; j2++) {
					eScalar += probabilityKsi[i] * paiValue[i][j2] * vectorH[i][j2];
				}
			}
			
			upperBound = eScalar;
			lowerBound = thetaValue;
			temp_index1 = 0;
			for (int i = 0; i < input.numBay; i++) {
				for (int j = 0; j < input.numQC; j++) {
					upperBound -= EVector[temp_index1] * XhkValue[i][j];
					temp_index1++;
				}
			}
			for (int i = 0; i < input.numQC; i++) {
				upperBound -= EVector[temp_index1] * alphaValue[i];
				temp_index1++;
			}
			for (int i = 0; i < input.numQC; i++) {
				upperBound -= EVector[temp_index1] * betaValue[i];
				temp_index1++;
			}
			
			total_obj = 0;
			for (int i = 0; i < num_ksi; i++) {
				total_obj += probabilityKsi[i] * CmaxKsiLowerValue[i];
			}
			input.writerOverallStochastic.write("Iteration " + Iteration + "upperBound:" + upperBound + "\r\n");
			input.writerOverallStochastic.write("Iteration " + Iteration + "lowerBound:" + lowerBound + "\r\n");
			input.writerOverallStochastic.write("Obj value:" + total_obj + "\r\n");
			input.writerOverallStochastic.write("\r\n");
			System.out.println("�Ͻ�Ϊ" + upperBound);
			System.out.println("�½�Ϊ" + lowerBound);
			System.out.println("Ŀ�꺯��ֵΪ" + total_obj);
		}
		
		// ��¼X��ȡֵ��Ҳ���ǲ�λ�Ͱ��ŵķ����ϵ������alpha��beta��ȡֵ����Ҫ�ȶ����ݽ������
		input.resetWriterStochasticAssignment();
		assignmentRecord(input);
		
		endTime = System.currentTimeMillis();
		long totalComputingTime = endTime - startTime;
		input.writerOverallStochastic.write("Total computing time:" + totalComputingTime + "\r\n");
		input.writerOverallStochastic.write("Iteration:" + Iteration);
		input.writerOverallStochastic.close();
		input.writerStochasticAssignment.close();
	}
	
	// ���ڼ�¼�ϲ����������ָ�ɷ�����Ҳ���Ǳ���X��alpha��beta��ȡֵ
	public void assignmentRecord(Input input) {
		try {
			for (int i = 0; i < XhkValue.length; i++) {
				for (int j = 0; j < XhkValue[0].length; j++) {
					input.writerStochasticAssignment.write("X" + i + "," + j + ":" + (int) Math.round(XhkValue[i][j]) + "\n");
				}
			}

			for (int i = 0; i < alphaValue.length; i++) {
				input.writerStochasticAssignment.write("alpha_" + i + ":" + (int) Math.round(alphaValue[i]) + "\n");
			}

			for (int i = 0; i < betaValue.length; i++) {
				input.writerStochasticAssignment.write("beta_" + i + ":" + (int) Math.round(betaValue[i]) + "\n");
			}

		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

}
