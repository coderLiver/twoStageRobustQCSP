package ccg_unidirectional;

import java.io.IOException;

import ilog.concert.IloException;

public class Implement {

	public static void main(String[] args) {
		if (args.length != 1 || args[0].charAt(0) != '-') {
			usage();
			return;
		}

		int startInstanceID = 1;
		int endInstanceID = 10;
		String savePath = "./KimPark2004instance/reformate"; // 按bay改造后的benchmark的保存地址
		double[] uncertain = new double[]{0.1, 0.2, 0.3, 0.4, 0.5, 1.0};
		double[] budgetRatio = new double[] {0.1, 0.2, 0.3, 0.4};

		switch (args[0].charAt(1)) {
		case 'a':
			for (int i = startInstanceID; i <= endInstanceID; i++) {
				String path = "KimPark2004instance/test_" + i + ".txt";
				Input data = new Input();
				data.dataInput(path, i);
				data.resetWriterOverall();
				data.resetWriterAssignment();
				data.resetWriterSenario();
				Model model = new Model();
				
				try {
					// 基于强对偶进行迭代
					model.upperModel(data, false, false);
					model.lowerModelSD(data, true, false, false);
					model.ccgIteration(data, "SD",false, false);
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				
				System.out.println("程序运行时间：" + (model.endTime - model.startTime) + "ms"); // 输出程序运行时间
			}
			break;
		case 'b':
			for (int i = startInstanceID; i <= endInstanceID; i++) {
				String path = "KimPark2004instance/test_" + i + ".txt";
				Input data = new Input();
				data.dataInput(path, i);
				data.resetWriterOverall();
				data.resetWriterAssignment();
				data.resetWriterSenario();
				Model model = new Model();
				
				try {
					// 基于KKT进行迭代
					model.upperModel(data, false, false);
					model.lowerModelKKT(data, false, false);
					model.ccgIteration(data, "KKT", false, false);
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				
				System.out.println("程序运行时间：" + (model.endTime - model.startTime) + "ms"); // 输出程序运行时间
			}
			break;
		case 'c':
			
			for (int i = startInstanceID; i <= endInstanceID; i++) {
				String path = "KimPark2004instance/test_" + i + ".txt";
				Input data = new Input();
				data.dataInput(path, i);
				data.resetWriterOverall();
				data.resetWriterAssignment();
				data.resetWriterSenario();
				Model model = new Model();
				
				try {
					// Benders Decomposition
					model.upperModel(data, false, false);
					model.lowerModelSD(data, true, false, false);
					model.BDIteration(data, true, false);
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				
				System.out.println("程序运行时间：" + (model.endTime - model.startTime) + "ms"); // 输出程序运行时间
			}
			break;
		case 'd':
			StaticMethod.clearFolder("./result of benchmark/uncertaintyRecord");
			
			for (int i = startInstanceID; i <= endInstanceID; i++) {
				String path = "KimPark2004instance/test_" + i + ".txt";
				Input data = new Input();
				data.dataInput(path, i);
				data.resetWriterOverall();
				data.resetWriterAssignment();
				data.resetWriterSenario();
				Model model = new Model();
				
				try {
					// Benders与C&CG集合
					model.upperModel(data, false, false);
					model.lowerModelSD(data, true, false, false);
					model.CCGBDIteration(data, false, false);
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				
				System.out.println("程序运行时间：" + (model.endTime - model.startTime) + "ms"); // 输出程序运行时间
			}
			break;
		case 'e':
			
			for (int i = startInstanceID; i <= endInstanceID; i++) {
				String path = "KimPark2004instance/test_" + i + ".txt";
				Input data = new Input();
				data.dataInput(path, i);
				data.resetWriterOverall();
				data.resetWriterAssignment();
				data.resetWriterSenario();
				Model model = new Model();
				
				try {
					// 以确定性为初始解的BDC&CG
					model.deterministicInitial(data, false, false);
					model.lowerModelSD(data, true, false, false);
					model.CCGBDIteration(data, false, false);
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				
				System.out.println("程序运行时间：" + (model.endTime - model.startTime) + "ms"); // 输出程序运行时间
			}
			break;
		case 'f':
			for (int i = startInstanceID; i <= endInstanceID; i++) {
				String path = "KimPark2004instance/test_" + i + ".txt";
				Input data = new Input();
				data.dataInput(path, i);
				Model model = new Model();
				
				try {
					// 确定性模型迭代
					model.deterministicModel(data, false, false);
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				
				System.out.println("程序运行时间：" + (model.endTime - model.startTime) + "ms"); // 输出程序运行时间
			}
			break;
		case 'g':
			
			StaticMethod.clearFolder("./determinstic");
			for (int i = startInstanceID; i <= endInstanceID; i++) {
				String path = "KimPark2004instance/test_" + i + ".txt";
				Input data = new Input();
				data.dataInput(path, i);
				Model model = new Model();
				
				try {
					// 确定性模型迭代 + Attack
					model.deterministicModel(data, false, false);
					model.deterministicAttackModel(data, false);
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				
				System.out.println("程序运行时间：" + (model.endTime - model.startTime) + "ms"); // 输出程序运行时间
			}
			break;
		case 'h':
			
			for (int i = startInstanceID; i <= endInstanceID; i++) {
				String path = "KimPark2004instance/test_" + i + ".txt";
				String saveFileName = "reformateTest_" + i + ".txt";
				Input data = new Input();
				data.dataInput(path, i);
				Model model = new Model();
				
				// 按bay改造benchmark
				data.instanceReformate(savePath, saveFileName);
				
				System.out.println("程序运行时间：" + (model.endTime - model.startTime) + "ms"); // 输出程序运行时间
			}
			break;
		case 'i':
			for (int i = startInstanceID; i <= endInstanceID; i++) {
				String path = "KimPark2004instance/test_" + i + ".txt";
				Input data = new Input();
				data.dataInput(path, i);
				Model model = new Model();
				
				// 固定X和Z变量后进行求解
				try {
					model.assignmentUncertaintyFix(data);
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			break;
		case 'j':
			// 基于BD进行uncertainty和budgetRatio的灵敏度分析
			for (int k = 0; k < uncertain.length; k++) {
				for (int j = 0; j < budgetRatio.length; j++) {
					for (int i = startInstanceID; i <= endInstanceID; i++) {
						String path = "KimPark2004instance/test_" + i + ".txt";
						Input data = new Input();
						data.uncertain = uncertain[k];
						data.budgetRatio = budgetRatio[j];
						data.dataInput(path, i);
						data.resetWriterOverall();
						data.resetWriterAssignment();
						data.resetWriterSenario();
						Model model = new Model();
						
						try {
							// Benders Decomposition
							model.upperModel(data, false, false);
							model.lowerModelSD(data, true, false, false);
							model.BDIteration(data, true, false);
						} catch (IOException e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
						}
						
						System.out.println("程序运行时间：" + (model.endTime - model.startTime) + "ms"); // 输出程序运行时间
					}
				}
			}
			break;
		case 'k':
//			for (int i = startInstanceID; i <= endInstanceID; i++) {
			for (int i = 1; i <= 1; i++) {
				String path = "KimPark2004instance/test_" + i + ".txt";
				Input data = new Input();
				data.dataInput(path, i);
				data.resetWriterOverall();
				data.resetWriterAssignment();
				data.resetWriterSenario();
				StochasticProgrammingModel model = new StochasticProgrammingModel();
				// 基于Stochastic Programming进行建模
				try {
					model.upperModel(data, false, true);
					model.lowerModelLShape(data, false, true);
					model.lowerModelIterate(data, false, true);
				} catch (IOException | IloException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				
				System.out.println("程序运行时间：" + (model.endTime - model.startTime) + "ms"); // 输出程序运行时间
			}
			break;
		default:
			usage();
			return;
		}
	}

	static void usage() {
		System.out.println("usage:   Implement <option>");
		System.out.println("options:       -a   Solve using C&CG based on Strong Duality");
		System.out.println("options:       -b   Solve using C&CG baesd on KKT condition");
		System.out.println("options:       -c   Solve using Benders Decomposition");
		System.out.println("options:       -d   Solve using the Combination of Benders Decomposition and C&CG (BDC&CG)");
		System.out.println("options:       -e   Solve using the result of deterministic model as initial solution of BDC&CG");
		System.out.println("options:       -f   Solve deterministic model");
		System.out.println("options:       -g   Solve deterministic model and attack it");
		System.out.println("options:       -h   reformate the benchmark from container group to complete bay");
		System.out.println("options:       -i   Given Assignment and uncertainty and solve");
	}
}
