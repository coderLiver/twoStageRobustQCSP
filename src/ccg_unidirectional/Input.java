package ccg_unidirectional;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;

public class Input {
	int numBay;
	int numQC;
	double traverseTime;
	double[] processTime; // 贝位所需加工时间
	double uncertain = 1.0; //任务的最大偏离程度
	int budget;
	double budgetRatio = 0.1; // 偏移任务的最大比例
	int safetyMargin;
	FileWriter writerOverall; // 记录迭代过程中的上下界，即主问题子问题相关的目标函数
	FileWriter writerAssignment; // 记录鲁棒模型中贝位岸桥的指派关系
	FileWriter writerStochasticAssignment; // 记录随机模型中贝位岸桥的指派关系
	FileWriter writerSenario; // 记录发生加工时间偏移的贝位
	FileWriter writerOverallStochastic; // 记录随机优化模型迭代过程中的上下界
	int testID;
	// 闅忔満浼樺寲閮ㄥ垎
	int scenarioEachBay = 2;
	int num_ksi; // the amount of scenarios
	double[][] stochasticProcessTime; // i-scenario_index, j-processing time of bay 
	double multiple = 1.5; // multiple of stochastic processing time 
	
	public void dataInput(String pathname, Integer instanceID) {
		int numGroupTask;
		double[] processTimeGroup = null;
		int [] taskLoaction = null;
		testID = instanceID;

		File file = new File("./result of benchmark/overall" + uncertain + "_" + budgetRatio);
		if(!file.exists()){
			file.mkdirs();
		}
		
		file = new File("./result of benchmark/overall_stochastic");
		if(!file.exists()){
			file.mkdirs();
		}
		
		file = new File("./result of benchmark/assignment" + uncertain + "_" + budgetRatio);
		if(!file.exists()){
			file.mkdir();
		}
		
		file = new File("./result of benchmark/assignment_stochastic");
		if(!file.exists()){
			file.mkdir();
		}
		
		file = new File("./result of benchmark/senario" + uncertain + "_" + budgetRatio);
		if(!file.exists()){
			file.mkdir();
		}
		
		try {
			file = new File(pathname);
			FileReader fileReader = new FileReader(file);
			BufferedReader bufferedReader = new BufferedReader(fileReader);
			String line = null;
			
			if((line = bufferedReader.readLine()) != null) {
				numGroupTask = Integer.parseInt(line);
			}
			
			if ((line = bufferedReader.readLine()) != null) {
				numBay = Integer.parseInt(line);
			}
			
			if ((line = bufferedReader.readLine()) != null) {
				numQC = Integer.parseInt(line);
			}
			
			if ((line = bufferedReader.readLine()) != null) {
				traverseTime = Integer.parseInt(line);
			}
			
			if((line = bufferedReader.readLine()) != null) {
				safetyMargin = Integer.parseInt(line);
			}
			
			if ((line = bufferedReader.readLine()) != null) {
				line = line.substring(1, line.length() - 1); //去除首尾字符
				String[] strings = line.split(",");
				processTimeGroup = new double[strings.length];
				for (int i = 0; i < strings.length; i++) {
					processTimeGroup[i] = Double.parseDouble(strings[i]);
				}
			}
			
			if ((line = bufferedReader.readLine()) != null) {
				line = line.substring(1, line.length() - 1); //去除首尾字符
				String[] strings = line.split(",");
				taskLoaction = new int[strings.length];
				for (int i = 0; i < strings.length; i++) {
					taskLoaction[i] = Integer.parseInt(strings[i]);
				}
			}
			
			double totalProcessTime = 0;
			for (int i = 0; i < processTimeGroup.length; i++) {
				totalProcessTime += processTimeGroup[i];
			}
			System.out.println("Total processing time of all the bays:" + totalProcessTime);
			bufferedReader.close();

			// 将group为粒度的task转化为bay为粒度的task
			ArrayList<Integer> taskBay = new ArrayList<Integer>(); // 记录task对应的贝位位置
			ArrayList<Double> processTimeBay = new ArrayList<Double>(); // 记录有task贝位的任务时间和
			
			for (int i = 0; i < taskLoaction.length; i++) {
				int j = taskBay.size();
				Boolean temp = false;
				for (int j2 = 0; j2 < j; j2++) {
					if (taskLoaction[i] == taskBay.get(j2)) {
						processTimeBay.set(j2, processTimeBay.get(j2) + processTimeGroup[i]);
						temp = true;
						break;
					}
				}
				if (!temp) {
					taskBay.add(taskLoaction[i]);
					processTimeBay.add(processTimeGroup[i]);
				}
			}
			
//			budget = 0;
			// 可以发生uncertainty的贝位数量取有任务的贝位数*budgetRatio，再向上取整。
			budget = (int) Math.ceil(taskBay.size() * budgetRatio);

			processTime = new double[numBay];
			for (int i = 0; i < taskBay.size(); i++) {
				processTime[taskBay.get(i) - 1] = processTimeBay.get(i);
			}
			
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	
	// 灏咾im and Park benchmark浠巆ontainer group杞寲鎴恈omplete bay鐨勫舰寮� 
	// 璇ユ柟娉曢渶瑕佸湪dataInput()鏂规硶鍚庤皟鐢�
	// savePath涓轰繚瀛樼殑璺緞锛屼互鏂囦欢澶圭殑鍚嶅瓧缁撳熬   saveFileName涓鸿淇濈暀鐨勬枃浠跺悕绉�
	public void instanceReformate(String savePath, String saveFileName) {
		File file = new File(savePath);
		if(!file.exists()){//濡傛灉鏂囦欢澶逛笉瀛樺湪
			file.mkdir();//鍒涘缓鏂囦欢澶�
		}
		
		try {
			FileWriter writer = new FileWriter(savePath + "/" + saveFileName);
			writer.write("[");
			
			for (int i = 0; i < processTime.length - 1; i++) {
				writer.write(processTime[i] + ",");
			}
			
			writer.write(processTime[processTime.length - 1] + "]");
			writer.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	// 鐢ㄤ簬娓呯┖鏂囦欢鍐呭
	public void resetWriterAssignment() {
		String writerAssignmentPath = "./result of benchmark/assignment" + uncertain + "_" + budgetRatio + "/best assignment of test_" + testID + ".txt";
		
		try {
			writerAssignment = new FileWriter(writerAssignmentPath);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public void resetWriterStochasticAssignment() {
		String writerStochasticAssignmentPath = "./result of benchmark/assignment_stochastic" + "/best assignment of test_" + testID + ".txt";
		
		try {
			writerStochasticAssignment = new FileWriter(writerStochasticAssignmentPath);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public void resetWriterSenario() {
		String writerSenarioPath = "./result of benchmark/senario" + uncertain + "_" + budgetRatio + "/uncertainty senario of test_" + testID + ".txt";
		
		try {
			writerSenario = new FileWriter(writerSenarioPath);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public void resetWriterOverall() {
		String writerOverallPath = "./result of benchmark/overall" + uncertain + "_" + budgetRatio + "/test_" + testID + ".txt";
		
		try {
			writerOverall = new FileWriter(writerOverallPath);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public void resetWriterOverallStochastic() {
		String writerOverallPath = "./result of benchmark/overall_stochastic" + "/test_" + testID + ".txt";
		
		try {
			writerOverallStochastic = new FileWriter(writerOverallPath);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	// all the bays can be stochastic
	public void stochasticProcessingGeneration() {
		num_ksi = (int) Math.pow(scenarioEachBay, numBay);
		stochasticProcessTime = new double[num_ksi][numBay];
		String binary_i;
		for (int i = 0; i < num_ksi; i++) {
			// Convert i to binary. Each bit represents that the bay is stochastic or not.
			binary_i = Integer.toString(i, 2);
			String formatString = "%" + numBay + "s";
			binary_i = String.format(formatString, binary_i).replace(" ", "0");
			for (int j = 0; j < binary_i.length(); j++) {
				if (binary_i.charAt(j) == '0') {
					stochasticProcessTime[i][j] = processTime[j];
				} else if (binary_i.charAt(j) == '1') {
					stochasticProcessTime[i][j] = processTime[j] * multiple;
				}else {
					System.out.println("Unexpected Value: not 0 and 1!");
				}
			}
		}
	}
	
	/**
	 *  从加工时间不为0的贝位中随机抽取最多不超过selectedBay个，并保存结果到txt中
	 * @param selectedBay
	 */
	public void generateStochasticBaySave(Integer selectedBay, String readTxtPath, String writeTxtPath) {
		ArrayList<Integer> processBayArrayList = new ArrayList<Integer>();
		File file = new File(readTxtPath);
		FileReader fileReader;
		try {
			fileReader = new FileReader(file);
			BufferedReader bufferedReader = new BufferedReader(fileReader);
			String content = bufferedReader.readLine();
			// 查找加工时间不为0的贝位
			content = content.substring(1, content.length() - 1);
			String[] tempProcessTime = content.split(",");
			
			for (int j = 0; j < tempProcessTime.length; j++) {
				double temp_double = Double.parseDouble(tempProcessTime[j]);
				if (Math.abs(temp_double - 0) > 0.001) {
					processBayArrayList.add(j);
				}
			}
			
			bufferedReader.close();
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		// 随机抽取最多k个
		Random random = new Random();
		ArrayList<Integer> resultArrayList = new ArrayList<Integer>();
		int list_size = processBayArrayList.size();
		if (list_size <= selectedBay) {
			resultArrayList = processBayArrayList;
		}else {
			for (int i = 0; i < selectedBay; i++) {
				int index = random.nextInt(list_size - i);
				resultArrayList.add(processBayArrayList.get(index));
				Collections.swap(processBayArrayList,index, list_size - i - 1);
			}
		}
		
		File writeFile = new File(writeTxtPath);
		
		if(!writeFile.exists()){//濡傛灉鏂囦欢澶逛笉瀛樺湪
			writeFile.getParentFile().mkdir();//鍒涘缓鏂囦欢澶�
		}
		
		try {
			FileWriter writer = new FileWriter(writeFile);
			
			for (int i = 0; i < resultArrayList.size(); i++) {
				writer.write(resultArrayList.get(i).toString());
				if (i != resultArrayList.size() - 1) {
					writer.write(',');
				}
			}
			
			writer.close();
		} catch (IOException e) {
			e.printStackTrace();
			// TODO: handle exception
		}
	}
	
	// 获取txt文件中保存的stochastic bay的位置，生成好stochasticProcessTime，用于建模
	public void getScenarioStochasticModel(String readTxtPath) {
		File file = new File(readTxtPath);
		FileReader fileReader;
		try {
			fileReader = new FileReader(file);
			BufferedReader bufferedReader = new BufferedReader(fileReader);
			String content = bufferedReader.readLine();
			String[] stochasticBaysStrings = content.split(",");
			int[] stochasticBaysLocation = new int[stochasticBaysStrings.length];
			bufferedReader.close();
			num_ksi = (int) Math.pow(2, stochasticBaysLocation.length);
			
			for (int i = 0; i < stochasticBaysStrings.length; i++) {
				stochasticBaysLocation[i] = Integer.parseInt(stochasticBaysStrings[i]);
			}
			
			// innitial stochasticProcessTime
			stochasticProcessTime = new double[num_ksi][processTime.length];
			
			for (int i = 0; i < num_ksi; i++) {
				stochasticProcessTime[i] = processTime.clone();
			}
			
			String binary_i;
			for (int i = 0; i < num_ksi; i++) {
				binary_i = Integer.toString(i, 2);
				String formatString = "%" + stochasticBaysLocation.length + "s";
				binary_i = String.format(formatString, binary_i).replace(" ", "0");
				for (int j = 0; j < binary_i.length(); j++) {
					if (binary_i.charAt(j) == '1') {
						stochasticProcessTime[i][stochasticBaysLocation[j]] = processTime[stochasticBaysLocation[j]] * multiple;
					}
				}
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
}
