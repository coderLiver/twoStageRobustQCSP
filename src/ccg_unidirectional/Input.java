package ccg_unidirectional;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

public class Input {
	int numBay;
	int numQC;
	double traverseTime;
	double[] processTime; // 每个bay的任务所需时间和
	double uncertain = 1.0; //处理时间的最大波动比率
	int budget;
	double budgetRatio = 0.1; //表示有百分之10的贝位可能取最糟糕的情况
	int safetyMargin;
	FileWriter writerOverall; // 保存算例结果文件的相关写入类，保存的是每步迭代的上下界结果、运行时间
	FileWriter writerAssignment; // 保存最后求得的Robust性最好的方案结果，这里不是保存最后一次迭代的输出结果，而是应该选择迭代过程中上界最小的一次
	FileWriter writerSenario; // 记录writerSchedule对应的迭代次数时下层问题的具体变量Y值和发生attack的贝位
	int testID;
	
	public void dataInput(String pathname, Integer instanceID) {
		int numGroupTask;
		double[] processTimeGroup = null;
		int [] taskLoaction = null;
		testID = instanceID;

		File file = new File("./result of benchmark/overall" + uncertain + "_" + budgetRatio);
		if(!file.exists()){//如果文件夹不存在
			file.mkdirs();//创建文件夹
		}
		
		file = new File("./result of benchmark/assignment" + uncertain + "_" + budgetRatio);
		if(!file.exists()){//如果文件夹不存在
			file.mkdir();//创建文件夹
		}
		
		file = new File("./result of benchmark/senario" + uncertain + "_" + budgetRatio);
		if(!file.exists()){//如果文件夹不存在
			file.mkdir();//创建文件夹
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
				line = line.substring(1, line.length() - 1); //删除首尾字符
				String[] strings = line.split(",");
				processTimeGroup = new double[strings.length];
				for (int i = 0; i < strings.length; i++) {
					processTimeGroup[i] = Double.parseDouble(strings[i]);
				}
			}
			
			if ((line = bufferedReader.readLine()) != null) {
				line = line.substring(1, line.length() - 1); //删除首尾字符
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
			System.out.println("所有任务的总处理时间为：" + totalProcessTime);
			// 后三行数据没用，暂时抛弃
			bufferedReader.close();

			// 进行整贝位转换
			ArrayList<Integer> taskBay = new ArrayList<Integer>(); // 依次不重复地记录每个task对应的bay索引
			ArrayList<Double> processTimeBay = new ArrayList<Double>(); // 与taskBay一一对应的记录每个bay上的任务时间和
			
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
			// 根据设置的比例，取相应的有任务的bay数乘以比例再向上取整
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
	
	// 将Kim and Park benchmark从container group转化成complete bay的形式 
	// 该方法需要在dataInput()方法后调用
	// savePath为保存的路径，以文件夹的名字结尾   saveFileName为要保留的文件名称
	public void instanceReformate(String savePath, String saveFileName) {
		File file = new File(savePath);
		if(!file.exists()){//如果文件夹不存在
			file.mkdir();//创建文件夹
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
	
	// 用于清空文件内容
	public void resetWriterAssignment() {
		String writerAssignmentPath = "./result of benchmark/assignment" + uncertain + "_" + budgetRatio + "/best assignment of test_" + testID + ".txt";
		
		try {
			writerAssignment = new FileWriter(writerAssignmentPath);
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
}
