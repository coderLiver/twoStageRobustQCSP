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
	double[] processTime; // ÿ��bay����������ʱ���
	double uncertain = 1.0; //����ʱ�����󲨶�����
	int budget;
	double budgetRatio = 0.1; //��ʾ�аٷ�֮10�ı�λ����ȡ���������
	int safetyMargin;
	FileWriter writerOverall; // ������������ļ������д���࣬�������ÿ�����������½���������ʱ��
	FileWriter writerAssignment; // ���������õ�Robust����õķ�����������ﲻ�Ǳ������һ�ε�����������������Ӧ��ѡ������������Ͻ���С��һ��
	FileWriter writerSenario; // ��¼writerSchedule��Ӧ�ĵ�������ʱ�²�����ľ������Yֵ�ͷ���attack�ı�λ
	int testID;
	
	public void dataInput(String pathname, Integer instanceID) {
		int numGroupTask;
		double[] processTimeGroup = null;
		int [] taskLoaction = null;
		testID = instanceID;

		File file = new File("./result of benchmark/overall" + uncertain + "_" + budgetRatio);
		if(!file.exists()){//����ļ��в�����
			file.mkdirs();//�����ļ���
		}
		
		file = new File("./result of benchmark/assignment" + uncertain + "_" + budgetRatio);
		if(!file.exists()){//����ļ��в�����
			file.mkdir();//�����ļ���
		}
		
		file = new File("./result of benchmark/senario" + uncertain + "_" + budgetRatio);
		if(!file.exists()){//����ļ��в�����
			file.mkdir();//�����ļ���
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
				line = line.substring(1, line.length() - 1); //ɾ����β�ַ�
				String[] strings = line.split(",");
				processTimeGroup = new double[strings.length];
				for (int i = 0; i < strings.length; i++) {
					processTimeGroup[i] = Double.parseDouble(strings[i]);
				}
			}
			
			if ((line = bufferedReader.readLine()) != null) {
				line = line.substring(1, line.length() - 1); //ɾ����β�ַ�
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
			System.out.println("����������ܴ���ʱ��Ϊ��" + totalProcessTime);
			// ����������û�ã���ʱ����
			bufferedReader.close();

			// ��������λת��
			ArrayList<Integer> taskBay = new ArrayList<Integer>(); // ���β��ظ��ؼ�¼ÿ��task��Ӧ��bay����
			ArrayList<Double> processTimeBay = new ArrayList<Double>(); // ��taskBayһһ��Ӧ�ļ�¼ÿ��bay�ϵ�����ʱ���
			
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
			// �������õı�����ȡ��Ӧ���������bay�����Ա���������ȡ��
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
	
	// ��Kim and Park benchmark��container groupת����complete bay����ʽ 
	// �÷�����Ҫ��dataInput()���������
	// savePathΪ�����·�������ļ��е����ֽ�β   saveFileNameΪҪ�������ļ�����
	public void instanceReformate(String savePath, String saveFileName) {
		File file = new File(savePath);
		if(!file.exists()){//����ļ��в�����
			file.mkdir();//�����ļ���
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
	
	// ��������ļ�����
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
