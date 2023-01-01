package ccg_unidirectional;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

// 本类用于录入岸桥与贝位指派关系的X变量和代表uncertainty取值的Z变量的取值
public class Schedule {
	double[] ZhValue;
	double[][] XhkValue;
	double[] alphaValue;
	double[] betaValue;
	int numBay;
	int numQC;
	
	public Schedule(int numBay, int numQC){
		this.numBay = numBay;
		this.numQC = numQC;
		ZhValue = new double[numBay];
		XhkValue = new double[numBay][numQC];
		alphaValue = new double[numQC];
		betaValue = new double[numQC];
	}
	
	// assignmentPath是保留X变量相关取值文件的路径，uncertaintyPath是保留Z变量相关取值文件的路径
	public void dataInput(String assignmentPath, String uncertaintyPath) throws IOException {
		File assignmentFile = new File(assignmentPath);
		FileReader fileReader = new FileReader(assignmentFile);
		BufferedReader bufferedReader = new BufferedReader(fileReader);

		
		for (int i = 0; i < numBay; i++) {
			for (int j = 0; j < numQC; j++) {
				String line = bufferedReader.readLine();
				line = line.split(":")[1];
				XhkValue[i][j] = Double.parseDouble(line);
			}
		}
		
		for (int i = 0; i < numQC; i++) {
			String line = bufferedReader.readLine();
			line = line.split(":")[1];
			alphaValue[i] = Double.parseDouble(line);
		}
		
		for (int i = 0; i < numQC; i++) {
			String line = bufferedReader.readLine();
			line = line.split(":")[1];
			betaValue[i] = Double.parseDouble(line);
		}
		
		bufferedReader.close();
		
		File uncertaintyFile = new File(uncertaintyPath);
		FileReader fileReader1 = new FileReader(uncertaintyFile);
		BufferedReader bufferedReader1 = new BufferedReader(fileReader1);
		
		for (int i = 0; i < numBay; i++) {
			String line = bufferedReader1.readLine();
			line = line.split(":")[1];
			ZhValue[i] = Double.parseDouble(line);
		}
		
		bufferedReader1.close();
	}
}
