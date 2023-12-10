package ccg_unidirectional;

import java.io.IOException;

public class Test {
	public static void main(String[] args) {
		//  Solve deterministic model and attack it or attack given a specific assignment
		StaticMethod.clearFolder("./determinstic");
		for (int i = 90; i <= 90; i++) {
			String path = "KimPark2004instance/test_" + i + ".txt";
			Input data = new Input();
			data.dataInput(path, i);
			Model model = new Model();
			
			String path1 = String.format("result of benchmark/assignment_stochastic/best assignment of test_%d.txt", i);
			model.readAssignmentFromTxt(data, path1);
			
	//		System.out.println("程序运行时间：" + (model.endTime - model.startTime) + "ms"); // 输出程序运行时间
		}
	}
}
