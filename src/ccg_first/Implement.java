package ccg_first;

public class Implement {
	public static void main(String[] args) {
		String string = new String("D:\\Better\\Future\\�廪IEѧϰ\\��ҵ���\\testinstance.txt");
		Input data = new Input();
		data.dataInput(string);
		TwoStageRoubustKKT test = new TwoStageRoubustKKT();
		test.upperModel(data.numBay, data.numQC);
		test.lowerModel(data);
		test.ccgIteration(data);
		test.crossForbidCut(data);
	}
}
