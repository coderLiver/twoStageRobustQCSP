package ccg_unidirectional;

import java.util.ArrayList;

import ilog.concert.IloRange;

public class UncertainSenario {
	int[] ZhSet; // ��ȷ����Zh�����
	int frequency; // ���ֵĴ���
	ArrayList<IloRange> upper1_15Iteration = new ArrayList<IloRange>(); // ��¼�ò�ȷ��������cutʱ��Ӧ��Լ�������ñ���
}
