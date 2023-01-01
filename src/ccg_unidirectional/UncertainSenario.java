package ccg_unidirectional;

import java.util.ArrayList;

import ilog.concert.IloRange;

public class UncertainSenario {
	int[] ZhSet; // 不确定集Zh的组合
	int frequency; // 出现的次数
	ArrayList<IloRange> upper1_15Iteration = new ArrayList<IloRange>(); // 记录该不确定集加入cut时对应的约束的引用变量
}
