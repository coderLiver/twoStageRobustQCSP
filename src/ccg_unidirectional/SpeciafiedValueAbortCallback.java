package ccg_unidirectional;

import ilog.concert.IloException;
import ilog.cplex.IloCplex;

public class SpeciafiedValueAbortCallback extends IloCplex.MIPInfoCallback{
	IloCplex _cplex;
	double _abortValue;
	private static final double ACCURACY = 0.0001;
	
	public SpeciafiedValueAbortCallback(IloCplex cplex, double abortValue) {
		// TODO Auto-generated constructor stub
		_cplex = cplex;
		_abortValue = abortValue;
	}
	protected void main() throws IloException {
		// TODO Auto-generated method stub
		if (getIncumbentObjValue() <= _abortValue + ACCURACY) {
			abort();
		}
	}

}
