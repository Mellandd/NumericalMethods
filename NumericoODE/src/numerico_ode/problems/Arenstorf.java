package numerico_ode.problems;

import numerico_ode.methods.AdaptativeStepMethod;
import numerico_ode.methods.AdaptativeStepRKF45Method;
import numerico_ode.ode.InitialValueProblem;
import numerico_ode.tools.DisplaySolution;

public class Arenstorf implements InitialValueProblem{
	
	private double mMass = 0.01227741;
	private double eMass = 1-mMass;

	@Override
	public double getInitialTime() {
		return 0;
	}

	@Override
	public double[] getInitialState() {
		return new double[] { 0.994, 0, 0, -2.001585106 };
	}

	@Override
	public double[] getDerivative(double t, double[] x) {
		double d1 = Math.pow(Math.pow(x[0]+mMass,2)+x[2]*x[2],3./2.);
		double d2 = Math.pow(Math.pow(x[0]-eMass,2)+x[2]*x[2],3./2.);
        return new double[] { x[1], 
        		x[0]+2*x[3]-eMass*(x[0]+mMass)/d1-mMass*(x[0]-eMass)/d2, 
                x[3], 
        		x[2]-2*x[1]-eMass*(x[2])/d1-mMass*(x[2])/d2
        };
	}
	
	public static void main(String[] args) {
		InitialValueProblem problem = new Arenstorf();
		AdaptativeStepMethod method = new AdaptativeStepRKF45Method(problem, 0.00005,1e-8,1,1e-50);
		method.solve(18);
        DisplaySolution.statePlot(method.getSolution(), 0, 2,100);
	}

}
