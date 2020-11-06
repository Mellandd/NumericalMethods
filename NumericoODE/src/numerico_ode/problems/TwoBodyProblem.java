package numerico_ode.problems;

import numerico_ode.methods.*;
import numerico_ode.ode.*;
import numerico_ode.tools.*;

public class TwoBodyProblem implements InitialValueProblem{
	private final double constantG = 8.649828e-4;
	private final double mTierra = 0.0059724;
	private final double mSol = 1988.500;
	private final double dAfelio = 152.1;
	private final double vAfelio = 0.105444;
	

	@Override
	public double getInitialTime() {
		return 0;
	}

	@Override
	public double[] getInitialState() { //x,vx,y,vy
		return new double[] { dAfelio, 0, 0, vAfelio };
	}

	@Override
	public double[] getDerivative(double t, double[] x) {
        return new double[] { x[1], 
        		-constantG*(mTierra+mSol)*(x[0]/Math.pow(x[0]*x[0]+x[2]*x[2],3./2.)), 
                x[3], 
                -constantG*(mTierra+mSol)*(x[2]/Math.pow(x[0]*x[0]+x[2]*x[2],3./2.))};
	}
	
	public static void main(String[] args) {
		InitialValueProblem problem = new TwoBodyProblem();
		FixedStepMethod method = new FixedStepRungeKuttaOrder4Method(problem, 1);
		FixedStepMethod method2 = new FixedStepEulerMethod(problem, 1);
		FixedStepMethod method3 = new FixedStepModifiedEulerMethod(problem, 1);
		
        NumericalSolutionPoint previousPoint, currentPoint;
        //method.solve(9000);
        previousPoint = currentPoint = method.getSolution().getLastPoint();
        double x0 = currentPoint.getState(0);
        currentPoint = method.step();
        double tolerance = 5e-2;
        //double y = previousPoint.getState(2);
        while (Math.pow(currentPoint.getState(0)-x0,2) + Math.pow(currentPoint.getState(2),2)>tolerance) {
            previousPoint = currentPoint;
            currentPoint = method.step();
        }
        previousPoint.println();
        currentPoint.println();
        DisplaySolution.statePlot(method.getSolution(), 0, 2);
	}

}
