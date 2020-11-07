package numerico_ode.problems;

import numerico_ode.interpolation.*;
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
        
        //Obtenemos las órbitas para euler y euler modificado.
        method2.solve(9000); 
        method3.solve(9000);
        DisplaySolution.statePlot(method2.getSolution(), 0, 2);
        DisplaySolution.statePlot(method3.getSolution(), 0, 2);
        
        //Empezamos a calcular el perifelio con runge-kutta.
        previousPoint = currentPoint = method.getSolution().getLastPoint();
        while (currentPoint.getState(2)>=0) {
            previousPoint = currentPoint;
            currentPoint = method.step();
        }
        previousPoint.println();
        currentPoint.println();
        
        //Ahora interpolamos para hallar exactamente el 0.
        StateFunction interpolator = new EulerMethodInterpolator(problem, previousPoint);
        double zeroYAt = BisectionMethod.findZero(interpolator, previousPoint.getTime(), currentPoint.getTime(), 1.0e-8, 2);
        if (Double.isNaN(zeroYAt)) {
        	System.out.println("Perihelion not found!");
        } else {
        	double periphelio = interpolator.getState(zeroYAt, 0);
        	System.out.println("Perihelion at x= "+periphelio);
        }
        
        //DisplaySolution.statePlot(method.getSolution(), 0, 2); //Gráfica que muestra de afelio a perifelio
        
        //Seguimos ahora para obtener la vuelta completa y calcular el año.
        while(currentPoint.getState(2)<0) {
            previousPoint = currentPoint;
            currentPoint = method.step();
        }
        previousPoint.println();
        currentPoint.println();
        
        //Interpolamos para obtenerlo exactamente.
        StateFunction interpolator2 = new EulerMethodInterpolator(problem, previousPoint);
        double zeroYAt2 = BisectionMethod.findZero(interpolator2, previousPoint.getTime(), currentPoint.getTime(), 1.0e-8, 2);
        if (Double.isNaN(zeroYAt2)) {
        	System.out.println("Orbit not found");
        } else {
        	System.out.println("Year at: "+ zeroYAt2/24. +" days");
        }
        
        DisplaySolution.statePlot(method.getSolution(), 0, 2); //Obtenemos la órbita completa con Runge-Kutta
        
	}

}
