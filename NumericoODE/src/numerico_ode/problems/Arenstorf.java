package numerico_ode.problems;

import numerico_ode.interpolation.HermiteInterpolator;
import numerico_ode.interpolation.StateFunction;
import numerico_ode.methods.AdaptativeStepMethod;
import numerico_ode.methods.AdaptativeStepRK4Method;
import numerico_ode.methods.AdaptativeStepRKF45Method;
//import numerico_ode.methods.FixedStepEulerMethod;
import numerico_ode.methods.FixedStepMethod;
import numerico_ode.methods.FixedStepModifiedEulerMethod;
import numerico_ode.methods.FixedStepRungeKuttaOrder4Method;
import numerico_ode.ode.InitialValueProblem;
import numerico_ode.ode.NumericalSolutionPoint;
import numerico_ode.tools.BisectionMethod;
import numerico_ode.tools.DisplaySolution;

public class Arenstorf implements InitialValueProblem{
	
	private double mMass = 0.012277471;
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
        //Declaramos los métodos que vamos a usar.
		
		AdaptativeStepMethod method = new AdaptativeStepRKF45Method(problem, 0.001,1e-7, 1, 1e-20);
		AdaptativeStepMethod method2 = new AdaptativeStepRK4Method(problem, 0.001, 1e-6, 1, 1e-20);
		FixedStepMethod method3 = new FixedStepModifiedEulerMethod(problem, 5e-6);
		FixedStepMethod method4 = new FixedStepRungeKuttaOrder4Method(problem, 1e-5);
//		FixedStepMethod method5 = new FixedStepEulerMethod(problem, 1e-6);
		
		
		//Para dar una vuelta, ya que sabemos que el periodo esta alrededor de los 17 segundos
		//Vamos a resolver unos segundos, y a partir de ahí, poner de condición de parada
		//Que se acerque al punto inicial
		NumericalSolutionPoint initialPoint, currentPoint, lastPoint;
		
		System.out.println("METODO DE RUNGE-KUTTA-FEHLBERG");

		double diff = 1;
		initialPoint = currentPoint =  lastPoint = method.getSolution().getLastPoint();
		method.solve(5);
		while (diff > 1e-5 || currentPoint.getState(2) > 0) {
			diff = 0;
			lastPoint = currentPoint;
			currentPoint = method.step();
			diff+= Math.pow(initialPoint.getState(0) - currentPoint.getState(0),2.);
			diff+= Math.pow(initialPoint.getState(2) - currentPoint.getState(2),2.);
		}
		
		// Vamos a interpolar para sacar exactamente el 0 y tener el periodo
		
		StateFunction interpolator = new HermiteInterpolator(problem, lastPoint, currentPoint);
		double periodo = BisectionMethod.findZero(interpolator, lastPoint.getTime(), currentPoint.getTime(), 1e-4, 2);
		System.out.println("Periodo = "+periodo);
		
		// Una vez tenemos el periodo, vamos a dar otra vuelta para ver como cierra.
		method.solve(2*periodo);
		currentPoint = method.getSolution().getLastPoint();
		diff = Math.pow(initialPoint.getState(0) - currentPoint.getState(0),2.) + Math.pow(initialPoint.getState(2) - currentPoint.getState(2),2.);
		diff = Math.sqrt(diff);
		System.out.println("Diferencia al dar 2 vueltas: "+diff);
		System.out.println("Evaluaciones: "+method.getEvaluationCounter());
		
		method.stepPlot();
		
		// Con esto tenemos el periodo. Vamos a probar cada metodo dando 2 vueltas con este periodo para sacar los pasos y las tolerancias.

		/**
		 * 
		 * 
		 * 
		 * 
		 */
		
		System.out.println("METODO DE RUNGE-KUTTA ADAPTATIVO CON EXTRAPOLACIÓN DE RICHARDSON");
		
		diff = 1;
		initialPoint = currentPoint =  lastPoint = method2.getSolution().getLastPoint();
		method2.solve(5);
		while (diff > 1e-5 || currentPoint.getState(2) > 0) {
			diff = 0;
			lastPoint = currentPoint;
			currentPoint = method2.step();
			diff+= Math.pow(initialPoint.getState(0) - currentPoint.getState(0),2.);
			diff+= Math.pow(initialPoint.getState(2) - currentPoint.getState(2),2.);
		}
		
		// Vamos a interpolar para sacar exactamente el 0 y tener el periodo
		
		StateFunction interpolator2 = new HermiteInterpolator(problem, lastPoint, currentPoint);
		periodo = BisectionMethod.findZero(interpolator2, lastPoint.getTime(), currentPoint.getTime(), 1e-4, 2);
		System.out.println("Periodo = "+periodo);
		
		// Una vez tenemos el periodo, vamos a dar otra vuelta para ver como cierra.
		method2.solve(2*periodo);
		currentPoint = method2.getSolution().getLastPoint();
		diff = Math.pow(initialPoint.getState(0) - currentPoint.getState(0),2.) + Math.pow(initialPoint.getState(2) - currentPoint.getState(2),2.);
		diff = Math.sqrt(diff);
		System.out.println("Diferencia al dar 2 vueltas: "+diff);
        
        System.out.println("Evaluaciones: " +method2.getEvaluationCounter());
        
        method2.stepPlot();
        
        /**
         * 
         * 
         * 
         * 
         * 
         */
        
		System.out.println("METODO DE EULER MODIFICADO");

        
		diff = 1;
		initialPoint = currentPoint =  lastPoint = method3.getSolution().getLastPoint();
		method3.solve(5);
		while (diff > 1e-5 || currentPoint.getState(2) > 0) {
			diff = 0;
			lastPoint = currentPoint;
			currentPoint = method3.step();
			diff+= Math.pow(initialPoint.getState(0) - currentPoint.getState(0),2.);
			diff+= Math.pow(initialPoint.getState(2) - currentPoint.getState(2),2.);
		}
		
		// Vamos a interpolar para sacar exactamente el 0 y tener el periodo
		
		StateFunction interpolator3 = new HermiteInterpolator(problem, lastPoint, currentPoint);
		periodo = BisectionMethod.findZero(interpolator3, lastPoint.getTime(), currentPoint.getTime(), 1e-4, 2);
		System.out.println("Periodo = "+periodo);
		
		// Una vez tenemos el periodo, vamos a dar otra vuelta para ver como cierra.
		method3.solve(2*periodo);
		currentPoint = method3.getSolution().getLastPoint();
		diff = Math.pow(initialPoint.getState(0) - currentPoint.getState(0),2.) + Math.pow(initialPoint.getState(2) - currentPoint.getState(2),2.);
		diff = Math.sqrt(diff);
		System.out.println("Diferencia al dar 2 vueltas: "+diff);

		System.out.println("Evaluaciones: "+method3.getEvaluationCounter());

		
		/*
		 * 
		 * 
		 * 
		 */
		
		System.out.println("METODO DE RUNGE-KUTTA DE ORDEN 4");

		
		diff = 1;
		initialPoint = currentPoint =  lastPoint = method4.getSolution().getLastPoint();
		method4.solve(5);
		while (diff > 1e-5 || currentPoint.getState(2) > 0) {
			diff = 0;
			lastPoint = currentPoint;
			currentPoint = method4.step();
			diff+= Math.pow(initialPoint.getState(0) - currentPoint.getState(0),2.);
			diff+= Math.pow(initialPoint.getState(2) - currentPoint.getState(2),2.);
		}
		
		// Vamos a interpolar para sacar exactamente el 0 y tener el periodo
		
		StateFunction interpolator4 = new HermiteInterpolator(problem, lastPoint, currentPoint);
		periodo = BisectionMethod.findZero(interpolator4, lastPoint.getTime(), currentPoint.getTime(), 1e-4, 2);
		System.out.println("Periodo = "+periodo);
		
		// Una vez tenemos el periodo, vamos a dar otra vuelta para ver como cierra.
		method4.solve(2*periodo);
		currentPoint = method4.getSolution().getLastPoint();
		diff = Math.pow(initialPoint.getState(0) - currentPoint.getState(0),2.) + Math.pow(initialPoint.getState(2) - currentPoint.getState(2),2.);
		diff = Math.sqrt(diff);
		System.out.println("Diferencia al dar 2 vueltas: "+diff);
		
		System.out.println("Evaluaciones: "+method4.getEvaluationCounter());
		
		/**
		 * 
		 * 
		 * 
		 * 
		 * 
		 */
		DisplaySolution.statePlot(method.getSolution(), 0, 2);
		
		/**
		 * Método de Euler normal
		 * 
		 */
		
		//method5.solve(2*periodo);
        //DisplaySolution.statePlot(method5.getSolution(), 0, 2,4000);
		//No funciona, el metodo es un desastre para este problema. Con un paso inferior
		//salta error de heap space.
		

	}

}
