package numerico_ode.problems;

import numerico_ode.interpolation.StateFunction;
import numerico_ode.methods.*;
import numerico_ode.ode.*;
import numerico_ode.tools.*;
import java.lang.Math;

public class Sho implements InitialValueProblem{

	private final double mass = 1;
	private final double k = 1.5;
	private final double l = 0.7;
	private final double a = 0.4;
	private final double w;
	private final double b;
	

	
	public Sho(double w, double b) {
		this.w = w;
		this.b = b;
	}
	
	@Override
	public double getInitialTime() {
		return 0;
		
	}

	@Override
	public double[] getInitialState() {
		return new double[] {1.5,0,0,0}; 
	}

	@Override
	public double[] getDerivative(double t, double[] x) {
		double f = a* Math.sin(w*t);
		return new double[] {
				x[1], (-k*(x[0]-l)- b*x[1] + f)/mass, 0 ,0 };
		}
	
	
	/*
	 * La resolvemos análiticamente para el caso sin rozamiento ni forzado.
	 * 
	 * x(t) = A cos(wt + phi) donde A es la amplitud máxima, w = sqrt(k/m),
	 * phi es la fase inicial.
	 * 
	 * v(t) = -wA sen(wt + phi).
	 * 
	 */
    static private class TrueSol implements StateFunction {
    	double amplitud = 0.8;
    	double faseInicial = 0;
    	double wAn = Math.sqrt(1.5);
        public double[] getState(double time) {
            return new double[] { 0.7 + amplitud*Math.cos(wAn*time + faseInicial), 
            		-wAn * amplitud * Math.sin(wAn*time + faseInicial), 0, 0};
        }
        public double getState(double time, int index) {
            switch (index) {
                case 0 : return 0.7 + amplitud*Math.cos(wAn*time + faseInicial);
                case 1 : return -wAn * amplitud * Math.sin(wAn*time + faseInicial);
                case 2 : return 0;
                case 3 : return 0;
                default : return Double.NaN;
            }
        }
}
	
	public static void main(String[] args) {
		
		/**
		 * Euler para ecuacion 3, y dibujo de la gráfica
		 * de x y vx.
		 */
		InitialValueProblem problem = new Sho(1.3, 0.3);
		FixedStepMethod method = new FixedStepEulerMethod(problem,1.0e-2);
		method.solve(40);
        DisplaySolution.timePlot(method.getSolution(), new int[]{1});
        DisplaySolution.timePlot(method.getSolution(), new int[]{0});
        
        /*
         * Valido el método para el caso sin rozamiento ni forzado.
         */
        
        InitialValueProblem problem2 = new Sho(0, 0);
        FixedStepMethod method2 = new FixedStepEulerMethod(problem2,1.0e-4);
        method2.solve(40);
        
        //DisplaySolution.listError(method2.getSolution(), new TrueSol(), new int[]{0, 1});
        
        /*
         * Con un paso de 1.0e-4 e genera un error de magnitud 1e-5.
         */
        
        
        /**
         * Frecuencia del movimiento no forzado.
         * 
         * Sacamos el primer momento a parte del inicial donde la velocidad es 0
         * y x > l.
         */
        
        InitialValueProblem problem3 = new Sho(0, 0.3);
        FixedStepMethod method3 = new FixedStepEulerMethod(problem3,1.0e-4);
        NumericalSolutionPoint previousPoint, currentPoint;
        double tolerance = 1e-5;
        previousPoint = currentPoint = method3.getSolution().getLastPoint();
        currentPoint = method3.step();
        while (currentPoint.getState(0) <= 0.7 || Math.abs(currentPoint.getState(1)) > tolerance) {
            previousPoint = currentPoint;
            currentPoint = method3.step();
        }
        currentPoint.println();
        
        double period = currentPoint.getTime();
        double frecuence = 1 / period;
        
        /*
         * Usamos la frecuencia nueva como w para forzarlo.
         *
         */
        
        InitialValueProblem problem4 = new Sho(frecuence, 0.3);

        FixedStepMethod method4 = new FixedStepEulerMethod(problem4,1.0e-2);
        method4.solve(40);
        DisplaySolution.timePlot(method4.getSolution(), new int[]{1});
        DisplaySolution.timePlot(method4.getSolution(), new int[]{0});
        
	}
		
		
}
