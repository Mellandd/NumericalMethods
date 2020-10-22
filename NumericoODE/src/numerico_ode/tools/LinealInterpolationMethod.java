package numerico_ode.tools;

import numerico_ode.interpolation.StateFunction;

/**
 * 
 * @author jose
 *
 */


public class LinealInterpolationMethod {
	
	
	/*
	 * We use a linear interpolation based on the Euler method itself. We
	 * calculate the step from the initial time at which the zero is using
	 * the Euler method formula.
	 * 
	 * 0 = y_n + delta * f(t_n, y_n)
	 * delta = -y_n / f(t_n, y_n)
	 * 
	 * And then, the time is t_n + delta
	 */
	
	/**
	 * 
	 * @param function 
	 * @param initialTime
	 * @param index
	 * @param indexDer
	 * @return the time at which zero is found
	 */
	
	
	static public double findZero(StateFunction function, double initialTime, int index, int indexDer) {
		
        double initialState = function.getState(initialTime, index);
        double derivative = function.getState(initialTime, indexDer);
        
        double step = - (initialState / derivative);
        return  initialTime + step;
	}
	
	

}
