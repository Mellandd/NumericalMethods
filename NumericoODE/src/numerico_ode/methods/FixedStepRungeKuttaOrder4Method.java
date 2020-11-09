package numerico_ode.methods;

import java.util.stream.IntStream;

import numerico_ode.ode.InitialValueProblem;

public class FixedStepRungeKuttaOrder4Method  extends FixedStepMethod {
	
    /**
     * Initializes the method for a given InitialValueProblem
     * @param InitialValueProblem problem 
     * @param step the fixed step to take. If negative, we'd solve backwards in time
     */
    public FixedStepRungeKuttaOrder4Method(InitialValueProblem problem, double step) {
        super(problem,step);
    }

	@Override
    /**
     * Runge-Kutta method implementation of order 4.
     * @param deltaTime the step to take
     * @param time the current time
     * @param state the current state
     * @return the value of time of the step taken, state will contain the updated state
     */
	protected double doStep(double deltaTime, double time, double[] state) {
		super.addToEvaluationCounter(4);
		double[] k1 = mProblem.getDerivative(time, state);
        double[] stateAux = new double[state.length];
        // We use Streams, from Java 8 and above, to do the vectorial operation
        // since we need state + 1/2 k1 to do k2.
        IntStream.range(0, state.length).forEach(i-> stateAux[i] = state[i]+ 0.5*k1[i]);
		double[] k2 = mProblem.getDerivative(time + deltaTime/2., stateAux);
        IntStream.range(0, state.length).forEach(i-> stateAux[i] = state[i]+ 0.5*k2[i]);
		double[] k3 = mProblem.getDerivative(time + deltaTime/2., stateAux);
        IntStream.range(0, state.length).forEach(i-> stateAux[i] = state[i]+ k3[i]);
		double[] k4 = mProblem.getDerivative(time + deltaTime, stateAux);
		// With everything calculated, we calculate our need step with the formula given
		// by the Runge-Kutta method.
		IntStream.range(0, state.length).forEach(i->state[i] = state[i] + (1./6.)*deltaTime*(k1[i]+2*k2[i]+2*k3[i]+k4[i]));
		return time + deltaTime;
	}

}
