package numerico_ode.methods;

import numerico_ode.ode.InitialValueProblem;

public class FixedStepRungeKuttaOrder4Method  extends FixedStepMethod {
	protected double[] auxState;
	
    /**
     * Initializes the method for a given InitialValueProblem
     * @param InitialValueProblem problem 
     * @param step the fixed step to take. If negative, we'd solve backwards in time
     */
    public FixedStepRungeKuttaOrder4Method(InitialValueProblem problem, double step) {
        super(problem,step);
        auxState = problem.getInitialState();
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
        double h2 = deltaTime/2.0;
        double[] k1 = mProblem.getDerivative(time, state);
        for (int i=0; i<state.length; i++) {
            auxState[i] = state[i] + h2 * k1[i];
        }
        double[] k2 = mProblem.getDerivative(time+h2, auxState);
        for (int i=0; i<state.length; i++) {
            auxState[i] = state[i] + h2 * k2[i];
        }
        double[] k3 = mProblem.getDerivative(time+h2, auxState);
        for (int i=0; i<state.length; i++) {
            auxState[i] = state[i] + deltaTime * k3[i];
        }
        double[] k4 = mProblem.getDerivative(time+deltaTime, auxState);
        double h6 = deltaTime/6;
        for (int i=0; i<state.length; i++) {
            state[i] += h6 * (k1[i]+2*k2[i]+2*k3[i]+k4[i]);
        }
        return time+deltaTime;

	}

}
