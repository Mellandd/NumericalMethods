package numerico_ode.methods;
import java.util.stream.IntStream;

import numerico_ode.ode.InitialValueProblem;

public class FixedStepModifiedEulerMethod extends FixedStepMethod{
    
    /**
     * Initializes the method for a given InitialValueProblem
     * @param InitialValueProblem problem 
     * @param step the fixed step to take. If negative, we'd solve backwards in time
     */
    public FixedStepModifiedEulerMethod(InitialValueProblem problem, double step) {
        super(problem,step);
    }

    
    @Override
    /**
     * Euler method implementation modified to do 2 evaluations
     * @param deltaTime the step to take
     * @param time the current time
     * @param state the current state
     * @return the value of time of the step taken, state will contain the updated state
     */
    public double doStep(double deltaTime, double time, double[] state) {
        double[] derivative = mProblem.getDerivative(time, state);
        double[] state2 = new double[state.length];
        // We use Streams, from Java 8 and above, to do the vectorial operation
        // since we need state + 1/2 k1 to do k2.
        IntStream.range(0, state.length).forEach(i-> state2[i] = state[i]+ deltaTime*derivative[i]);
        double[] derivative2 = mProblem.getDerivative(time + deltaTime, state2);
        super.addToEvaluationCounter(2);
        for (int i=0; i<state.length; i++) {
            state[i] = state[i] + (deltaTime/2.) * (derivative[i]+derivative2[i]);	
        }
        return time+deltaTime;
    }
    
}