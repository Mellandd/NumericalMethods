/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package numerico_ode.methods;

import numerico_ode.ode.InitialValueProblem;
import numerico_ode.ode.NumericalSolution;
import numerico_ode.ode.NumericalSolutionPoint;

/**
 * Abstract class for a Fixed Step Method to solve an InitialValueProblem
 * 
 * @author F. Esquembre
 * @version September 2020
 */
abstract public class FixedStepMethod {
 
    protected InitialValueProblem mProblem;
    protected double mStep;
    protected NumericalSolution mSolution;
    protected long mEvaluationCounter = 0;
    
    /**
     * Initializes the method for a given InitialValueProblem
     * @param InitialValueProblem problem 
     * @param step the fixed step to take. If negative, we'd solve backwards in time
     */
    public FixedStepMethod(InitialValueProblem problem, double step) {
        mProblem = problem;
        mStep = step;
        mSolution = new NumericalSolution(problem);
    }
    
    /**
     * Particular method implementation
     * @param deltaTime the step to take
     * @param time the current time
     * @param state the current state
     * @return the value of time of the step taken, state will contain the updated state
     */
    abstract public double doStep(double deltaTime, double time, double[] state);
    
    /**
     * Steps the problem once
     * @return the newly computed solution point
     */
    public NumericalSolutionPoint step() {
        NumericalSolutionPoint lastPoint = mSolution.getLastPoint();
        double time = lastPoint.getTime();
        double[] state = lastPoint.getState();
        time = doStep(mStep,time,state);
        if (Double.isNaN(time)) return null;
        return mSolution.add(time, state);
    }
    
    /**
     * Iteratively steps the problem until time equals or exceeds finalTime
     * @param finalTime the time which we want to reach or exceed
     * @return the actual time of the last computed solution point (may differ -exceed- the requested finalTime)
     */
    public double solve(double finalTime) {
        NumericalSolutionPoint lastPoint = mSolution.getLastPoint();
        double time = lastPoint.getTime();
        double[] state = lastPoint.getState();
        if (mStep>0) {
            while (time<finalTime) {
                time = doStep(mStep,time,state);
                if (Double.isNaN(time)) return Double.NaN;
                mSolution.add(time, state);
            }
        } 
        else if (mStep<0) {
            while (time>finalTime) {
                time = doStep(mStep,time,state);
                if (Double.isNaN(time)) return Double.NaN;
                mSolution.add(time, state);
            }
        } // does nothing if mStep = 0
        return time;
    }
    
    /**
     * Gets the solution computed so far
     * @return an instance of NumericalSolution
     */
    public NumericalSolution getSolution() { return mSolution; }
    
    public void resetEvaluationCounter() { 
        mEvaluationCounter = 0;
    }
    
    public long getEvaluationCounter() {
        return mEvaluationCounter;
    }

    protected void addToEvaluationCounter(int add) {
        mEvaluationCounter += add;
    }
    

}
