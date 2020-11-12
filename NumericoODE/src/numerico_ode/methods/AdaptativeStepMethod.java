package numerico_ode.methods;

import java.util.ArrayList;
import java.util.Iterator;

import javax.swing.JFrame;

import org.opensourcephysics.display.Dataset;
import org.opensourcephysics.frames.PlotFrame;

import numerico_ode.ode.InitialValueProblem;
import numerico_ode.ode.NumericalSolution;
import numerico_ode.ode.NumericalSolutionPoint;

public abstract class AdaptativeStepMethod {
	/**
	 * Abstract class for a Fixed Step Method to solve an InitialValueProblem
	 * 
	 * @author Jose Luis
	 * @version November 2020
	 */
	 
	    protected InitialValueProblem mProblem;
	    protected double mStep;
	    protected double mTol;
	    protected double maxStep;
	    protected double minStep;
	    protected NumericalSolution mSolution;
	    protected long mEvaluationCounter = 0;
	    protected ArrayList<Double> steps;
	    protected ArrayList<Double> timesStepChanges;
	    
	    /**
	     * Initializes the method for a given InitialValueProblem
	     * @param InitialValueProblem problem 
	     * @param step the fixed step to take. If negative, we'd solve backwards in time
	     */
	    protected AdaptativeStepMethod(InitialValueProblem problem, double step, double tol, double maxStep, double minStep) {
	        mProblem = problem;
	        mStep = step;
	        mTol = tol;
	        this.minStep = minStep;
	        this.maxStep = maxStep;
	        mSolution = new NumericalSolution(problem);
	        steps = new ArrayList<>();
	        steps.add(mStep);
	        timesStepChanges = new ArrayList<>();
	        timesStepChanges.add(0.0);
	    }
	    
	    /**
	     * Particular method implementation
	     * @param deltaTime the step to take
	     * @param time the current time
	     * @param state the current state
	     * @return the value of time of the step taken, state will contain the updated state
	     */
	    abstract protected double doStep(double time, double[] state);
	    
	    /**
	     * Steps the problem once
	     * @return the newly computed solution point
	     */
	    public NumericalSolutionPoint step() {
	        NumericalSolutionPoint lastPoint = mSolution.getLastPoint();
	        double time = lastPoint.getTime();
	        double[] state = lastPoint.getState();
	        time = doStep(time,state);
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
	                time = doStep(time,state);
	                if (Double.isNaN(time)) return Double.NaN;
	                mSolution.add(time, state);
	            }
	        } 
	        else if (mStep<0) {
	            while (time>finalTime) {
	                time = doStep(time,state);
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
	    
	    public void stepPlot() {
	        PlotFrame frame = new PlotFrame ( "Time" ,"Step" , "Step change frame") ;
	        frame.setConnected(true); // sets default to connect dataset points
	        frame.setSize(800,600);
	        frame.setMarkerShape(0, Dataset.NO_MARKER);
	        frame.setMarkerColor(0,java.awt.Color.BLUE);
	        frame.setXYColumnNames (0 , "Time" , "Step") ; // sets names for each dataset
	        Iterator<Double> iterator = steps.iterator();
	        Iterator<Double> iterator2 = timesStepChanges.iterator();
	        double point=0;
	        double time = 0;
	        while (iterator.hasNext()) {
	            point = iterator.next();
	            time = iterator2.next();
	            frame.append(0, time, point);
	        }
	        frame.setVisible ( true ) ;
	        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE) ;
	    }
	    
}
