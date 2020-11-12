package numerico_ode.methods;

import numerico_ode.ode.InitialValueProblem;

public class AdaptativeStepRK4Method extends AdaptativeStepMethod{
	
	public AdaptativeStepRK4Method(InitialValueProblem problem, double step, double tol, double maxStep, double minStep) {
		super(problem, step, tol, maxStep, minStep);
	}
	
	@Override
	protected double doStep(double time, double[] state) {
		Boolean accepted = false;
		double[] rk41 = new double[state.length];
		double[] rk42 = new double[state.length];
		double[] rk43 = new double[state.length];
		double[] k1 = new double[state.length];
		double[] k2 = new double[state.length];
		double[] k3 = new double[state.length];
		double[] k4 = new double[state.length];
		double[] auxState = new double[state.length];
		double q, normError, h2;
		while(!accepted) {
			this.addToEvaluationCounter(11);			
			/**
			 * CALCULAMOS LOS 11 K
			 */
			h2 = mStep/2.;
	        k1 = mProblem.getDerivative(time, state);
	        for (int i=0; i<state.length; i++) 
	            auxState[i] = state[i] + h2 * k1[i];
	        k2 = mProblem.getDerivative(time+h2, auxState);
	        for (int i=0; i<state.length; i++) 
	            auxState[i] = state[i] + h2 * k2[i];
	        k3 = mProblem.getDerivative(time+h2, auxState);
	        for (int i=0; i<state.length; i++) 
	            auxState[i] = state[i] + mStep * k3[i];
	        k4 = mProblem.getDerivative(time+mStep, auxState);
	        double h6 = mStep/6;
	        for (int i=0; i<state.length; i++) 
	            rk41[i] = state[i]+ h6*(k1[i]+2*k2[i]+2*k3[i]+k4[i]);
	        double h = mStep/2.;
	        h2 = mStep/4.;
	        h6 = mStep/12.;
			
	        for (int i=0; i<state.length; i++) 
	            auxState[i] = state[i] + h2 * k1[i];
	        k2 = mProblem.getDerivative(time+h2, auxState);
	        for (int i=0; i<state.length; i++) 
	            auxState[i] = state[i] + h2 * k2[i];
	        k3 = mProblem.getDerivative(time+h2, auxState);
	        for (int i=0; i<state.length; i++) 
	            auxState[i] = state[i] + h * k3[i];
	        k4 = mProblem.getDerivative(time+h, auxState);
	        for (int i=0; i<state.length; i++) 
	            rk42[i] =state[i]+ h6*(k1[i]+2*k2[i]+2*k3[i]+k4[i]);
			
	        k1 = mProblem.getDerivative(time + h, rk42);
	        for (int i=0; i<state.length; i++) 
	            auxState[i] = rk42[i] + h2 * k1[i];
	        k2 = mProblem.getDerivative(time+h+h2, auxState);
	        for (int i=0; i<state.length; i++) 
	            auxState[i] = rk42[i] + h2 * k2[i];
	        k3 = mProblem.getDerivative(time+h+h2, auxState);
	        for (int i=0; i<state.length; i++) 
	            auxState[i] = rk42[i] + h* k3[i];
	        k4 = mProblem.getDerivative(time+mStep, auxState);
	        for (int i=0; i<state.length; i++) 
	            rk43[i] = rk42[i]+ h6*(k1[i]+2*k2[i]+2*k3[i]+k4[i]);

			//Calculate the error
			normError = 0;
			for (int j = 0; j < state.length; j++) {
				normError+=Math.pow(rk43[j]-rk41[j],2.);
			}
			normError = (16./15.)*Math.sqrt(normError);
			if (normError < mTol*mStep) {
				accepted=true;
			}
			q = Math.pow((mTol*mStep)/(2*normError),1./4.);
			q = Math.min(4, Math.max(q, 0.1));
			mStep = q* mStep;
			if (mStep < minStep) {
				return -1;
			}
			if (mStep > maxStep) mStep = maxStep;
			steps.add(mStep);
			timesStepChanges.add(time);
		}
		for (int j = 0; j < rk43.length; j++) {
			state[j]=(16.*rk43[j]-rk41[j])/15.;
		}
		return time + mStep;
	}

}
