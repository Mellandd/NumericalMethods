package numerico_ode.methods;

import java.util.stream.IntStream;

import numerico_ode.ode.InitialValueProblem;

public class AdaptativeStepRKF45Method extends AdaptativeStepMethod{
	
	public AdaptativeStepRKF45Method(InitialValueProblem problem, double step, double tol, double maxStep, double minStep) {
		super(problem, step, tol, maxStep, minStep);
	}
	
	@Override
	protected double doStep(double time, double[] state) {
		Boolean accepted = false;
		double[] stateRK4 = new double[state.length];
		double[] stateRK5 = new double[state.length];
		while(!accepted) {
			double[] k1 = mProblem.getDerivative(time, state);
	        double[] stateAux = new double[state.length];
	        // We use Streams, from Java 8 and above, to do the vectorial operation
	        IntStream.range(0, state.length).forEach(i-> stateAux[i] = state[i]+ 0.25*k1[i]);
			double[] k2 = mProblem.getDerivative(time + mStep/4., stateAux);
	        IntStream.range(0, state.length).forEach(i-> stateAux[i] = state[i]+ (3./32.)*k1[i]+(9./32.)*k2[i]);
			double[] k3 = mProblem.getDerivative(time + (3./8.)*mStep, stateAux);
	        IntStream.range(0, state.length).forEach(i-> stateAux[i] = state[i]+ (1932./2197.)*k1[i]-(7200./2197.)*k2[i]+(7296./2197.)*k3[i]);
			double[] k4 = mProblem.getDerivative(time + (12./13.)*mStep, stateAux);
	        IntStream.range(0, state.length).forEach(i-> stateAux[i] = state[i]+ (439./216.)*k1[i]-8*k2[i]+(3680./513.)*k3[i]-(845./4104.)*k4[i]);
			double[] k5 = mProblem.getDerivative(time + mStep, stateAux);
	        IntStream.range(0, state.length).forEach(i-> stateAux[i] = state[i]-(8./27.)*k1[i]+2*k2[i]-(3544./2565.)*k3[i]+(1859./4104.)*k4[i]-(11./40.)*k5[i]);
			double[] k6 = mProblem.getDerivative(time + mStep/2., stateAux);
			IntStream.range(0, state.length).forEach(i->stateRK4[i] = state[i]+(25./216.)*k1[i]+(1408./2565.)*k3[i]+(2197./4101.)*k4[i]-(1./5.)*k5[i]);
			IntStream.range(0, state.length).forEach(i->stateRK5[i] = state[i]+(16./135.)*k1[i]+(6656./12825.)*k3[i]
					+(28561./56430.)*k4[i]-(9./50.)*k5[i]+(2./55.)*k6[i]);
			//Calculate the error
			double normError = 0;
			for (int j = 0; j < state.length; j++) {
				normError+=Math.pow(stateRK5[j]-stateRK4[j],2.);
			}
			normError = (16./15.)*Math.sqrt(normError);
			if (normError < mTol*mStep) {
				accepted=true;
			}
			double q = Math.pow((mTol*mStep)/(2*normError),1./4.);
			q = Math.min(4, Math.max(q, 0.1));
			mStep = q* mStep;
			if (mStep < minStep) return -1;
			if (mStep > maxStep) mStep = maxStep;
		}
			IntStream.range(0, state.length).forEach(i->state[i]=(16*stateRK5[i]-stateRK4[i])/15.);
			return time + mStep;
	}

}
