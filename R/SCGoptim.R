
SCGoptim = function( x, options, data, X, nEffectGenes, R, C) {

#function [x, options, flog, pointlog, scalelog] = scg(f, x, options, gradf, varargin)
#%SCG	Scaled conjugate gradient optimization.
#%
#%	Description
#%	[X, OPTIONS] = SCG(F, X, OPTIONS, GRADF) uses a scaled conjugate
#%	gradients algorithm to find a local minimum of the function F(X)
#%	whose gradient is given by GRADF(X).  Here X is a row vector and F
#%	returns a scalar value. The point at which F has a local minimum is
#%	returned as X.  The function value at that point is returned in
#%	OPTIONS(8).
#%
#%	[X, OPTIONS, FLOG, POINTLOG, SCALELOG] = SCG(F, X, OPTIONS, GRADF)
#%	also returns (optionally) a log of the function values after each
#%	cycle in FLOG, a log of the points visited in POINTLOG, and a log of
#%	the scale values in the algorithm in SCALELOG.
#%
#%	SCG(F, X, OPTIONS, GRADF, P1, P2, ...) allows additional arguments to
#%	be passed to F() and GRADF().     The optional parameters have the
#%	following interpretations.
#%
#%	OPTIONS(1) is set to 1 to display error values; also logs error
#%	values in the return argument ERRLOG, and the points visited in the
#%	return argument POINTSLOG.  If OPTIONS(1) is set to 0, then only
#%	warning messages are displayed.  If OPTIONS(1) is -1, then nothing is
#%	displayed.
#%
#%	OPTIONS(2) is a measure of the absolute precision required for the
#%	value of X at the solution.  If the absolute difference between the
#%	values of X between two successive steps is less than OPTIONS(2),
#%	then this condition is satisfied.
#%
#%	OPTIONS(3) is a measure of the precision required of the objective
#%	function at the solution.  If the absolute difference between the
#%	objective function values between two successive steps is less than
#%	OPTIONS(3), then this condition is satisfied. Both this and the
#%	previous condition must be satisfied for termination.
#%
#%	OPTIONS(9) is set to 1 to check the user defined gradient function.
#%
#%	OPTIONS(10) returns the total number of function evaluations
#%	(including those in any line searches).
#%
#%	OPTIONS(11) returns the total number of gradient evaluations.
#%
#%	OPTIONS(14) is the maximum number of iterations; default 100.
#%
#%	See also
#%	CONJGRAD, QUASINEW
#%

#%	Copyright (c) Ian T Nabney (1996-2001)


eps= 2.2204e-16

#  Set up the options.
if (length(options) < 18){
  stop('Options vector too short')
}

if(options[14]!=0){
	niters = options[14];
} else {
	niters = 10;
}


display = options[1];
gradcheck = options[9];

nparams = length(x)

sigma0 = 1.0e-4;
fold = chipDynoLikeStat(x,data, X,nEffectGenes,R,C);	# Initial function value.
fnow = fold;
options[10] = options[10] + 1;		# Increment function evaluation counter.


gradnew = chipDynoLikeStatGrad(x, data, X, nEffectGenes, R, C);	# Initial gradient.

gradold = gradnew;
options[11] = options[11] + 1;		# Increment gradient evaluation counter.

d = -gradnew;				# Initial search direction.
success = 1;				# Force calculation of directional derivs.
nsuccess = 0;				# nsuccess counts number of successes.
beta = 1.0;				# Initial scale parameter.
betamin = 1.0e-15; 			# Lower bound on scale.
betamax = 1.0e100;			# Upper bound on scale.
j = 1;					# j counts number of iterations.

flog = array(0, dim=c(niters,1))
pointlog = array(0, dim=c(niters,length(params)))
scalelog = array(0, dim=c(niters,1))

flog[j,] = fold
pointlog[j,] = x;


cat("SCG optmization is running- This process will require some hours!\n");

while (j <= niters){

	if (success==1){
	mu = d %*% t(gradnew)
		if (mu >= 0){
			d= -gradnew
			mu = d %*% t(gradnew)
		}
	kappa = d%*%t(d)
	# eps= 2.2204e-16
		if (kappa < eps){
			options[8] = fnow;
			return(x) ##??
		}

	sigma = sigma0/sqrt(kappa)
	xplus = x + sigma%*%d

	gplus = chipDynoLikeStatGrad(xplus, data, X, nEffectGenes, R, C)
	
	options[11]= options[11] + 1
	theta = (d%*%(t(gplus) - t(gradnew)))/sigma
	}

	delta = theta + beta*kappa
	if (delta <= 0){ 
		delta = beta*kappa;
		beta = beta - theta/kappa;
	}
	alpha = - mu/delta;
	xnew = x + alpha %*% d;
	fnew = chipDynoLikeStat(xnew, data, X, nEffectGenes, R, C)
	options[10] = options[10] + 1;
	Delta = 2*(fnew - fold)/(alpha*mu);

	if (Delta  >= 0){
		success = 1;
		nsuccess = nsuccess + 1;
		x = xnew;
		fnow = fnew;
	} else {
		success = 0;
		fnow = fold;
	}

	flog[j] = fnow
	pointlog[j,] = x
	scalelog[j] = beta
	####

	if (display > 0){
		#sprintf("Cycle %4d  Error %11.6f  Scale %e", j, fnow, beta);
		cat("Cycle :", j, "\tError= ",  fnow, "\tScale= ", beta, "\n");
	}
    
	if (success == 1){

		if (max(abs(alpha%*%d)) < options[2] & max(abs(fnew-fold)) < options[3]){
			options[8] = fnew;
			return(x);

		} else {
			fold = fnew;
			gradold = gradnew;
			gradnew = chipDynoLikeStatGrad(x, data, X, nEffectGenes, R, C)
			options[11] = options[11] + 1

			if (gradnew %*% t(gradnew) == 0){
				options[8] = fnew;
				return(x);
			}
		}
	}

	if (Delta < 0.25){
		beta = min(4.0*beta, betamax);
		}
	if (Delta > 0.75){
		beta = max(0.5*beta, betamin);
		}

	if (nsuccess == nparams){
		d = -gradnew;
		nsuccess = 0;
	} else {
		if (success == 1){
			gamma = (gradold - gradnew)%*%t(gradnew)/(mu);
			d = gamma %*% d - gradnew;
		}
	}	

	j = j + 1;
}

options[8] = fold;
if (options[1] >= 0){
	print('Warning: Maximum number of iterations has been exceeded');
}

return(x)
}
