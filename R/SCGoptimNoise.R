# SCGOPTIMNOISE perform Scaled conjugate gradient optimization with noise parameter
# CHIPDYNO toolbox
# SCGoptimNoise.R version 1.0.1
# FORMAT SCGoptimNoise <- SCGoptimNoise(params, options, data, precs, X, nEffectGenes, R, C)
# DESC Optimize the gradient based on Scaled conjugate gradient described as in Netlab 
# ARG params: concatenated vector of multiple parameters(beta, gamma, 
# initial mean of the transcription factors, and 
# a vector to create diagonal matrix used to reduce the sparsity of covariance)
# ARG options : 
#	options[1] is set to 1 to display error values; also logs error
#	values in the return argument ERRLOG, and the points visited in the
#	return argument POINTSLOG.  If options[1] is set to 0, then only
#	warning messages are displayed.  If options[1] is -1, then nothing is
#	displayed.
#	options[2] is a measure of the absolute precision required for the
#	value of X at the solution.  If the absolute difference between the
#	values of X between two successive steps is less than options[2],
#	then this condition is satisfied.
#	options[3] is a measure of the precision required of the objective
#	function at the solution.  If the absolute difference between the
#	objective function values between two successive steps is less than
#	options[3], then this condition is satisfied. Both this and the
#	previous condition must be satisfied for termination.
#	options[9] is set to 1 to check the user defined gradient function.
#	options[10] returns the total number of function evaluations
#	(including those in any line searches).
#	options[11] returns the total number of gradient evaluations.
#	options[14] is the maximum number of iterations; default 100.
# ARG data : point estimate of the expression level
# ARG precs : uncertainty of the expression level
# ARG X : connectivity measurement between genes and transcription factors
# ARG nEffectGenes : Number of effectice genes
# ARG R, C : same length integer vectors specifying the row and column 
# indices of the non-zero entries of the sparce matrix
# RETURN x : Optimized value of the parameter 'params'
# COPYRIGHT : Neil D. Lawrence, 2006
# COPYRIGHT : Guido Sanguinetti, 2006
# MODIFICATIONS : Muhammad A. Rahman, 2013
# SEEALSO : SCGoptim
#
# Copyright (c) Ian T Nabney (1996-2001)
# COPYRIGHT : Ian T Nabney, 1996-2001, (the matlab version)
# MODIFICATIONS : Muhammad A. Rahman, 2013, (the R version)

SCGoptimNoise <- function(x, options, data, precs, X, nEffectGenes, R, C) {

eps= 2.2204e-16

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
fold = chipDynoLikeStatNoise(x,data,precs, X,nEffectGenes,R,C);	# Initial function value.
fnow = fold;
options[10] = options[10] + 1;		# Increment function evaluation counter.


gradnew = chipDynoLikeStatNoiseGrad(x, data, precs, X, nEffectGenes, R, C);	# Initial gradient.

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

while (j <= niters){
	if (success==1){
	mu = d %*% t(gradnew)
		if (mu >= 0){
			d= -gradnew
			mu = d %*% t(gradnew)
		}
	kappa = d%*%t(d)

		if (kappa < eps){
			options[8] = fnow;
			return(x) ##?? gradnew?
		}

	sigma = sigma0/sqrt(kappa)
	xplus = x + sigma%*%d
	gplus = chipDynoLikeStatNoiseGrad(xplus, data, precs, X, nEffectGenes, R, C)
	
	options[11]= options[11] + 1
	theta = (d%*%(t(gplus) - t(gradnew)))/sigma
	}

#  % Increase effective curvature and evaluate step size alpha.
	delta = theta + beta*kappa
	if (delta <= 0){ 
		delta = beta*kappa;
		beta = beta - theta/kappa;
	}
	alpha = - mu/delta;

#  % Increase effective curvature and evaluate step size alpha.

  	#% Calculate the comparison ratio.
	xnew = x + alpha %*% d;
	fnew = chipDynoLikeStatNoise(xnew, data, precs, X, nEffectGenes, R, C)
	options[10] = options[10] + 1;
	Delta = 2*(fnew - fold)/(alpha*mu);


#  Calculate the comparison ratio.

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

	if (display > 0){

		cat("Cycle :", j, "\tError= ",  fnow, "\tScale= ", beta, "\n");
	}

	if (success == 1){

		if (max(abs(alpha%*%d)) < options[2] & max(abs(fnew-fold)) < options[3]){
			options[8] = fnew;

			return(x);

		} else {
			fold = fnew;
			gradold = gradnew;

			gradnew = chipDynoLikeStatNoiseGrad(x, data, precs, X, nEffectGenes, R, C)

			options[11] = options[11] + 1

#      % Update variables for new position

			#      % If the gradient is zero then we are done.
			if (gradnew %*% t(gradnew) == 0){
				options[8] = fnew;

				return(x);
			}
		}
	}

#      % If the gradient is zero then we are done.

	#  % Adjust beta according to comparison ratio.
	if (Delta < 0.25){
		beta = min(4.0*beta, betamax);
		}
	if (Delta > 0.75){
		beta = max(0.5*beta, betamin);
		}

	#  % Update search direction using Polak-Ribiere formula, or re-start 
	#  % in direction of negative gradient after nparams steps.
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

#% If we get here, then we haven't terminated in the given number of 
#% iterations.

options[8] = fold;
if (options[1] >= 0){
	print('Warning: Maximum number of iterations has been exceeded');
}

#% If we get here, then we haven't terminated in the given number of 
#% iterations.

return(x)
}
