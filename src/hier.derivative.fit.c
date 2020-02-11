/*************************************************************
 * Copyright (c) 2012 Garritt Leland Page
 *
 * This file contains C code for an MCMC algorithm constructed
 * to fit a hierarchical q-degree (b)-spline model to a
 * constructed sequence of empirical derivatives.  A  fixed
 * number of knots is used and penalization is carried out
 * using P-splines.
 *
 *************************************************************/

#include "matrix.h"
#include "Rutil.h"

#include <R_ext/Lapack.h>
#include <R.h>
#include <Rmath.h>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/*****************************************************************************************
* The following are the inputs of the function that are read from R
*
* draws = total number of MCMC draws
* burn = number of MCMC draws discarded as burn-in
* thin = indicates how much MCMC chain should be thinned
*
* nsubject = number of subjects in data set
* nobs = vector whose entries indicate number of observations per subject
* y = vector containing nobs production values for each player
* tm = vector containing nobs time (time is incremented by one) values for each player
* nb = scalar indication ncolumns in H matrix
* K = smoothing matrix for penalized B-splines.
* gvec = nsubject x 1 vector indicating to which group subject belongs
* ng = integer indicating the number of groups
*
* Aparm = scalar associated with upper bound prior on beta_j
*
* modelPriors = vector containing prior values from the hierarchical model
*
* p_lag = Prior distribution for D
* Hmat
* Htheta
* nrHt = scalar indication number of rows in Htheta matrix
* balanced = scalar indicating if design in balanced or not
* ztype = integer indicating what sequence of quotient differences to use
	0 - Mine (employing both D and u)
	1 - Cote (employing only u)
	2 - De Brabantar.
* verbose = integer determining whether to print to screen information regarding run.
*****************************************************************************************/


void mcmcloop(int *draws, int *burn, int *thin, int *nsubject, int *nobs, double *y,
			        double *tm, int *nb, double *K, int *gvec, int *ng, double *Aparm,
			        double *p_lag, double *au, double *bu, double *modelPriors, double *Hmat,
			        double* Htheta, int *nrHt, int *balanced, int *ztype, int *verbose,
			        double *beta, double *theta, double *lam, double *tau2, double *u, int *D,
			        double *zout, double *sig2, double *fprime, double *fgprime){




	// i - MCMC iterate
	// ii - MCMC iterate that is saved
	// j - player iterate
	// jj - second player iterate
	// t - time (cumulative minutes played) per player iterate
	// b - beta and thetah iterate
	// bb - second beta and thetah iterate
	// k - cluster iterate
	// kk - knot iterate

	int i, j, t, jj, k, kk, b, bb, h;
	int ii = 0;
	int csobs = 0, csHobs=0;
	int N ;


	int nout = (*draws - *burn)/(*thin);


	// create max_obs and n_group vectors etc.
	int  max_nobs;
	int n_group[*ng];
	for(k = 0; k < *ng; k++) n_group[k] = 0;

	double *max_tm = R_Vector(*nsubject);
	double *min_tm = R_Vector(*nsubject);

	double *ypy = R_VectorInit(*nsubject, 0.0);
	double *sumy = R_VectorInit(*nsubject,0.0);

	max_nobs = nobs[0];
	csobs = 0;
	for(j = 0; j < (*nsubject); j++){
		max_tm[j] = tm[csobs];
		min_tm[j] = tm[csobs];

		for(t = 0; t < nobs[j]; t++){

			ypy[j] = ypy[j] + y[csobs]*y[csobs];
			sumy[j] = sumy[j] + y[csobs];

			if(max_tm[j] < tm[csobs]) max_tm[j] = tm[csobs];
			if(min_tm[j] > tm[csobs]) min_tm[j] = tm[csobs];

			csobs = csobs + 1;

		}

		if(max_nobs < nobs[j]) max_nobs = nobs[j];

		n_group[gvec[j]-1] = n_group[gvec[j]-1] + 1;

	}


	N = csobs;

	double max_D;
	max_D = ceil((double) max_nobs/2.0);


	if(*verbose){
		Rprintf("nsubject = %d\n", *nsubject);
		Rprintf("nb = %d\n", *nb);
		Rprintf("nout = %d\n", nout);

		RprintIVecAsMat("nobs", nobs, 1, *nsubject);

		RprintIVecAsMat("n_group" , n_group, 1, *ng);
		Rprintf("max_nobs = %d\n", max_nobs);
		Rprintf("N = %d\n", N);
		Rprintf("max_D = %f\n", max_D);
	}








	// ===================================================================================
	//
	// Memory vectors to hold MCMC iterates for non cluster specific parameters
	//
	// ===================================================================================


	double *sig2_iter = R_VectorInit(*nsubject, 1.0);

	double *beta_iter = R_VectorInit((*nb)*(*nsubject),0.0);

	double *u_iter = R_VectorInit(*nsubject, 1.0);
	double *z_iter = R_VectorInit((*nsubject)*max_nobs, 0.0);
	double *fprime_iter = R_VectorInit((*nsubject)*max_nobs, 0.0);
	double *fgprime_iter = R_VectorInit((*ng)*(*nrHt), 0.0);

	int D_iter[*nsubject]; for(j=0;j<*nsubject;j++)D_iter[j]=1;



	// ===================================================================================
	//
	// Memory vectors to hold MCMC iterates for cluster specific parameters
	//
	// ===================================================================================

	double *tau2h = R_VectorInit(*ng, 1.0);
	double *lamh = R_VectorInit(*ng, 0.5*(*Aparm));
	double *thetah = R_VectorInit((*nb)*(*ng), 0.0);




	// ===================================================================================
	//
	// scratch vectors of memory needed to update parameters
	//
	// ===================================================================================

	// These are made particularly big to make sure there is enough memory
	double *scr1 = R_Vector(((*nb)*(*nb)));
	double *scr2 = R_Vector(((*nb)*(*nb)));
	double *scr3 = R_Vector(((*nb)*max_nobs));

	// stuff needed to update u and/or D;
	int Do, Dn, Dup;
	double uo, un, pidtmp, totpid;
	double *zo = R_VectorInit(max_nobs*2,0.0);
	double *zn = R_VectorInit(max_nobs*2,0.0);
	double *ytmp = R_VectorInit(max_nobs*2,0.0);
	double *ttmp = R_VectorInit(max_nobs*2,0.0);
	double *pid = R_VectorInit(max_nobs*2,1.0);


	// stuff I need to update sig2 (player specific), sig2b0;
	double astar, bstar, sumsq;

	// stuff that I need for player specific beta's;
	double *H = R_VectorInit((*nsubject)*((*nb))*(max_nobs),0.0);
	double *tH = R_VectorInit((*nsubject)*((*nb))*(max_nobs),0.0);
	double *HtH = R_Vector(((*nb)*(*nb)));
	double *Htz = R_Vector(((*nb)*max_nobs));

	double *z_tmp = R_VectorInit(max_nobs*2,0.0);

	double sumz_Hb;
	double *z_b0 = R_Vector(max_nobs*2);
	double *Hb = R_Vector(max_nobs*2);
	double *Ht = R_Vector((*nrHt));


	// Create the inverse of the K penalty matrix
	double ld;
	double *Kinv = R_Vector((*nb)*(*nb));
	for(b = 0; b < (*nb); b++){for(bb = 0; bb < (*nb); bb++){Kinv[b*(*nb)+bb] = K[b*(*nb)+bb];}}
	cholesky(Kinv, (*nb), &ld);
	inverse_from_cholesky(Kinv, scr1, scr2, (*nb));

	// stuff I need for tau2h;
	double *thtmp = R_Vector((*nb));

	// stuff I need for thetah
	double *sumbeta = R_Vector((*nb));
	double *Mstar = R_Vector((*nb));
	double *Sstar = R_Vector((*nb)*(*nb));
	double *outrmvnorm = R_Vector((*nb));

	//stuff I need for lamh
	double olam, nlam, lln, llo, ldo, ldn, llr, uu;
	double *btmp = R_Vector((*nb)*(*nsubject));
	double *nV = R_Vector((*nb)*(*nb));
	double *oV = R_Vector((*nb)*(*nb));



	// Recall Hmat is read in as a long contiguous array,
	// so care must be taken to read in the FORTRAN format.
	if(*balanced==1){
		for(t=0; t < nobs[0]; t++){
			for(kk=0; kk<(*nb); kk++){
				H[t*((*nb)) + kk] = Hmat[kk*((nobs[0])) + t];
			}
		}
	}



	// ===================================================================================
	//
	// Prior parameter values
	//
	// ===================================================================================

	// prior values for sig2
	double asig = modelPriors[0];
	double bsig = modelPriors[1];


	// IG parameters for tau2
	double at = modelPriors[2]; double bt = modelPriors[3];

	double a_u = *au;
	double b_u = *bu;

	// Uniform prior for lam
	double A = *Aparm;

	// Candidate density sd's for M-H stuff
	double csigLAM = 0.5, csigU=0.5;


	if(*verbose){

		Rprintf("at = %f\n", at);
		Rprintf("bt = %f\n", bt);

		Rprintf("A = %f\n", A);
		Rprintf("au = %f\n", a_u);
		Rprintf("bu = %f\n", b_u);
		Rprintf("p = %f\n", *p_lag);


	}


	GetRNGstate();


	// ===================================================================================
	//
	// start of the mcmc algorithm;
	//
	// ===================================================================================

	for(i = 0; i < *draws; i++){

		if(*verbose){
			if((i+1) % 10000 == 0){
				time_t now;
				time(&now);

				Rprintf("mcmc iter = %d =========================================== \n", i+1);
				Rprintf("%s", ctime(&now));

			}
		}


		//////////////////////////////////////////////////////////////////////////////////
		//
		// Update parameters associated with individual within on for-loop.  These include
		// transformation from Y to Z (i.e., u and D), the beta vector, data variance
		// (sig2), and intercept.
		//
		//////////////////////////////////////////////////////////////////////////////////
		csobs=0;
		csHobs=0;
		for(j = 0; j < *nsubject; j++){


			for(t = 0; t < nobs[j]; t++){
				ytmp[t] = y[csobs];
				ttmp[t] = tm[csobs];
				csobs = csobs+1;
			}



			if(*balanced!=1){

				for(t=0; t < nobs[j]*(*nb); t++){
					H[t] = Hmat[csHobs + t];
				}
				csHobs = csHobs + nobs[j]*(*nb);

			}



			for(b = 0; b < (*nb); b++){

				btmp[b] = beta_iter[b*(*nsubject) + j];

			}

			matrix_product(H, btmp, Hb, nobs[j], 1, (*nb));


			///////////////////////////////////////////////////////
			//									   				 //
			// update u and transformation within the same loop.  //
			//									   				 //
			///////////////////////////////////////////////////////

			if(*ztype != 2){
				uo = u_iter[j];


				un = rnorm(uo, csigU);


				if(un > 0.000001){



					// Recall second to last argument indicates whether we use
					// my transformation (0) or Cote's (1).  I will use Cote's for now

					emp_derivative(ytmp, ttmp, D_iter[j], uo, nobs[j], *ztype, zo);
					emp_derivative(ytmp, ttmp, D_iter[j], un, nobs[j], *ztype, zn);


					llo = 0.0, lln = 0.0;
					for(t=0; t<nobs[j]; t++){
						llo = llo + dnorm(zo[t], Hb[t], sqrt(sig2_iter[j]), 1);
						lln = lln + dnorm(zn[t], Hb[t], sqrt(sig2_iter[j]), 1);

					}


					llo = llo + dgamma(uo,a_u,b_u, 1); // second argument is scale
					lln = lln + dgamma(un,a_u,b_u, 1);


					llr = lln - llo;

					uu = runif(0,1);

					if(llr > log(uu)) u_iter[j] = un;


				}

			}




			///////////////////////////////////////////////////////
			//									   				 //
			// update D and transformation within the same loop.  //
			//									   				 //
			///////////////////////////////////////////////////////


			if(*ztype != 1){
				Do = D_iter[j];


				Dn = ran_discrete_unif(Do-3, Do+3);

				Dup = (int) floor((double) nobs[j]/2);


				if((Dn > 0) & (Dn < Dup)){



					// Recall second to last argument indicates whether we use
					// my transformation (0) or Cote's (1).

					emp_derivative(ytmp, ttmp, Do, u_iter[j], nobs[j], *ztype, zo);

					emp_derivative(ytmp, ttmp, Dn, u_iter[j], nobs[j], *ztype, zn);


					llo = 0.0, lln = 0.0;
					for(t=0; t<nobs[j]; t++){
						llo = llo + dnorm(zo[t], Hb[t], sqrt(sig2_iter[j]), 1);
						lln = lln + dnorm(zn[t], Hb[t], sqrt(sig2_iter[j]), 1);

					}


					// If not balanced, I need to create pid for each subject
					// separately.  This is not too expensive so I do it
					// for each subject

					pid[0] = *p_lag;
					pidtmp = 1;
					totpid = pid[0];

					for(h = 1; h < floor((double) nobs[j]/2); h++){
						pidtmp = pidtmp*(1-(*p_lag));
						pid[h] = (*p_lag)*pidtmp;
						totpid = totpid + pid[h];
					}

					for(h = 0; h < floor((double) nobs[j]/2); h++){
						pid[h] = pid[h]/totpid;
					}

					llo = llo + log(pid[Do-1]);
					lln = lln + log(pid[Dn-1]);

					llr = lln - llo;

					uu = runif(0,1);

					if(llr > log(uu)) D_iter[j] = Dn;

				}

			}

			emp_derivative(ytmp, ttmp, D_iter[j], u_iter[j], nobs[j], *ztype, z_tmp);


			for(t = 0; t < nobs[j]; t++){
				z_iter[j*max_nobs + t] = z_tmp[t];
			}


			/////////////////////////////////////////
			//									   //
			// udate beta within the same loop.    //
			//									   //
			/////////////////////////////////////////
			for(t = 0; t < nobs[j]; t++){
				z_b0[t] = z_tmp[t];
			}




			mat_transpose(H, tH, nobs[j], (*nb));

			matrix_product(tH, H, HtH, (*nb), (*nb), nobs[j]);
			matrix_product(tH, z_tmp, Htz, (*nb), 1, nobs[j]);



			for(b = 0; b < (*nb); b++){
				for(bb = 0; bb < (*nb); bb++){

					if(*nsubject > 1){

						Sstar[b*(*nb)+bb] = (1/sig2_iter[j])*HtH[b*(*nb)+bb];

						if(b == bb){
							Sstar[b*(*nb)+bb] = (1/sig2_iter[j])*HtH[b*(*nb)+bb] +
											 (1/(lamh[gvec[j]-1]*lamh[gvec[j]-1]));
						}
					}

					if(*nsubject==1){

						Sstar[b*(*nb)+bb] = (1/sig2_iter[j])*HtH[b*(*nb)+bb] +
						                    (1/tau2h[0])*K[b*(*nb)+bb];
					}



				}

			}

			cholesky(Sstar, (*nb), &ld);
			inverse_from_cholesky(Sstar, scr1, scr2, (*nb));


			for(b = 0; b < (*nb); b++){
				if(*nsubject > 1){
					scr3[b] = (1/sig2_iter[j])*Htz[b] +
							 (1/(lamh[gvec[j]-1]*lamh[gvec[j]-1]))*
							 thetah[(gvec[j]-1)*((*nb)) + b];
				}

				if(*nsubject == 1){
					scr3[b] = (1/sig2_iter[j])*Htz[b];

				}

			}


			matrix_product(Sstar, scr3, Mstar, (*nb), 1, (*nb));


			cholesky(Sstar, (*nb) , &ld);

			ran_mvnorm(Mstar, Sstar, (*nb), scr1, outrmvnorm);





			for(b = 0; b < (*nb); b++){

				beta_iter[b*(*nsubject) + j] = outrmvnorm[b];
				btmp[b] = beta_iter[b*(*nsubject) + j];

			}


			matrix_product(H, btmp, Hb, nobs[j], 1, (*nb));



			/////////////////////////////////////////
			//									   //
			// udate sigma2 within the same loop.  //
			//									   //
			/////////////////////////////////////////


			sumz_Hb = 0.0;

			for(jj = 0; jj < nobs[j]; jj++){

				scr3[jj] = z_b0[jj] - Hb[jj];
				sumz_Hb = sumz_Hb + (z_tmp[jj] - Hb[jj]);
			}

			sumsq = inner_product(scr3, 1, scr3, 1, nobs[j]);

			astar = 0.5*nobs[j] + asig;
			bstar = 0.5*sumsq + 1/bsig;

			//bstar is rate and rgamma requires scale hence inverse

			sig2_iter[j] = 1/rgamma(astar, 1/bstar);

			for(t = 0; t < nobs[j]; t++){
				fprime_iter[j*max_nobs + t] = Hb[t];
			}


		}







		//////////////////////////////////////////////////////////////////////////////////
		//
		// Update lam (using MH-step) this is the standard deviation not the variance
		//
		//////////////////////////////////////////////////////////////////////////////////
		if(*nsubject > 1){

			for(k = 0; k < *ng; k++){


				olam = lamh[k];
				nlam = rnorm(olam,csigLAM);

				if((nlam > 0) & (nlam < A)){

					for(b = 0; b < (*nb); b++){

						thtmp[b] = thetah[k*(*nb) + b];

						for(bb = 0; bb < (*nb); bb++){

							oV[b*(*nb)+bb] = 0.0;
							nV[b*(*nb)+bb] = 0.0;

							if(b == bb){

								oV[b*(*nb)+bb] = 1/(olam*olam);
								nV[b*(*nb)+bb] = 1/(nlam*nlam);
							}
						}
					}

					ldo = 2.0*(*nb)*log(olam);
					ldn = 2.0*(*nb)*log(nlam);

					lln = 0.0;
					llo = 0.0;
					for(j = 0; j < *nsubject; j++){

						if(gvec[j] == k+1){

							for(b = 0; b < (*nb); b++){

								btmp[b] = beta_iter[b*(*nsubject) + j];

							}

							llo = llo + dmvnorm(btmp, thtmp, oV, (*nb), ldo, scr1, 1);
							lln = lln + dmvnorm(btmp, thtmp, nV, (*nb), ldn, scr1, 1);

						}

					}


					llo = llo + dunif(olam, 0.0, A, 1);
					lln = lln + dunif(nlam, 0.0, A, 1);


					llr = lln - llo;
					uu = runif(0.0,1.0);

					if(log(uu) < llr) lamh[k] = nlam;
				}





				//////////////////////////////////////////////////////////////////////////////
				//																			//
				// udpate thetah each of the cluster specific coefficients;					//
				//																			//
				//////////////////////////////////////////////////////////////////////////////

				for(b = 0; b < (*nb); b++){
					for(bb = 0; bb < (*nb); bb++){

				 		Sstar[b*(*nb)+bb] = (1/tau2h[k])*K[b*(*nb)+bb];

				 		if(b == bb){ Sstar[b*(*nb)+bb] = ((double) n_group[k]/(lamh[k]*lamh[k])) +
				 							  (1/tau2h[k])*K[b*(*nb)+bb];}

					}

					sumbeta[b] = 0.0;

				}

				cholesky(Sstar, (*nb), &ld);
				inverse_from_cholesky(Sstar, scr1, scr2, (*nb));

				for(j = 0; j < *nsubject; j++){

					if(gvec[j] == k+1){

						for(b = 0; b < (*nb); b++){

							sumbeta[b] = sumbeta[b] + (1/(lamh[k]*lamh[k]))*
														beta_iter[b*(*nsubject) + j];
						}
					}
				}


				matrix_product(Sstar, sumbeta, Mstar, (*nb), 1, (*nb));


				cholesky(Sstar, (*nb) , &ld);



				ran_mvnorm(Mstar, Sstar, (*nb), scr1, outrmvnorm);


				for(b = 0; b < (*nb); b++){

					thetah[k*(*nb) + b] = outrmvnorm[b];

				}


 				matrix_product(Htheta, outrmvnorm, Ht, *nrHt, 1, (*nb));


				for(t = 0; t < *nrHt; t++){
					fgprime_iter[k*(*nrHt) + t] = Ht[t];
				}

			}
		}

		//////////////////////////////////////////////////////////////////////////////
		//
		// Update tau2 for each of the clusters (P-spline smoothing parameter)
		//
		//////////////////////////////////////////////////////////////////////////////
		for(k = 0; k < *ng; k++){


			for(b = 0; b < (*nb); b++){

				if(*nsubject > 1){
					thtmp[b] = thetah[k*(*nb) + b];
				}
				if(*nsubject == 1){
					thtmp[b] = beta_iter[b*(*nsubject) + 0];
				}
			}


			sumsq = quform(thtmp,K,(*nb));


			astar = 0.5*((*nb)-2) + at;
			bstar = 1/bt + 0.5*sumsq;


			tau2h[k] = 1/rgamma(astar, 1/bstar);// E(tau2) = astarbstar for gamma.  bstar is scale




		}



		//////////////////////////////////////////////////////////////////////////////////
		//																				//
		// Save MCMC iterates															//
		//																				//
		//////////////////////////////////////////////////////////////////////////////////
		if((i > (*burn-1)) & ((i+1) % *thin == 0)){

			for(j = 0; j < *nsubject; j++){

				for(b = 0; b < (*nb); b++){
					beta[(ii*(*nb) + b)*(*nsubject) + j] = beta_iter[b*(*nsubject) + j];
				}



				for(t=0; t < max_nobs; t++){
					zout[(ii*(*nsubject) + j)*max_nobs + t] = z_iter[j*max_nobs + t];
					fprime[(ii*(*nsubject) + j)*max_nobs + t] = fprime_iter[j*max_nobs + t];
				}



				sig2[ii*(*nsubject) + j] = sig2_iter[j];
				u[ii*(*nsubject) + j] = u_iter[j];
				D[ii*(*nsubject) + j] = D_iter[j];
			}


			for(k = 0; k < *ng; k++){

				tau2[ii*(*ng) + k] = tau2h[k];
				lam[ii*(*ng) + k] = lamh[k];
				for(b = 0; b < (*nb); b++){
					theta[(ii*(*nb) + b)*(*ng) + k] = thetah[(k)*(*nb) + b] ;
				}

				for(t=0; t < *nrHt; t++){
					fgprime[(ii*(*ng) + k)*(*nrHt) + t] = fgprime_iter[k*(*nrHt) + t];
				}

			}


			ii = ii+1;
		}

/**/
	}

	PutRNGstate();




}


