/* 
   Name: mpmm.c
   Author: Christopher DuBois
   Date: 2/5/2010
*/
#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>

int draw (double probs[], int K) // returns index in [0,K-1]
{
  int k;
  double cuml = 0;
  double sum = 0;
  for (k=0;k<K;k++) {
    sum += probs[k];
  }
  double r = unif_rand();
  for (k=0;k<K;k++) {
    cuml += probs[k]/sum;
    if (r < cuml) {
      break;
    } 
  }
  return k;
}

SEXP mpmm_cgs (SEXP s, SEXP r, SEXP a,
	       SEXP mode_counts,
               SEXP sender_counts,
               SEXP receiver_counts,
	       SEXP action_counts,
               SEXP alpha, SEXP beta, SEXP gamma, SEXP delta,
	       SEXP num_iterations, SEXP assignments)
{
  GetRNGstate();
  int nProtected = 0;
  int T, K, S, R, A, i,j,k,m,iter,niter;
  niter = INTEGER(num_iterations)[0];
  SEXP retval, n0, n1, n2, n3, zh, prb;
  PROTECT(prb = allocVector(REALSXP, K));  ++nProtected;
  T = length(assignments);
  K = length(mode_counts);
  S = ncols(sender_counts);
  R = ncols(receiver_counts);
  A = ncols(action_counts);
  // TODO: Check dimensions of other count matrices.
  // Rprintf("K:%i, S:%i, R:%i, A:%i\n", K,S,R,A);

  // Copy variables  to new matrices.
  PROTECT(zh = allocVector(INTSXP, T));  ++nProtected;
  memcpy(INTEGER(zh), INTEGER(assignments), sizeof(int) * T);

  PROTECT(n0 = allocVector(INTSXP, K));  ++nProtected;
  memcpy(INTEGER(n0), INTEGER(mode_counts), sizeof(int) * K);

  PROTECT(n1 = allocMatrix(INTSXP, K, S));  ++nProtected;
  memcpy(INTEGER(n1), INTEGER(sender_counts), sizeof(int)*K*S);

  PROTECT(n2 = allocMatrix(INTSXP, K, R));  ++nProtected;
  memcpy(INTEGER(n2), INTEGER(receiver_counts), sizeof(int)*K*R);

  PROTECT(n3 = allocMatrix(INTSXP, K, A));  ++nProtected;
  memcpy(INTEGER(n3), INTEGER(action_counts), sizeof(int)*K*A);

  /* Get row sums of count matrices */
  int sender_rowsums[K];
  int receiver_rowsums[K];
  int action_rowsums[K];
  for (k=0;k<K;++k){
    sender_rowsums[k] = 0;
    for (m=0;m<S;++m){
      sender_rowsums[k] += INTEGER(sender_counts)[k + K*m];
    }
    receiver_rowsums[k] = 0;
    for (m=0;m<R;++m){
      receiver_rowsums[k] += INTEGER(receiver_counts)[k + K*m];
    }
    action_rowsums[k] = 0;
    for (m=0;m<A;++m){
      action_rowsums[k] += INTEGER(action_counts)[k + K*m];
    }
  }

  // Get constants ready.
  // TODO: Allow for vectors of the Dirichlet parameters
  double probs[K];
  double ralpha = REAL(alpha)[0];
  double rbeta = REAL(beta)[0];
  double rgamma = REAL(gamma)[0];
  double rdelta = REAL(delta)[0];

  for (i=0;i<T;++i) {
    INTEGER(zh)[i] -= 1;
    INTEGER(s)[i] -= 1;
    INTEGER(r)[i] -= 1;
    INTEGER(a)[i] -= 1;
  }

  for (iter=0;iter<niter;++iter){ // Number of Gibbs scans
    /* Rprintf("Iteration: %i\n",iter); */
    if (iter%10 == 0) {
      Rprintf(".");
    }
    for (i=0;i<T;i++){ // Each event // switch to T
      
      // Decrement counts
      INTEGER(n0)[INTEGER(zh)[i]] -= 1;
      INTEGER(n1)[INTEGER(zh)[i] + K * INTEGER(s)[i]] -= 1;
      INTEGER(n2)[INTEGER(zh)[i] + K * INTEGER(r)[i]] -= 1;
      INTEGER(n3)[INTEGER(zh)[i] + K * INTEGER(a)[i]] -= 1;
      sender_rowsums[INTEGER(zh)[i]] -= 1;
      receiver_rowsums[INTEGER(zh)[i]] -= 1;
      action_rowsums[INTEGER(zh)[i]] -= 1;

      // Probability vector for new assignment
      for (k=0;k<K;++k) {  
  	probs[k] =
	(INTEGER(n0)[k] + ralpha) *
	(INTEGER(n1)[k + K * INTEGER(s)[i]] + rbeta) *
	(INTEGER(n2)[k + K * INTEGER(r)[i]] + rgamma) *
	(INTEGER(n3)[k + K * INTEGER(a)[i]] + rdelta) /
	(sender_rowsums[k] + S * rbeta) /
	(receiver_rowsums[k] + R * rgamma) /
  	(action_rowsums[k] + A * rdelta);
      }
      
      INTEGER(zh)[i] = draw(probs,K);

      // Increment counts
      INTEGER(n0)[INTEGER(zh)[i]] += 1;
      INTEGER(n1)[INTEGER(zh)[i] + K * INTEGER(s)[i]] += 1;
      INTEGER(n2)[INTEGER(zh)[i] + K * INTEGER(r)[i]] += 1;
      INTEGER(n3)[INTEGER(zh)[i] + K * INTEGER(a)[i]] += 1;
      sender_rowsums[INTEGER(zh)[i]] += 1;
      receiver_rowsums[INTEGER(zh)[i]] += 1;
      action_rowsums[INTEGER(zh)[i]] += 1;
    }
  }
  for (i=0;i<T;++i) {
    INTEGER(zh)[i] += 1;
  }

  PROTECT(retval = allocVector(VECSXP, 5)); ++nProtected;
  SET_VECTOR_ELT(retval, 0, n0);
  SET_VECTOR_ELT(retval, 1, n1);
  SET_VECTOR_ELT(retval, 2, n2);
  SET_VECTOR_ELT(retval, 3, n3);
  SET_VECTOR_ELT(retval, 4, zh);

  SEXP names;
  PROTECT(names = allocVector(STRSXP, 5));
  ++nProtected;
  SET_STRING_ELT(names, 0, mkChar("mode_counts"));
  SET_STRING_ELT(names, 1, mkChar("sender_counts"));
  SET_STRING_ELT(names, 2, mkChar("receiver_counts"));
  SET_STRING_ELT(names, 3, mkChar("action_counts"));
  SET_STRING_ELT(names, 4, mkChar("assignments"));
  setAttrib(retval, R_NamesSymbol, names);

  UNPROTECT(nProtected);
  PutRNGstate();
  return retval; 
}

