#include "mex.h"

void mexFunction (int nlhs, mxArray* plhs[],
		  int nrhs, const mxArray* prhs[]){

  if ( nrhs != 2 )
    mexErrMsgTxt("You must provide two complex arrays \n");
  
  mwSize lveca = mxGetNumberOfElements(prhs[0]);
  mwSize lvecb = mxGetNumberOfElements(prhs[1]);

  if ( lveca != lvecb )
    mexErrMsgTxt("Arrays must have same length \n");

  mwSize lvec = lveca;
  
  double *ar = mxGetPr(prhs[0]);
  double *ai = mxGetPi(prhs[0]);

  double *br = mxGetPr(prhs[1]);
  double *bi = mxGetPi(prhs[1]);
  
  plhs[0] = mxCreateDoubleMatrix(lvec, 1, mxREAL);
  plhs[1] = mxCreateDoubleMatrix(lvec, 1, mxREAL);

  double *cr = mxGetPr(plhs[0]);
  double *ci = mxGetPr(plhs[1]);

  
  for ( mwSize n = 0; n<lvec; n++)  {
    cr[n] = 0.0f;
    ci[n] = 0.0f;

    for ( mwSize nn=0; nn<lvec-n; nn++ ){
      cr[n] += ar[nn]*br[nn+n] - ai[nn]*bi[nn+n];  
      ci[n] += ar[nn]*bi[nn+n] + ai[nn]*br[nn+n];
    }
  }
  
  for ( mwSize n=0; n<lvec; n++ ){
    double fac = 1.0/(lvec-n);
    
    cr[n] *= fac;
    ci[n] *= fac;
  }

  
} 



  
/*
void mexFunction (int nlhs, mxArray* plhs[],
		  int nrhs, const mxArray* prhs[]) {

  if ( nrhs != 2 )
    mexErrMsgTxt("You must provide two complex arrays \n");    

  if ( !mxIsComplex(prhs[0]) || !mxIsComplex(prhs[1]) )
    mexErrMsgTxt("You must provide two complex arrays \n");    
    
  mxComplexDouble *a = mxGetComplexDoubles(prhs[0]);
  mxComplexDouble *b = mxGetComplexDoubles(prhs[1]);

  unsigned lvec = mxGetM(prhs[0]);
  if ( lvec == 1 ) lvec = mxGetN(prhs[0]);
  
  plhs[0] = mxCreateDoubleMatrix(lvec, 1, mxCOMPLEX);
  mxComplexDouble *ptrl= mxGetComplexDoubles(plhs[0]);
  
  for ( unsigned n=0; n<lvec; n++ ){
    ptrl[n].real = 0.0;
    ptrl[n].imag = 0.0;
 
    for ( unsigned nn=0; nn<lvec-n; nn++ ){
      ptrl[n].real += a[nn].real*b[nn+n].real - a[nn].imag*b[nn+n].imag;  
      ptrl[n].imag += a[nn].real*b[nn+n].imag + a[nn].imag*b[nn+n].real;
    }
  }
  
  for ( unsigned n=0; n<lvec; n++ ){
    double fac = 1.0/(lvec-n);

    ptrl[n].real *= fac;
    ptrl[n].imag *= fac;
  }
  
  
}
*/
