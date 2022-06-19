#include "Math/WrappedTF1.h"
using namespace ROOT::Math;

// here we can do error propagation for an arbitrary function
// matrix and vector notation is used
// return error on function sigma_f(x;par,COV) evaluated at x
double GetCovError(TF1 &f, Double_t x, Double_t *pars, TMatrixD &COV){
  auto wf=WrappedTF1(f);
  Int_t npar=COV.GetNrows();
  Double_t *grad=new Double_t[npar];
  wf.ParameterGradient(x,pars,grad);
  TVectorD g(npar,grad);
  double var = g*(COV*g);
  delete[] grad;
  return TMath::Sqrt(var);
}


TMatrixD MtimesM(TMatrixD &M1,TMatrixD &M2){
    return M1*M2;
}

TVectorD MtimesV(TMatrixD &M,TVectorD &v){
    return M*v;
}


