// Perform a fit and calculate the error bands based on 
// the error matrix
// usage:
// .L errorbands_demo.cc+
// errorbands_demo()

#include "TH1F.h"
#include "TF1.h"
#include "Math/WrappedTF1.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TVirtualFitter.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TGClient.h"

using namespace ROOT::Math;


Double_t lnfcn(Double_t *x, Double_t *par){
  return par[0] + par[1]*TMath::Log(x[0]) 
    + par[2]*TMath::Log(x[0])*TMath::Log(x[0]);
}

// not used, just for show
// calculate 1 sigma uncertainy of fcn model at 'x' wrt parameter uncertainties
// inputs
// grad: array of parameter gradents df/dpar evaluated at 'x'
// COV: covariance matrix for parameters
// npar: # of params
double sigmaFit(Double_t *grad, TMatrixD &COV, Int_t npar){
  double sigma=0;

  // s^2 = sum_i sum_j df/dp_i df/dp_j COV_ij
  for (int pi=0; pi<npar; pi++){
    for (int pj=0; pj<npar; pj++){
      sigma+=grad[pi]*grad[pj]*COV[pi][pj];
    }
  }
  //return sigma;
  return TMath::Sqrt(sigma);
}



// here we can do error propagation for an arbitrary function
// matrix and vector notation is used
// return error on function evaluated at x, sigma_f(x;par)
double GetError(WrappedTF1 &f, Double_t x, Double_t *pars, TMatrixD &COV, Int_t npar){
  Double_t *grad=new Double_t[npar];
  f.ParameterGradient(x,pars,grad);
  TVectorD g(npar,grad);
  double var = g*(COV*g);
  delete[] grad;
  return TMath::Sqrt(var);
}


void errorbands_demo(){
  UInt_t dh = gClient->GetDisplayHeight();
  TCanvas *c1=new TCanvas("c1","fit",50,60,dh/2,dh/3);

  const Int_t nPars=3;
  Double_t pars[nPars], grad[nPars];  
  TF1 model("model",lnfcn,1,500,nPars);
  model.SetParameter(0,1.23);
  model.SetParameter(1,.5);
  model.SetParameter(2,.2);
  TH1F *data=new TH1F("data","data;x;y",100,0,500);
  data->Sumw2();
  data->FillRandom("model",2000);
  data->Fit(&model); // alternative to referencing the TF1 by name: data->Fit("model");
  data->Draw();
  model.GetParameters(pars); // retrieve post fit parameters
  
  WrappedTF1 wmodel(model); // provides additional methods to our TF1

  // Get the Covariance matrix
  TVirtualFitter *fitter = TVirtualFitter::GetFitter();  // interface to the extract fitter info
  assert (nPars == fitter->GetNumberFreeParameters());
  TMatrixD* COV = new TMatrixD( nPars, nPars, fitter->GetCovarianceMatrix() );  
  TH1F *up=new TH1F(*data);  up->SetLineColor(kBlue);  up->SetLineWidth(2);
  TH1F *dn=new TH1F(*data);  dn->SetLineColor(kBlue);  dn->SetLineWidth(2);
  // calculate error bands by hand
  for (int ib=1; ib<=data->GetNbinsX(); ib++){
    double x=data->GetBinCenter(ib);
    double sigma=GetError(wmodel,x,pars,*COV,nPars);
    up->SetBinContent(ib,model(x)+sigma);  // multipy sigma by 2.0 to get 2 sigma, etc
    dn->SetBinContent(ib,model(x)-sigma);
  }
  up->Draw("same,hist,c");
  dn->Draw("same,hist,c");
  c1->Print("fit.png");

  TCanvas *c2=new TCanvas("c2","bands",dh/2+80,60,dh/2,dh/3);

  TH1F *uprat=new TH1F(*up); uprat->Divide(&model);
  TH1F *dnrat=new TH1F(*dn); dnrat->Divide(&model);
  uprat->SetLineWidth(2);  
  dnrat->SetLineWidth(2);  

  // data/fit
  TH1F *data2=(TH1F*)data->Clone("ratio");
  data2->SetTitle("Data to fit ratio;x;data/fit");
  data2->Divide(&model);
  data2->SetMaximum(data2->GetMaximum()*1.1);
  data2->Draw();
  
  uprat->Draw("hist,c,same");
  dnrat->Draw("hist,c,same");

  c2->Print("bands.png");
}


