// Perform a fit and calculate the error bands  
// using the GetConfidenceIntervals methods
// usage:
// .L errorbands2.cc+
// errorbands2()

#include "TH1F.h"
#include "TF1.h"
#include "Math/WrappedTF1.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TVirtualFitter.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TGClient.h"

using namespace ROOT::Math;

Double_t fcn(Double_t *x, Double_t *par){
  return par[0] + par[1]*TMath::Log(x[0]) 
    + par[2]*TMath::Log(x[0])*TMath::Log(x[0]);
}




void errorbands2(){
  UInt_t dh = gClient->GetDisplayHeight();

  TCanvas *c1=new TCanvas("c1","fit",50,60,dh/2,dh/3);

  TFile *tf=new TFile("mydata.root");
  TH1F *data=(TH1F*)tf->Get("data");
  const Int_t nPars=3;
  Double_t pars[nPars], grad[nPars];
  
  TF1 model("model",fcn,1,500,nPars);
  model.SetParameter(0,1.23);
  model.SetParameter(1,.5);
  model.SetParameter(2,.2);
  data->Fit(&model); // alternative to referencing the TF1 by name: data->Fit("model");
  data->SetTitle("Fit with 68% and 95% CL bands");
  data->Draw();
  model.GetParameters(pars);
  
  // Covariance matrix
  TVirtualFitter *fitter = TVirtualFitter::GetFitter();  // interface to the extract fitter info
  // Create TGraphErrors to hold the confidence intervals
  Int_t nbins = data->GetNbinsX();
  TGraphErrors *tg68 = new TGraphErrors(nbins);  // 68% CL band
  TGraphErrors *tg95 = new TGraphErrors(nbins);  // 95% CL band
  //for (int i=1; i<=nbins; i++){
  //tg68->SetPoint(i-1, data->GetBinCenter(i), 0);
  //tg95->SetPoint(i-1, data->GetBinCenter(i), 0);
  //}
  fitter->GetConfidenceIntervals(tg68,0.68);
  fitter->GetConfidenceIntervals(tg95,0.95);
  tg68->SetFillColor(kYellow);
  tg95->SetFillColor(kGreen);
  tg95->Draw("3");
  tg68->Draw("3");
  data->Draw("same");

  return;
}


