// Perform a fit and calculate the error bands  
// using the GetConfidenceIntervals methods
// usage:
// .L errorbandsCL_demo.cc(+)
// errorbandsCL_demo()

#include "TH1F.h"
#include "TF1.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TVirtualFitter.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TGClient.h"

using namespace ROOT::Math;

Double_t fcn(Double_t *x, Double_t *par){
  return par[0] + par[1]*TMath::Log(x[0]) 
    + par[2]*TMath::Log(x[0])*TMath::Log(x[0]);
}


void errorbandsCL_demo(){
  UInt_t dh = gClient->GetDisplayHeight();
  TCanvas *c1=new TCanvas("errorbandsCL_demo","fit",50,60,dh/2,dh/3);

  const Int_t nPars=3;
  TF1 model("model",fcn,1,500,nPars);
  model.SetParameter(0,1.23);
  model.SetParameter(1,.5);
  model.SetParameter(2,.2);
  TH1F *data=new TH1F("data","data;x;y",100,0,500);
  data->Sumw2();
  data->FillRandom("model",2000);

  
  data->Fit(&model); // alternative to referencing the TF1 by name: data->Fit("model");
  data->SetTitle("Fit with 68% and 95% CL bands");
  data->Draw();
  
  // Covariance matrix
  TVirtualFitter *fitter = TVirtualFitter::GetFitter();  // interface to the extract fitter info
  // Create TGraphErrors to hold the confidence intervals
  Int_t nbins = data->GetNbinsX();
  TGraphErrors *tg68 = new TGraphErrors(nbins);  // 68% CL band
  TGraphErrors *tg95 = new TGraphErrors(nbins);  // 95% CL band
  for (int i=1; i<=nbins; i++){
    tg68->SetPoint(i-1, data->GetBinCenter(i), 0);
    tg95->SetPoint(i-1, data->GetBinCenter(i), 0);
  }
  // error propogation here
  fitter->GetConfidenceIntervals(tg68,0.68);
  fitter->GetConfidenceIntervals(tg95,0.95);
  tg68->SetFillColor(kYellow);
  tg95->SetFillColor(kGreen);
  tg95->Draw("3");
  tg68->Draw("3");
  data->Draw("same");

  return;
}


