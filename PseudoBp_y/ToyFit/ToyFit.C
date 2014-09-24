#include <TH1D.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TNtuple.h>
#include <TStyle.h>
void ToyFit(int bin = 1){
  gStyle->SetOptStat("mr");
  gStyle->SetStatFontSize(0.055);

  int binlow = 0;
  int binhigh = 0;
  if(bin==0){
    binlow = 0; binhigh = 300;
  }
  if(bin==1){
    binlow = 0; binhigh = 300;
  }
  if(bin==2){
    binlow = 0; binhigh = 300;
  }
  if(bin==3){
    binlow = 0; binhigh = 300;
  }

  TFile *inf = new TFile("../ToMerge/SigmaBplus.root");//
  TTree *nt = (TTree*) inf->Get("pseudo");
  TH1D *toys;
  TH1D *true_y;
  toys = new TH1D("toys","yields",20,binlow,binhigh);   
  true_y = new TH1D("true_y","yields",20,binlow,binhigh);  
//	toys = new TH1D("toys","yields",20,150,500);
//	true_y = new TH1D("true_y","yields",20,150,500);                                   
  TCanvas *c=  new TCanvas("c","",600,600);                              
  c->cd();
  nt->Project("toys","yield",Form("bin==%d",bin));
  nt->Project("true_y","true_yield",Form("bin==%d",bin));
  TF1 *f = new TF1("f","[0]*Gaus(x,[1],[2])/(sqrt(2*3.14159)*[2])");
  f->SetParameter(0,1000);
  f->SetParameter(1,150);
  f->SetParameter(2,50);
  f->SetMarkerSize(5);
  toys->Fit("f","L m","",binlow,binhigh);
  toys->Draw("");
  toys->Draw("pe");
  f->SetLineColor(2);
  f->SetLineWidth(5);
  f->Draw("same");

  TLegend *leg = new TLegend(0.08,0.97,0.86,0.99);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.035);
  if(bin==0)  leg->AddEntry(toys,Form("y_{CM} -1.93~-1, true#: %.2f", true_y->GetMean()),"");
  if(bin==1)  leg->AddEntry(toys,Form("y_{CM} -1~0, true#: %.2f", true_y->GetMean()),"");
  if(bin==2)  leg->AddEntry(toys,Form("y_{CM} 0~1, true#: %.2f", true_y->GetMean()),"");
  if(bin==3)  leg->AddEntry(toys,Form("y_{CM} 1~1.93, true#: %.2f", true_y->GetMean()),"");
  leg->Draw();
  c->SaveAs(Form("bin%d.pdf",bin));
  c->SaveAs(Form("bin%d.gif",bin));

  TCanvas *c2=  new TCanvas("c2","",600,600);
  c2->cd();
  float nobs = 0;
  float yield;
  float yield_err;
  float true_yield;
  float b;
  float data_obs;
  float obs_err;
  nt->SetBranchAddress("bin", &b);
  nt->SetBranchAddress("yield", &yield);
  nt->SetBranchAddress("yield_err", &yield_err);
  nt->SetBranchAddress("true_yield", &true_yield);
  nt->SetBranchAddress("data_obs", &data_obs);
  nt->SetBranchAddress("obs_err", &obs_err);
  TH1D* pull = new TH1D("pull","pull",40,-10,10);
  TH1D* diff = new TH1D("diff","diff",50,-100,100);
  int nevt = nt->GetEntries();
  for(int ev=0; ev<nevt; ev++){
    nt->GetEntry(ev);
    if(b==bin) {
        if(nobs==0) nobs=data_obs;
      pull->Fill((yield-data_obs)/yield_err);
      diff->Fill((yield-true_yield));
    }
  }
  TF1 *fpull = new TF1("fpull","[0]*Gaus(x,[1],[2])/(sqrt(2*3.14159)*[2])");             
  fpull->SetParameter(0,1000);
  fpull->SetParameter(1,0);
  fpull->SetParameter(2,1);
  pull->Fit("fpull","L m","",-10,10);
  fpull->SetMarkerSize(5);
  fpull->SetLineColor(2);
//  fpull->SetLineWidth(5);
  pull->SetStats(1);
  pull->Draw("pe");
  fpull->Draw("same");
  leg->Draw("same");
//  diff->Draw("pe");
  c2->SaveAs(Form("pullbin%d.pdf",bin));
  c2->SaveAs(Form("pullbin%d.gif",bin));

  cout<<"true yield mean: :"<<true_y->GetMean()<<endl;
  cout<<"obs: :"<<nobs<<endl;
}
void doall()
{
ToyFit(0);
ToyFit(1);
ToyFit(2);
ToyFit(3);
}
