#include "utilities.h"
#define MAX_XB 8192
#define MAX_GEN 4096
using namespace std;

// current date/time based on current system
time_t now = time(0);   
char* dt = ctime(&now);

double luminosity=34.8*1e-3;

double setparam0=5.85733e+01;
double setparam1=5.28008e+00;
double setparam2=5.00000e-02;
double setparam3=1.02740e+02;
double setparam4=-1.62289e+01;
double setparam5=0.00000e+00;
double setparam6=1.43809e-01;
double setparam7=1.86128e-01;
double setparam8=1.42964e-02;

double fixparam1=5.279;
double fixparam2=0.04;

TString inputdata="/mnt/hadoop/cms/store/user/jwang/nt_20140403_PAMuon_HIRun2013_PromptrecoAndRereco_v1_MuonMatching_EvtBase_skim.root";
//TString inputdata="/net/hidsk0001/d00/scratch/jwang/nt_201403019_PAMuon_HIRun2013_PromptrecoAndRereco_v1_MuonMatching_EvtBase.root";
//TString inputdata="/mnt/hadoop/cms/store/user/jwang/nt_201403019_PAMuon_HIRun2013_PromptrecoAndRereco_v1_MuonMatching_EvtBase.root";
//TString inputmc="/mnt/hadoop/cms/store/user/jwang/nt_BoostedMC_20140318_Kstar_TriggerMatchingMuon_EvtBase.root";
//TString inputmc="/net/hidsk0001/d00/scratch/jwang/nt_BoostedMC_20140318_Kstar_TriggerMatchingMuon_EvtBase.root";
TString inputmc="/mnt/hadoop/cms/store/user/jwang/nt_BoostedMC_20140403_Kstar_TriggerMatchingMuon_EvtBase_skim.root";

TString cut_kpi="(HLT_PAMu3_v1)&&abs(mumumass-3.096916)<0.15&&chi2cl>0.15&&(d0)/d0Err>8.1&&cos(dtheta)>-0.44&&TMath::Abs((trk2Dxy)/trk2D0Err)>0.81&&abs(tktkmass-0.89591)<0.14&&mass>5&&mass<6";
TString cut_pik="(HLT_PAMu3_v1)&&abs(mumumass-3.096916)<0.15&&chi2cl>0.15&&(d0)/d0Err>8.1&&cos(dtheta)>-0.44&&TMath::Abs((trk1Dxy)/trk1D0Err)>0.81&&abs(tktkmass-0.89591)<0.14&&mass>5&&mass<6";
//TString cut="chi2cl>0.15&&(d0)/d0Err>8.1&&cos(dtheta)>-0.44&&TMath::Abs((trk1Dxy)/trk1D0Err)>0.81&&abs(tktkmass-0.89591)<0.14&&mass>5&&mass<6";
//abs(y+0.465)<1.93
TString seldata_kpi=Form("abs(y+0.465)<1.93&&%s",cut_kpi.Data());
TString seldata_pik=Form("abs(y+0.465)<1.93&&%s",cut_pik.Data());
TString seldata_2y_kpi=Form("((Run>=210498&&Run<=211256&&abs(y+0.465)<1.93)||(Run>=211313&&Run<=211631&&abs(y-0.465)<1.93))&&%s",cut_kpi.Data());
TString seldata_2y_pik=Form("((Run>=210498&&Run<=211256&&abs(y+0.465)<1.93)||(Run>=211313&&Run<=211631&&abs(y-0.465)<1.93))&&%s",cut_pik.Data());

TString selmc_kpi=Form("abs(y+0.465)<1.93&&(gen==22233||gen==41000)&&%s",cut_kpi.Data());
TString selmc_pik=Form("abs(y+0.465)<1.93&&(gen==22233||gen==41000)&&%s",cut_pik.Data());

TString selmcgen="abs(y+0.465)<1.93&&abs(pdgId)==511&&isSignal!=0";

TString weight = "27.493+pt*(-0.218769)";

void clean0(TH1D *h)
{
   for (int i=1;i<=h->GetNbinsX();i++)
   {
      if (h->GetBinContent(i)==0) h->SetBinError(i,1);
   }
}

//TF1 *fit(TTree *nt,TTree *nt2, TTree *ntMC, TTree *ntMC2,double ptmin,double ptmax){   
//std::vector<TF1*> fit(TTree *nt,TTree *nt2, TTree *ntMC, TTree *ntMC2, double ptmin,double ptmax, bool dopseudo, TH1D *pseudo){
TF1* fit(TTree *nt,TTree *nt2, TTree *ntMC, TTree *ntMC2, double ptmin,double ptmax, bool dopseudo, TH1D *pseudo, TF1 *pseudotf1){
   //cout<<cut.Data()<<endl;
   static TF1 return_th1[3];
   static int count=0;
   count++;
   TCanvas *c= new TCanvas(Form("c%d",count),"",600,600);
   TH1D *h = new TH1D(Form("h%d",count),"",50,5,6);
   TH1D *h2 = new TH1D(Form("h2%d",count),"",50,5,6);
//   TH1D *hBck = new TH1D(Form("hBck%d",count),"",50,5,6);
   
   TH1D *hMC = new TH1D(Form("hMC%d",count),"",50,5,6);
   TH1D *hMC2 = new TH1D(Form("hMC2%d",count),"",50,5,6);
   // Fit function
   TF1 *f = new TF1(Form("f%d",count),"[0]*([7]*Gaus(x,[1],[2])/(sqrt(2*3.14159)*[2])+(1-[7])*Gaus(x,[1],[8])/(sqrt(2*3.14159)*[8]))+[3]+[4]*x+[6]*(38.42*Gaus(x,5.25,0.03473)+15.04*Gaus(x,5.25,0.1121)+104.3*Gaus(x,5.026,0.0935))");
   if(!dopseudo){
     nt->Project(Form("h%d",count),"mass",Form("%s&&pt>%f&&pt<%f",seldata_2y_kpi.Data(),ptmin,ptmax));   
     ntMC->Project(Form("hMC%d",count),"mass",Form("%s&&pt>%f&&pt<%f",seldata_kpi.Data(),ptmin,ptmax));   
     nt2->Project(Form("h2%d",count),"mass",Form("%s&&pt>%f&&pt<%f",seldata_2y_pik.Data(),ptmin,ptmax));   
     ntMC2->Project(Form("hMC2%d",count),"mass",Form("%s&&pt>%f&&pt<%f",seldata_pik.Data(),ptmin,ptmax));   
     h->Add(h2);
     hMC->Add(hMC2);
    }
	if(dopseudo) hMC = (TH1D*)pseudo->Clone(Form("h%d",count));
	if(dopseudo) h = (TH1D*)pseudo->Clone(Form("h%d",count));

//   nt->Project(Form("hBck%d",count),"mass",Form("%s&&pt>%f&&pt<%f&&(gen==22233||gen==41000)",seldata.Data(),ptmin,ptmax));   
//   nt2->Project(Form("hBck%d",count),"mass",Form("%s&&pt>%f&&pt<%f&&(gen==22233||gen==41000)",seldata.Data(),ptmin,ptmax));   

//   cout <<"nsig = "<<hBck->GetEntries();
//   clean0(h);
   h->Draw();
   f->SetParLimits(4,-1000,0);
   f->SetParLimits(2,0.01,0.05);
   f->SetParLimits(8,0.01,0.1);
   f->SetParLimits(7,0,1);

   f->SetParameter(0,setparam0);
   f->SetParameter(1,setparam1);
   f->SetParameter(2,setparam2);
   f->SetParameter(3,setparam3);
   f->SetParameter(4,setparam4);
   f->SetParameter(5,setparam5);
   f->SetParameter(6,setparam6);
   f->SetParameter(7,setparam7);
   f->SetParameter(8,setparam8);

if(dopseudo){
   f->SetParameter(0,pseudotf1->GetParameter(0));
   f->SetParameter(1,pseudotf1->GetParameter(1));
   f->SetParameter(2,pseudotf1->GetParameter(2));
   f->SetParameter(3,pseudotf1->GetParameter(3));
   f->SetParameter(4,pseudotf1->GetParameter(4));
   f->SetParameter(5,pseudotf1->GetParameter(5));
   f->SetParameter(6,pseudotf1->GetParameter(6));
   f->SetParameter(7,pseudotf1->GetParameter(7));
   f->SetParameter(8,pseudotf1->GetParameter(8));
}

   f->FixParameter(1,fixparam1);
   h->GetEntries();
   hMC->GetEntries();

if(!dopseudo){
//   hMC->Fit(Form("f%d",count),"q","",5,6);
   hMC->Fit(Form("f%d",count),"q","",5,6);
   f->ReleaseParameter(1);
   hMC->Fit(Form("f%d",count),"","",5,6);
//   hMC->Fit(Form("f%d",count),"L q","",5,6);
//   hMC->Fit(Form("f%d",count),"L q","",5,6);
//   hMC->Fit(Form("f%d",count),"L q","",5,6);
//   hMC->Fit(Form("f%d",count),"L m","",5,6);
}

   f->FixParameter(1,f->GetParameter(1));
//   f->FixParameter(2,f->GetParameter(2));
//   f->FixParameter(7,f->GetParameter(7));
//   f->FixParameter(8,f->GetParameter(8));
   
//   h->Fit(Form("f%d",count),"q","",5,6);
   h->Fit(Form("f%d",count),"q","",5,6);
   f->ReleaseParameter(1);
   h->Fit(Form("f%d",count),"","",5,6);
//   h->Fit(Form("f%d",count),"L q","",5,6);
//   h->Fit(Form("f%d",count),"L q","",5,6);
//   h->Fit(Form("f%d",count),"L q","",5,6);
//   h->Fit(Form("f%d",count),"L m","",5,6);

   h->SetMarkerSize(0.8);
   h->SetMarkerStyle(20);
   cout <<h->GetEntries()<<endl;
   cout <<hMC->GetEntries()<<endl;

   // function for background shape plotting. take the fit result from f
   TF1 *background = new TF1(Form("background%d",count),"[0]+[1]*x+[3]*(38.42*Gaus(x,5.25,0.03473)+15.04*Gaus(x,5.25,0.1121)+104.3*Gaus(x,5.026,0.0935))");
//   TF1 *background = new TF1(Form("background%d",count),"[0]+[1]*x+[2]*(1.24e2*Gaus(x,5.107,0.02987)+1.886e2*Gaus(x,5.0116,5.546e-2))");
   background->SetParameter(0,f->GetParameter(3));
   background->SetParameter(1,f->GetParameter(4));
   background->SetParameter(2,f->GetParameter(5));
   background->SetParameter(3,f->GetParameter(6));
   background->SetLineColor(4);
   background->SetRange(5,6);
   background->SetLineStyle(2);
   
   // function for signal shape plotting. take the fit result from f
   TF1 *Bkpi = new TF1(Form("fBkpi",count),"[0]*(38.42*Gaus(x,5.25,0.03473)+15.04*Gaus(x,5.25,0.1121)+104.3*Gaus(x,5.026,0.0935))");
   Bkpi->SetParameter(0,f->GetParameter(6));
   Bkpi->SetLineColor(kGreen+1);
   Bkpi->SetFillColor(kGreen+1);
   Bkpi->SetRange(5.00,5.45);
   Bkpi->SetLineStyle(1);
   Bkpi->SetFillStyle(3005);

   // function for signal shape plotting. take the fit result from f
   TF1 *mass = new TF1(Form("fmass",count),"[0]*([3]*Gaus(x,[1],[2])/(sqrt(2*3.14159)*[2])+(1-[3])*Gaus(x,[1],[4])/(sqrt(2*3.14159)*[4]))");
   mass->SetParameters(f->GetParameter(0),f->GetParameter(1),f->GetParameter(2),f->GetParameter(7),f->GetParameter(8));
   mass->SetParError(0,f->GetParError(0));
   mass->SetParError(1,f->GetParError(1));
   mass->SetParError(2,f->GetParError(2));
   mass->SetParError(7,f->GetParError(7));
   mass->SetParError(8,f->GetParError(8));
   mass->SetLineColor(2);
   mass->SetLineStyle(2);

//   cout <<mass->Integral(0,1.2)<<" "<<mass->IntegralError(0,1.2)<<endl;
   h->SetMarkerStyle(24);
   h->SetStats(0);
   h->Draw("e");
   h->SetXTitle("M_{B} (GeV/c^{2})");
   h->SetYTitle("Entries / (20 MeV/c^{2})");
   h->GetXaxis()->CenterTitle();
   h->GetYaxis()->CenterTitle();
   h->SetTitleOffset(1.4,"Y");
   h->SetAxisRange(0,h->GetMaximum()*1.2,"Y");

 //  hBck->Draw("hist same");

   Bkpi->Draw("same");
   background->Draw("same");   
   mass->SetRange(5,6);
   mass->Draw("same");
   mass->SetLineStyle(2);
   mass->SetFillStyle(3004);
   mass->SetFillColor(2);
   f->Draw("same");

   double yield = mass->Integral(5,6)/0.02;
   double yieldErr = mass->Integral(5,6)/0.02*mass->GetParError(0)/mass->GetParameter(0);
   cout <<"fit result:"<<mass->GetParameter(0)/h->GetBinWidth(1)<<" "<<mass->Integral(5,6)/h->GetBinWidth(1)<<endl;
   
   // Draw the legend:)   
   TLegend *leg = myLegend(0.50,0.5,0.86,0.92);
   leg->AddEntry(h,"CMS Preliminary","");
   leg->AddEntry(h,"p+Pb #sqrt{s_{NN}}= 5.02 TeV","");
   leg->AddEntry(h,Form("%.0f<p_{T}^{B}<%.0f GeV/c",ptmin,ptmax),"");
   leg->AddEntry(h,"Data","pl");
   leg->AddEntry(f,"Fit","l");
   leg->AddEntry(mass,"Signal","f");
   leg->AddEntry(background,"Combinatorial Background","l");
   leg->AddEntry(Bkpi,"Non-prompt J/#psi","f");
   leg->Draw();
   TLegend *leg2 = myLegend(0.44,0.33,0.89,0.50);
   leg2->AddEntry(h,"B meson","");
   leg2->AddEntry(h,Form("M_{B}=%.2f #pm %.2f MeV/c^{2}",f->GetParameter(1)*1000.,f->GetParError(1)*1000.),"");
   leg2->AddEntry(h,Form("N_{B}=%.0f #pm %.0f",yield,yieldErr),"");
   leg2->Draw();

   c->SaveAs(Form("ResultsBzero/BMass-%d.C",count));
   c->SaveAs(Form("ResultsBzero/BMass-%d.gif",count));
//   c->SaveAs(Form("ResultsBzero/BMass-%d.eps",count));

//   return mass;
//   std::vector<TF1*> return_th1;
//   return_th1.push_back(mass);
//   return_th1.push_back(background);
//   TF1 return_th1[2] = {*mass, *background};
   return_th1[0] = *mass;
   return_th1[1] = *background;
   return_th1[2] = *f;
   return return_th1;
}

//void fitB0(TString infname="",bool doweight = 1)
void fitB0(int iseed = 1)
{

//return;
  TString infname="";
//  if (doweight==0) weight="1";
  if (infname=="") infname=inputdata.Data();
  TFile *inf = new TFile(infname.Data());
  TTree *nt = (TTree*) inf->Get("ntKstar1");
  TTree *nt2 = (TTree*) inf->Get("ntKstar2");

  TFile *infMC = new TFile(inputmc.Data());
  TTree *ntGen = (TTree*)infMC->Get("ntGen");
  TTree *ntMC = (TTree*)infMC->Get("ntKstar1");
  TTree *ntMC2 = (TTree*)infMC->Get("ntKstar2");

/*
  string active_br[50];
  active_br[0] = "HLT_PAMu3_v1";
  active_br[1] = "mumumass";
  active_br[2] = "chi2cl";
  active_br[3] = "d0";
  active_br[4] = "d0Err";
  active_br[5] = "dtheta";
  active_br[6] = "trk1Dxy";
  active_br[7] = "trk1D0Err";
  active_br[8] = "trk2Dxy";
  active_br[9] = "trk2D0Err";
  active_br[10] = "tktkmass";
  active_br[11] = "mass";
  active_br[12] = "y";
  active_br[13] = "Run";
  active_br[14] = "gen";
  active_br[15] = "pt";
  nt->SetBranchStatus("*", 0);
  nt2->SetBranchStatus("*", 0);
  ntMC->SetBranchStatus("*", 0);
  ntMC2->SetBranchStatus("*", 0);
  for(int s=0;s<16;s++){
    nt->SetBranchStatus(active_br[s].c_str(), 1);
    nt2->SetBranchStatus(active_br[s].c_str(), 1);
    ntMC->SetBranchStatus(active_br[s].c_str(), 1);
    ntMC2->SetBranchStatus(active_br[s].c_str(), 1);
  }
*/

  //nt->SetAlias("LD",LDalias.Data());
  //nt2->SetAlias("LD",LDalias.Data());
  //ntMC->SetAlias("LD",LDalias.Data());
  //ntMC2->SetAlias("LD",LDalias.Data());
  
//  const int nBins = 3;
//  double ptBins[nBins+1] = {10,15,20,60};
  const int nBins = 1;
  double ptBins[nBins+1] = {10,60};
  TH1D *hPt = new TH1D("hPt","",nBins,ptBins);
  TH1D *hRecoTruth = new TH1D("hRecoTruth","",nBins,ptBins);
  TH1D *hRecoTruth2 = new TH1D("hRecoTruth2","",nBins,ptBins);
  TH1D *hPtMC = new TH1D("hPtMC","",nBins,ptBins);
  TH1D *hPtMC2 = new TH1D("hPtMC2","",nBins,ptBins);
  TH1D *hPtGen = new TH1D("hPtGen","",nBins,ptBins);

  TH1D* dummy = new TH1D("dummy","",50, 5, 6);
  TF1* dummytf1;;
//  std::vector<TF1*> bg;
  TF1 bg[3];
  TF1 sigAndbg[3];
//  TF1* f;
  TF1 f;
//  std::vector<TF1*> results_tf1;
  TF1* results_tf1;
//  TF1 results_tf1[2];
  for (int i=0;i<nBins;i++)
    {
//      TF1 *f = fit(nt,nt2,ntMC,ntMC2,ptBins[i],ptBins[i+1]);
	  now = time(0);dt = ctime(&now);cout << "The local date and time is: " << dt << endl;
      results_tf1 = fit(nt,nt2,ntMC,ntMC2,ptBins[i],ptBins[i+1], 0, dummy, dummytf1);
	  now = time(0);dt = ctime(&now);cout << "The local date and time is: " << dt << endl;
//      f = results_tf1[0];
//      f = (TF1*)results_tf1[0].Clone();
      f = results_tf1[0];
//      bg.push_back(results_tf1[1]);
      bg[i] = results_tf1[1];
      sigAndbg[i] = results_tf1[2];
//      double yield = f->Integral(5,6)/0.02;
      double yield = f.Integral(5,6)/0.02;
      double yieldErr = f.Integral(5,6)/0.02*f.GetParError(0)/f.GetParameter(0);
      hPt->SetBinContent(i+1,yield/(ptBins[i+1]-ptBins[i]));
      hPt->SetBinError(i+1,yieldErr/(ptBins[i+1]-ptBins[i]));
  }  
  
  TCanvas *c=  new TCanvas("cResult","",600,600);
  hPt->SetXTitle("B^{0} p_{T} (GeV/c)");
  hPt->SetYTitle("Uncorrected B^{0} dN/dp_{T}");
  hPt->Sumw2();
  hPt->Draw();

  ntMC->Project("hPtMC","pt",TCut(weight)*(TCut(seldata_kpi.Data())&&"(gen==22233||gen==41000)"));
  ntMC2->Project("hPtMC2","pt",TCut(weight)*(TCut(seldata_pik.Data())&&"(gen==22233||gen==41000)"));
  hPtMC->Add(hPtMC2);

  nt->Project("hRecoTruth","pt",TCut(seldata_kpi.Data())&&"(gen==22233||gen==41000)");
  nt2->Project("hRecoTruth2","pt",TCut(seldata_pik.Data())&&"(gen==22233||gen==41000)");
  hRecoTruth->Add(hRecoTruth2);
  ntGen->Project("hPtGen","pt",TCut(weight)*(selmcgen.Data()));
  divideBinWidth(hRecoTruth);
  
  hRecoTruth->Draw("same hist");
  divideBinWidth(hPtMC);
  divideBinWidth(hPtGen);
  
  hPtMC->Sumw2();
  TH1D *hEff = (TH1D*)hPtMC->Clone("hEff");
  hPtMC->Sumw2();
  hEff->Divide(hPtGen);
  
  TH1D *hPtCor = (TH1D*)hPt->Clone("hPtCor");
  hPtCor->Divide(hEff);
  TCanvas *cCor=  new TCanvas("cCorResult","",600,600);
  hPtCor->SetYTitle("Correctd B^{0} dN/dp_{T}");
  hPtCor->Draw();
  hPtGen->Draw("same hist");

  TH1D *hPtSigma= (TH1D*)hPtCor->Clone("hPtSigma");
  hPtSigma->Scale(1./(2*luminosity));
  hPtSigma->SetYTitle("d#sigma/dp_{T} (B^{0}) ");

  TCanvas *cSigma=  new TCanvas("cSigma","",600,600);

  hPtSigma->Draw();

  TFile *outf = new TFile("ResultsBzero/SigmaBzero.root","recreate");
  outf->cd();

//delete nt;
//delete nt2;
TTree *dummynt;
  
//Pseudo
TNtuple* pseudo = new TNtuple("pseudo","","yield:bin:true_yield:yield_err");                   
int ngen = ntGen->GetEntries();
int gensize;
float genpt[MAX_GEN];
float geny[MAX_GEN];
float genpdgId[MAX_GEN];
float genisSignal[MAX_GEN];
ntGen->SetBranchAddress("size", &gensize);
ntGen->SetBranchAddress("pt", genpt);
ntGen->SetBranchAddress("y", geny);
ntGen->SetBranchAddress("pdgId", genpdgId);
ntGen->SetBranchAddress("isSignal", genisSignal);
ntGen->SetBranchStatus("*",0);
ntGen->SetBranchStatus("size",1);
ntGen->SetBranchStatus("pt",1);
ntGen->SetBranchStatus("y",1);
ntGen->SetBranchStatus("pdgId",1);
ntGen->SetBranchStatus("isSignal",1);
TRandom3 *evtrandom = new TRandom3();
evtrandom->SetSeed(iseed);
TH1D* thisevt = new TH1D("thisevt","",50, 5, 6);
TH1D* thisevt2 = new TH1D("thisevt2","",50, 5, 6);
TH1D* GenPt = new TH1D("GenPt","",60, 0, 60);
TH1D *hyield[nBins];
for (int i=0;i<nBins;i++)
	{
//	std::vector<TF1*> pseudo_results;
	TF1* pseudo_results;
	double Nsig = hPtCor->GetBinContent(i+1)*(ptBins[i+1]-ptBins[i]);
	cout<<"Got Nsig:"<<Nsig<<"/"<<hPt->GetBinContent(i+1)*(ptBins[i+1]-ptBins[i])<<endl;
//	double Nbg = bg[i]->Integral(5,6)/0.02;
	double Nbg = bg[i].Integral(5,6)/0.02;
	hyield[i] = new TH1D(Form("hyield%d", i),"",100,Nsig-7*sqrt(Nsig),Nsig+7*sqrt(Nsig));
	TH1D *hMass = new TH1D("hMass","",50,5,6);
	for(int ex = 0; ex < 1; ex++){
//	for(int ex = 0; ex < 20; ex++){
//		for(int nsig = 0; nsig < 1; nsig++){
			cout<<"Pseudo exp:"<<ex<<endl;
			now = time(0);dt = ctime(&now);cout << "The local date and time is: " << dt << endl;
			bool gotsig = false;
			int getthisevt = -1;
			int count = 0;
			int true_yield = 0;
			int _count = 0;
			getthisevt = int(evtrandom->Uniform(ngen));
			do{
				_count++;
//                cout<<"Getting evt: "<<getthisevt<<" count: "<<count <<endl;
				ntGen->Project("GenPt","pt",Form("pdgId==511 && isSignal!=0 && abs(y+0.465)<1.93 && pt>%f && pt<%f", ptBins[i],ptBins[i+1]),"",1,getthisevt);
//				cout<<"pt: "<<GenPt->GetMaximumBin()<< " entry: " << GenPt->GetEntries()<<endl;
				if(GenPt->GetEntries() != 0) {
					gotsig = 1;
				}
/*
				ntGen->GetEntry(getthisevt);				
				for(int g=0; g<gensize; g++){
					if(abs(genpdgId[g])!=511) continue;
					if(genisSignal[g]!=0 && abs(geny[g]+0.465)<1.93 && genpt[g] > ptBins[i] && genpt[g] < ptBins[i+1]) {
						gotsig = 1;
						break;
					}
				}
*/
				//gotsig=1;
				if(gotsig){
					ntMC->Project("thisevt","mass",Form("%s&&pt>%f&&pt<%f",selmc_kpi.Data(),ptBins[i],ptBins[i+1]),"",1,getthisevt);
					ntMC2->Project("thisevt2","mass",Form("%s&&pt>%f&&pt<%f",selmc_pik.Data(),ptBins[i],ptBins[i+1]),"",1,getthisevt);
					hMass->Add(thisevt);
					hMass->Add(thisevt2);
//	                true_yield += thisevt->Integral(5,6);	
//	                true_yield += thisevt2->Integral(5,6);	
					count++;
					gotsig = 0;
					//cout<<"this"<<thisevt->Integral()<<endl;
				}
				//cout<<"count"<<count<<endl;

				//Get following
				getthisevt++;
				//Random
//                getthisevt = int(evtrandom->Uniform(ngen));
				if(getthisevt >= ngen) getthisevt=0;
//			}while(count<500);
			}while(count<(int)Nsig);
			cout<<"cout:"<<count<<endl;
			cout<<"_cout:"<<_count<<endl;
			true_yield = hMass->Integral();
			cout<<"true_yield:"<<true_yield<<endl;
//		}
		for(int nbg = 0; nbg < (int)Nbg; nbg++){
//    	  double bkg=bg[i]->GetRandom();
    	  double bkg=bg[i].GetRandom();
	      hMass->Fill(bkg);
		}
//	    pseudo_results.clear();
		cout<<"Start pseudo exp fit:"<<ex<<endl;
		now = time(0);dt = ctime(&now);cout << "The local date and time is: " << dt << endl;
//		pseudo_results = fit(nt,nt2,ntMC,ntMC2,ptBins[i],ptBins[i+1], 1, hMass, &(sigAndbg[i]));
		pseudo_results = fit(dummynt,dummynt,dummynt,dummynt,ptBins[i],ptBins[i+1], 1, hMass, &(sigAndbg[i]));
		now = time(0);dt = ctime(&now);cout << "The local date and time is: " << dt << endl;
//		cout<<"sig yield:"<<pseudo_results[0]->Integral(5,6)/0.02<<endl;
//		pseudo->Fill(pseudo_results[0]->Integral(5,6)/0.02,i);
		pseudo->Fill(pseudo_results[0].Integral(5,6)/0.02,i,true_yield,pseudo_results[2].GetParError(0)/0.02);
//		hyield[i]->Fill(pseudo_results[0]->Integral(5,6)/0.02);
		hMass->Reset();
	}//pseudo exp loop
//    hyield[i]->Write();
}//pt bin loo[
pseudo->Write();
//Pseudo

  hPt->Write();
  hEff->Write();
  hPtCor->Write();
  hPtSigma->Write();
  outf->Close();
  delete outf;
}
