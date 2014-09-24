#include "utilities.h"
#define MAX_XB 8192
#define MAX_GEN 4096
using namespace std;

// current date/time based on current system
time_t now = time(0);   
char* dt = ctime(&now);
//stdlib qsort
int compare (const void * a, const void * b) {return ( *(int*)a - *(int*)b );}
//my qsort
void quickSort(int arr[], int left, int right) {
      int i = left, j = right;
      int tmp;
      int pivot = arr[(left + right) / 2];
      /* partition */
      while (i <= j) {
            while (arr[i] < pivot)
                  i++;
            while (arr[j] > pivot)
                  j--;
            if (i <= j) {
                  tmp = arr[i];
                  arr[i] = arr[j];
                  arr[j] = tmp;
                  i++;
                  j--;
            }
      };
      /* recursion */
      if (left < j)
            quickSort(arr, left, j);
      if (i < right)
            quickSort(arr, i, right);
}

double luminosity=34.8*1e-3;

//double setparam0=5.85733e+01;
//double setparam1=5.28008e+00;
//double setparam2=5.00000e-02;
//double setparam3=1.02740e+02;
//double setparam4=-1.62289e+01;
//double setparam5=0.00000e+00;
//double setparam6=1.43809e-01;
//double setparam7=1.86128e-01;
//double setparam8=1.42964e-02;

double setparam0=100.;
double setparam1=5.28;
double setparam2=0.03;
double fixparam1=5.279;
double setparam3=0.03;
double fixparam2=0.04;

//TString inputdata="/mnt/hadoop/cms/store/user/jwang/nt_20140403_PAMuon_HIRun2013_PromptrecoAndRereco_v1_MuonMatching_EvtBase_skim.root";
//TString inputmc="/mnt/hadoop/cms/store/user/jwang/nt_BoostedMC_20140403_Kstar_TriggerMatchingMuon_EvtBase_skim.root";
//TString inputdata="/mnt/hadoop/cms/store/user/jwang/nt_20140411_PAMuon_HIRun2013_PromptrecoAndRereco_v1_MuonMatching_EvtBase_skim.root";
//TString inputmc="/mnt/hadoop/cms/store/user/jwang/nt_BoostedMC_20140411_Kstar_TriggerMatchingMuon_EvtBase_skim.root";
//TString inputdata="/mnt/hadoop/cms/store/user/jwang/nt_20140418_PAMuon_HIRun2013_PromptrecoAndRereco_v1_MuonMatching_EvtBase_skim.root";
//TString inputmc="/mnt/hadoop/cms/store/user/jwang/nt_BoostedMC_20140418_Kstar_TriggerMatchingMuon_EvtBase_skim.root";
TString inputdata="/mnt/hadoop/cms/store/user/jwang/nt_20140427_PAMuon_HIRun2013_PromptrecoAndRereco_v1_MuonMatching_EvtBase_skim.root";
TString inputmc="/mnt/hadoop/cms/store/user/jwang/nt_BoostedMC_20140427_Kstar_TriggerMatchingMuon_EvtBase_skim.root";

//TString cut_kpi="(HLT_PAMu3_v1)&&abs(mumumass-3.096916)<0.15&&chi2cl>0.15&&(d0)/d0Err>8.1&&cos(dtheta)>-0.44&&TMath::Abs((trk2Dxy)/trk2D0Err)>0.81&&abs(tktkmass-0.89591)<0.14&&mass>5&&mass<6";
//TString cut_kpi="(HLT_PAMu3_v1)&&abs(mumumass-3.096916)<0.15&&chi2cl>0.069&&(d0)/d0Err>3.2&&cos(dtheta)>-0.7&&abs(tktkmass-0.89591)<0.13&&mass>5&&mass<6&&isbesttktkmass";
TString cut_kpi="(HLT_PAMu3_v1)&&abs(mumumass-3.096916)<0.15&&mass>5&&mass<6&& isbestchi2&&trk1Pt>0.7&&trk2Pt>0.7&&chi2cl>1.65e-01&&(d0/d0Err)>4.16&&cos(dtheta)>7.50e-01&&abs(tktkmass-0.89594)<2.33e-01";

TString seldata_kpi=Form("abs(y+0.465)<1.93&&%s",cut_kpi.Data());
TString seldata_2y_kpi=Form("((Run>=210498&&Run<=211256&&abs(y+0.465)<1.93)||(Run>=211313&&Run<=211631&&abs(y-0.465)<1.93))&&%s",cut_kpi.Data());

TString selmc_kpi=Form("abs(y+0.465)<1.93&&(gen==23333||gen==41000)&&%s",cut_kpi.Data());

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
TF1* fit(TTree *nt,TTree *nt2, TTree *ntMC, TTree *ntMC2, double ptmin,double ptmax, bool dopseudo, TH1D *pseudo, TF1 *pseudotf1){//pe
   //cout<<cut.Data()<<endl;
   static int count=0;
   static TF1 return_th1[4];//pe
   TF1 fitfn_MC;//pe
   count++;
   TCanvas *c= new TCanvas(Form("c%d",count),"",600,600);
   TH1D *h = new TH1D(Form("h%d",count),"",50,5,6);
//   TH1D *h2 = new TH1D(Form("h2%d",count),"",50,5,6);
//   TH1D *hBck = new TH1D(Form("hBck%d",count),"",50,5,6);
   
   TH1D *hMC = new TH1D(Form("hMC%d",count),"",50,5,6);
//   TH1D *hMC2 = new TH1D(Form("hMC2%d",count),"",50,5,6);
   // Fit function
//   TF1 *f = new TF1(Form("f%d",count),"[0]*([7]*Gaus(x,[1],[2])/(sqrt(2*3.14159)*[2])+(1-[7])*Gaus(x,[1],[8])/(sqrt(2*3.14159)*[8]))+[3]+[4]*x+[6]*(38.42*Gaus(x,5.25,0.03473)+15.04*Gaus(x,5.25,0.1121)+104.3*Gaus(x,5.026,0.0935))");
   TString iNP="6.71675e+00*Gaus(x,5.30142e+00,8.42680e-02)/(sqrt(2*3.14159)*8.42680e-02)+4.06744e+01*Gaus(x,5.00954e+00,8.11305e-02)/(sqrt(2*3.14159)*8.11305e-02)+5.99974e-01*(2.376716*Gaus(x,5.640619,0.095530)/(sqrt(2*3.14159)*0.095530)+3.702342*Gaus(x,5.501706,0.046222)/(sqrt(2*3.14159)*0.046222))+1.31767e-01*(61.195688*Gaus(x,5.127566,0.087439)/(sqrt(2*3.14159)*0.087439)+58.943919*Gaus(x,5.246471,0.041983)/(sqrt(2*3.14159)*0.041983))";
   TF1 *f = new TF1(Form("f%d",count),"[0]*([7]*Gaus(x,[1],[2])/(sqrt(2*3.14159)*[2])+(1-[7])*Gaus(x,[1],[8])/(sqrt(2*3.14159)*[8]))+[3]+[4]*x+[6]*("+iNP+")");

   if(!dopseudo){//pe
     nt->Project(Form("h%d",count),"mass",Form("%s&&pt>%f&&pt<%f",seldata_2y_kpi.Data(),ptmin,ptmax));   
     ntMC->Project(Form("hMC%d",count),"mass",Form("%s&&pt>%f&&pt<%f",seldata_kpi.Data(),ptmin,ptmax));   
//     h->Add(h2);
//     hMC->Add(hMC2);
    }//pe
	if(dopseudo) hMC = (TH1D*)pseudo->Clone(Form("h%d",count));//pe
	if(dopseudo) h = (TH1D*)pseudo->Clone(Form("h%d",count));//pe

//   nt->Project(Form("hBck%d",count),"mass",Form("%s&&pt>%f&&pt<%f&&(gen==23333||gen==41000)",seldata.Data(),ptmin,ptmax));   
//   nt2->Project(Form("hBck%d",count),"mass",Form("%s&&pt>%f&&pt<%f&&(gen==23333||gen==41000)",seldata.Data(),ptmin,ptmax));   

//   cout <<"nsig = "<<hBck->GetEntries();
   clean0(h);
   h->Draw();
   f->SetParLimits(4,-1000,0);
   f->SetParLimits(2,0.01,0.05);
   f->SetParLimits(8,0.01,0.1);
   f->SetParLimits(7,0,1);
   f->SetParLimits(6,0,1000);

//   f->SetParameter(0,setparam0);
//   f->SetParameter(1,setparam1);
//   f->SetParameter(2,setparam2);
//   f->SetParameter(3,setparam3);
//   f->SetParameter(4,setparam4);
//   f->SetParameter(5,setparam5);
//   f->SetParameter(6,setparam6);
//   f->SetParameter(7,setparam7);
//   f->SetParameter(8,setparam8);

   f->SetParameter(0,setparam0);
   f->SetParameter(1,setparam1);
   f->SetParameter(2,setparam2);
   f->SetParameter(8,setparam3);
   f->FixParameter(1,fixparam1);

if(dopseudo){//pe
   f->SetParameter(0,pseudotf1->GetParameter(0));
   f->SetParameter(1,pseudotf1->GetParameter(1));
   f->SetParameter(2,pseudotf1->GetParameter(2));
   f->SetParameter(3,pseudotf1->GetParameter(3));
   f->SetParameter(4,pseudotf1->GetParameter(4));
   f->SetParameter(5,pseudotf1->GetParameter(5));
   f->SetParameter(6,pseudotf1->GetParameter(6));
   f->SetParameter(7,pseudotf1->GetParameter(7));
   f->SetParameter(8,pseudotf1->GetParameter(8));
}//pe

   f->FixParameter(1,fixparam1);
   h->GetEntries();

if(!dopseudo){//pe
   hMC->Fit(Form("f%d",count),"q","",5,6);
   hMC->Fit(Form("f%d",count),"q","",5,6);
   f->ReleaseParameter(1);
/*
   hMC->Fit(Form("f%d",count),"q","",5,6);
   hMC->Fit(Form("f%d",count),"q","",5,6);
   hMC->Fit(Form("f%d",count),"q","",5,6);
//   hMC->Fit(Form("f%d",count),"","",5,6);
   hMC->Fit(Form("f%d",count),"m","",5,6);
*/
   hMC->Fit(Form("f%d",count),"L q","",5,6);
   hMC->Fit(Form("f%d",count),"L q","",5,6);
   hMC->Fit(Form("f%d",count),"L q","",5,6);
//   hMC->Fit(Form("f%d",count),"","",5,6);
   hMC->Fit(Form("f%d",count),"L m","",5,6);
}//pe
   
   fitfn_MC = *f;
   f->FixParameter(1,f->GetParameter(1));
   f->FixParameter(2,f->GetParameter(2));
   f->FixParameter(7,f->GetParameter(7));
   f->FixParameter(8,f->GetParameter(8));

   h->Fit(Form("f%d",count),"q","",5,6);
   h->Fit(Form("f%d",count),"q","",5,6);
   f->ReleaseParameter(1);
/*
   h->Fit(Form("f%d",count),"q","",5,6);
   h->Fit(Form("f%d",count),"q","",5,6);
   h->Fit(Form("f%d",count),"q","",5,6);
//   h->Fit(Form("f%d",count),"","",5,6);
   h->Fit(Form("f%d",count),"m","",5,6);
*/
   h->Fit(Form("f%d",count),"L q","",5,6);
   h->Fit(Form("f%d",count),"L q","",5,6);
   h->Fit(Form("f%d",count),"L q","",5,6);
//   h->Fit(Form("f%d",count),"","",5,6);
   h->Fit(Form("f%d",count),"L m","",5,6);

   h->SetMarkerSize(0.8);
   h->SetMarkerStyle(20);
   cout <<h->GetEntries()<<endl;
   cout <<hMC->GetEntries()<<endl;

   TF1 *all_background = new TF1(Form("all_background%d",count),"[0]+[1]*x+[2]*("+iNP+")");
   all_background->SetParameter(0,f->GetParameter(3));
   all_background->SetParameter(1,f->GetParameter(4));
   all_background->SetParameter(2,f->GetParameter(6));
   all_background->SetRange(5,6);

   // function for background shape plotting. take the fit result from f
   TF1 *background = new TF1(Form("background%d",count),"[0]+[1]*x");
   background->SetParameter(0,f->GetParameter(3));
   background->SetParameter(1,f->GetParameter(4));
   background->SetParameter(2,f->GetParameter(5));
   background->SetParameter(3,f->GetParameter(6));
   background->SetLineColor(4);
   background->SetRange(5,6);
   background->SetLineStyle(2);
   
   // function for signal shape plotting. take the fit result from f
//   TF1 *Bkpi = new TF1(Form("fBkpi"),"[0]*(38.42*Gaus(x,5.25,0.03473)+15.04*Gaus(x,5.25,0.1121)+104.3*Gaus(x,5.026,0.0935))");
   TF1 *Bkpi = new TF1(Form("fBkpi",count),"[0]*("+iNP+")");
   Bkpi->SetParameter(0,f->GetParameter(6));
   Bkpi->SetLineColor(kGreen+1);
   Bkpi->SetFillColor(kGreen+1);
   Bkpi->SetRange(5.00,5.45);
   Bkpi->SetLineStyle(1);
   Bkpi->SetFillStyle(3005);

   // function for signal shape plotting. take the fit result from f
   TF1 *mass = new TF1(Form("fmass"),"[0]*([3]*Gaus(x,[1],[2])/(sqrt(2*3.14159)*[2])+(1-[3])*Gaus(x,[1],[4])/(sqrt(2*3.14159)*[4]))");
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
   return_th1[0] = *mass;//pe
   return_th1[1] = *all_background;//pe
   return_th1[2] = *f;//pe
   return_th1[3] = fitfn_MC;//pe
   return return_th1;//pe
}

//void fitB0(TString infname="",bool doweight = 1)
void fitB0(int iseed = 1)
{
  gRandom->SetSeed(iseed);
//return;
  TString infname="";
//  if (doweight==0) weight="1";
  if (infname=="") infname=inputdata.Data();
  TFile *inf = new TFile(infname.Data());
  TTree *nt = (TTree*) inf->Get("ntKstar");
  TTree *nt2 = (TTree*) inf->Get("ntKstar");

  TFile *infMC = new TFile(inputmc.Data());
  TTree *ntGen = (TTree*)infMC->Get("ntGen");
  TTree *ntMC = (TTree*)infMC->Get("ntKstar");
  TTree *ntMC2 = (TTree*)infMC->Get("ntKstar");

  //nt->SetAlias("LD",LDalias.Data());
  //nt2->SetAlias("LD",LDalias.Data());
  //ntMC->SetAlias("LD",LDalias.Data());
  //ntMC2->SetAlias("LD",LDalias.Data());
  
  const int nBins = 3;
  double ptBins[nBins+1] = {10,15,20,60};
//  const int nBins = 1;
//  double ptBins[nBins+1] = {10,60};
  TH1D *hPt = new TH1D("hPt","",nBins,ptBins);
  TH1D *hRecoTruth = new TH1D("hRecoTruth","",nBins,ptBins);
  TH1D *hRecoTruth2 = new TH1D("hRecoTruth2","",nBins,ptBins);
  TH1D *hPtMC = new TH1D("hPtMC","",nBins,ptBins);
  TH1D *hPtMC2 = new TH1D("hPtMC2","",nBins,ptBins);
  TH1D *hPtGen = new TH1D("hPtGen","",nBins,ptBins);

  TH1D* dummy = new TH1D("dummy","",50, 5, 6);//pe
  TF1* dummytf1;//pe
  TF1 mass[nBins];//pe
  TF1 bg[nBins];//pe
  TF1 sigAndbg[nBins];//pe
  TF1 sigAndbg_MC[nBins];//pe
  TF1 f;//pe
  TF1* results_tf1;//pe
  //std::vector<TF1*> bg;
  //TF1* f;
  //std::vector<TF1*> results_tf1;
  //TF1 results_tf1[2];
  for (int i=0;i<nBins;i++)
    {
//      TF1 *f = fit(nt,nt2,ntMC,ntMC2,ptBins[i],ptBins[i+1]);
	  now = time(0);dt = ctime(&now);cout << "The local date and time is: " << dt << endl;
      results_tf1 = fit(nt,nt2,ntMC,ntMC2,ptBins[i],ptBins[i+1], 0, dummy, dummytf1);//pe
	  now = time(0);dt = ctime(&now);cout << "The local date and time is: " << dt << endl;
//      f = results_tf1[0];
//      f = (TF1*)results_tf1[0].Clone();
//      bg.push_back(results_tf1[1]);
      f = results_tf1[0];//pe
      mass[i] = results_tf1[0];//pe
      bg[i] = results_tf1[1];//pe
      sigAndbg[i] = results_tf1[2];//pe
      sigAndbg_MC[i] = results_tf1[3];//pe
      double yield = f.Integral(5,6)/0.02;//pe
      double yieldErr = f.Integral(5,6)/0.02*f.GetParError(0)/f.GetParameter(0);//pe
      hPt->SetBinContent(i+1,yield/(ptBins[i+1]-ptBins[i]));
      hPt->SetBinError(i+1,yieldErr/(ptBins[i+1]-ptBins[i]));
  }  
  
//  TCanvas *c=  new TCanvas("cResult","",600,600);
  hPt->SetXTitle("B^{0} p_{T} (GeV/c)");
  hPt->SetYTitle("Uncorrected B^{0} dN/dp_{T}");
  hPt->Sumw2();
  hPt->Draw();

  ntMC->Project("hPtMC","pt",TCut(weight)*(TCut(seldata_kpi.Data())&&"(gen==23333||gen==41000)"));
  nt->Project("hRecoTruth","pt",TCut(seldata_kpi.Data())&&"(gen==23333||gen==41000)");
//  hPtMC->Add(hPtMC2);

  nt->Project("hRecoTruth","pt",TCut(seldata_kpi.Data())&&"(gen==23333||gen==41000)");
//  hRecoTruth->Add(hRecoTruth2);
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
//  TCanvas *cCor=  new TCanvas("cCorResult","",600,600);
  hPtCor->SetYTitle("Correctd B^{0} dN/dp_{T}");
  hPtCor->Draw();
  hPtGen->Draw("same hist");

  TH1D *hPtSigma= (TH1D*)hPtCor->Clone("hPtSigma");
  hPtSigma->Scale(1./(2*luminosity));
  hPtSigma->SetYTitle("d#sigma/dp_{T} (B^{0}) ");

//  TCanvas *cSigma=  new TCanvas("cSigma","",600,600);

  hPtSigma->Draw();

  TFile *outf = new TFile("ResultsBzero/SigmaBzero.root","recreate");
  outf->cd();

  //Pseudo
  TNtuple* pseudo = new TNtuple("pseudo","","yield:bin:true_yield:yield_err:data_obs:obs_err");                 
  int ngen = ntGen->GetEntries();
  TRandom3 *evtrandom = new TRandom3();
  evtrandom->SetSeed(iseed);
  TH1D* thisevt = new TH1D("thisevt","",50, 5, 6);
  TH1D* thisevt2 = new TH1D("thisevt2","",50, 5, 6);
  TH1D* GenPt = new TH1D("GenPt","",60, 0, 60);
  TH1D *hyield[nBins];
  bool GenSigFromPDF = 0;
  bool PoissonSig = 1;
  for (int i=0;i<nBins;i++)
  	{
  //	std::vector<TF1*> pseudo_results;
  	TF1* pseudo_results;
  	double extracted_sig = hPt->GetBinContent(i+1)*(ptBins[i+1]-ptBins[i]);
  	double Nsig = hPtCor->GetBinContent(i+1)*(ptBins[i+1]-ptBins[i]);
    if(PoissonSig) Nsig = evtrandom->Poisson(Nsig);
  	cout<<"Got Nsig:"<<Nsig<<"/"<<hPt->GetBinContent(i+1)*(ptBins[i+1]-ptBins[i])<<endl;
  	cout<<"Got extracted:"<<extracted_sig<<endl;
  //	double Nbg = bg[i]->Integral(5,6)/0.02;
  	double Nbg = bg[i].Integral(5,6)/0.02;
  	hyield[i] = new TH1D(Form("hyield%d", i),"",100,Nsig-7*sqrt(Nsig),Nsig+7*sqrt(Nsig));
  	TH1D *hMass = new TH1D("hMass","",50,5,6);
//  	for(int ex = 0; ex < 0; ex++){
  	for(int ex = 0; ex < 20; ex++){
  		cout<<"Pseudo exp:"<<ex<<endl;
  		now = time(0);dt = ctime(&now);cout << "The local date and time is: " << dt << endl;
  		bool gotsig = false;
  		int getthisevt = -1;
  		int count = 0;
  		int true_yield = 0;
  		int _count = 0;
  		int myran[50000];
  		if(!GenSigFromPDF){
  			do{
  				getthisevt = int(evtrandom->Uniform(ngen));
  				myran[0] = getthisevt;
  			    for(int s=1; s<50000; s++){
  		    	    myran[s] = int(evtrandom->Uniform(ngen));
  //					myran[s] = getthisevt++;
  	                if(getthisevt >= ngen) getthisevt=0;
  			    }
  //			    qsort (myran, 50000, sizeof(int), compare);
  				quickSort(myran, 0, 49999);
  
  				for(int s=0; s<50000; s++){
  					getthisevt = myran[s];
  					_count++;
  //					cout<<"Getting evt: "<<getthisevt<<" count: "<<count <<endl;
  					ntGen->Project("GenPt","pt",Form("abs(pdgId)==511 && isSignal!=0 && abs(y+0.465)<1.93 && pt>%f && pt<%f", ptBins[i],ptBins[i+1]),"",1,getthisevt);
  					if(GenPt->GetEntries() != 0) {
  						gotsig = 1;
  					}
  
  					//gotsig=1;
  					if(gotsig){
  						ntMC->Project("thisevt","mass",Form("%s&&pt>%f&&pt<%f",selmc_kpi.Data(),ptBins[i],ptBins[i+1]),"",1,getthisevt);
//  						ntMC->Project("thisevt","mass",Form("%s&&pt>%f&&pt<%f&&mass<5.38&&mass>5.18",selmc_kpi.Data(),ptBins[i],ptBins[i+1]),"",1,getthisevt);
  						hMass->Add(thisevt);
//  						hMass->Add(thisevt2);
  						count++;
  						gotsig = 0;
  						//cout<<"this"<<thisevt->Integral()<<endl;
  					}
  					if(count==(int)Nsig) break;
  				}
  //			}while(count<500);
  			}while(count<(int)Nsig);
  		}
  		cout<<"cout:"<<count<<endl;
  		cout<<"_cout:"<<_count<<endl;
  		true_yield = hMass->Integral();
  		cout<<"true_yield:"<<true_yield<<endl;
        gRandom = new TRandom3(iseed);                                                                                                                                                                      
  		if(GenSigFromPDF){
  			for(int nsig = 0; nsig < (int)extracted_sig; nsig++){
      		  double sig=mass[i].GetRandom();
  	    	  hMass->Fill(sig);
  			}
  		}
  		hMass->Write();
  		for(int nbg = 0; nbg < (int)Nbg; nbg++){
  //    	  double bkg=bg[i]->GetRandom();
      	  double bkg=bg[i].GetRandom();
  	      hMass->Fill(bkg);
  		}
  //	    pseudo_results.clear();
  		cout<<"Start pseudo exp fit:"<<ex<<endl;
  		now = time(0);dt = ctime(&now);cout << "The local date and time is: " << dt << endl;
  //		pseudo_results = fit(nt,nt2,ntMC,ntMC2,ptBins[i],ptBins[i+1], 1, hMass, &(sigAndbg[i]));
  		pseudo_results = fit(nt,nt2,ntMC,ntMC2,ptBins[i],ptBins[i+1], 1, hMass, &(sigAndbg_MC[i]));
  		now = time(0);dt = ctime(&now);cout << "The local date and time is: " << dt << endl;
  //		cout<<"sig yield:"<<pseudo_results[0]->Integral(5,6)/0.02<<endl;
  //		pseudo->Fill(pseudo_results[0]->Integral(5,6)/0.02,i);
  		pseudo->Fill(pseudo_results[0].Integral(5,6)/0.02,i,true_yield,pseudo_results[2].GetParError(0)/0.02,hPt->GetBinContent(i+1)*(ptBins[i+1]-ptBins[i]),hPt->GetBinError(i+1)*(ptBins[i+1]-ptBins[i]));
  //		hyield[i]->Fill(pseudo_results[0]->Integral(5,6)/0.02);
  		hMass->Reset();
  	}//pseudo exp loop
  //    hyield[i]->Write();
  }//pt bin loop
  pseudo->Write();
  //Pseudo

  hPt->Write();
  hEff->Write();
  hPtCor->Write();
  hPtSigma->Write();
  outf->Close();
  delete outf;
}
