#include "utilities.h"
// current date/time based on current system
time_t now = time(0);
char* dt = ctime(&now);
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

double luminosity=34*1e-3;
double setparam0=100.;
double setparam1=5.37;
double setparam2=0.02;
double setparam3=0.03;
double fixparam1=5.37;

//svmit02
//TString inputdata="/data/bmeson/data/nt_20140403_PAMuon_HIRun2013_PromptrecoAndRereco_v1_MuonMatching_EvtBase_skim.root";
//TString inputmc="/data/bmeson/MC/nt_BoostedMC_20140403_Phi_TriggerMatchingMuon_EvtBase_skim.root";
//cgate
//TString inputdata="/mnt/hadoop/cms/store/user/jwang/nt_20140411_PAMuon_HIRun2013_PromptrecoAndRereco_v1_MuonMatching_EvtBase_skim.root";
//TString inputmc="/mnt/hadoop/cms/store/user/jwang/nt_BoostedMC_20140411_Phi_TriggerMatchingMuon_EvtBase_skim.root";
//TString inputdata="/mnt/hadoop/cms/store/user/jwang/nt_20140418_PAMuon_HIRun2013_PromptrecoAndRereco_v1_MuonMatching_EvtBase_skim.root";
//TString inputmc="/mnt/hadoop/cms/store/user/jwang/nt_BoostedMC_20140424_Phi_TriggerMatchingMuon_EvtBase_skim.root";
TString inputdata="/mnt/hadoop/cms/store/user/jwang/nt_20140427_PAMuon_HIRun2013_PromptrecoAndRereco_v1_MuonMatching_EvtBase_skim.root";
TString inputmc="/mnt/hadoop/cms/store/user/jwang/nt_BoostedMC_20140427_Phi_TriggerMatchingMuon_EvtBase_skim.root";

//TString cut="(HLT_PAMu3_v1)&&chi2cl>0.135&&(d0)/d0Err>3.24&&dtheta<2.03&&TMath::Abs((trk1Dxy/trk1D0Err))>0.792&&TMath::Abs(trk2Dxy/trk2D0Err)>0.742&&TMath::Abs(tktkmass-1.020)<1.02e-02";
//TString cut="(HLT_PAMu3_v1)&&abs(mumumass-3.096916)<0.15&&chi2cl>3.8e-02 && (d0/d0Err)>3.2e+00 && cos(dtheta)>7.4e-01 && abs(trk2Dxy/trk2D0Err)>1.3e+00 && abs(tktkmass-1.01944)<1.5e-02&&mass>5&&mass<6";
//TString cut="(HLT_PAMu3_v1)&&abs(mumumass-3.096916)<0.15&&chi2cl>2.7e-02&&(d0/d0Err)>2.8&&cos(dtheta)>-0.8&&abs(tktkmass-1.01944)<0.014&&mass>5&&mass<6&&trk1Pt>0.7&&trk2Pt>0.7";
TString cut="(HLT_PAMu3_v1)&&abs(mumumass-3.096916)<0.15&&mass>5&&mass<6&& isbestchi2&&trk1Pt>0.7&&trk2Pt>0.7&& chi2cl>3.71e-02&&(d0/d0Err)>3.37&&cos(dtheta)>2.60e-01&&abs(tktkmass-1.019455)<1.55e-02";

TString seldata=Form("abs(y+0.465)<1.93&&%s",cut.Data());
TString seldata_2y=Form("((Run>=210498&&Run<=211256&&abs(y+0.465)<1.93)||(Run>=211313&&Run<=211631&&abs(y-0.465)<1.93))&&%s",cut.Data());
TString selmc=Form("abs(y+0.465)<1.93&&gen==23333&&%s",cut.Data());
TString selmcgen="abs(y+0.465)<1.93&&isSignal>0";

TString weight = "27.493+pt*(-0.218769)";

void clean0(TH1D *h)
{
   for (int i=1;i<=h->GetNbinsX();i++)
   {
      if (h->GetBinContent(i)==0) h->SetBinError(i,1);
   }
}


//TF1 *fit(TTree *nt,TTree *ntMC,double ptmin,double ptmax)
TF1* fit(TTree *nt, TTree *ntMC, double ptmin,double ptmax, bool dopseudo, TH1D *pseudo, TF1 *pseudotf1){//pe
   //cout<<cut.Data()<<endl;
   static int count=0;
   static TF1 return_th1[4];//pe
   TF1 fitfn_MC;//pe
   count++;
   TCanvas *c= new TCanvas(Form("c%d",count),"",600,600);
   TH1D *h = new TH1D(Form("h%d",count),"",50,5,6);
   TH1D *hMC = new TH1D(Form("hMC%d",count),"",50,5,6);
   // Fit function
//   TF1 *f = new TF1(Form("f%d",count),"[0]*([7]*Gaus(x,[1],[2])/(sqrt(2*3.14159)*[2])+(1-[7])*Gaus(x,[1],[8])/(sqrt(2*3.14159)*[8]))+[3]+[4]*x+[5]*x*x");
   TString iNP="Gaus(x,5.36800e+00,5.77122e-02)/(sqrt(2*3.14159)*abs(5.77122e-02))";
   TF1 *f = new TF1(Form("f%d",count),"[0]*([7]*Gaus(x,[1],[2])/(sqrt(2*3.14159)*[2])+(1-[7])*Gaus(x,[1],[8])/(sqrt(2*3.14159)*[8]))+ [3]+[4]*x+[5]*x*x + [11]*("+iNP+")");
   if(!dopseudo){//pe
     nt->Project(Form("h%d",count),"mass",Form("%s&&pt>%f&&pt<%f",seldata_2y.Data(),ptmin,ptmax));   
     ntMC->Project(Form("hMC%d",count),"mass",Form("%s&&pt>%f&&pt<%f",seldata.Data(),ptmin,ptmax));   
   }//pe
   if(dopseudo) hMC = (TH1D*)pseudo->Clone(Form("h%d",count));//pe
   if(dopseudo) h = (TH1D*)pseudo->Clone(Form("h%d",count));//pe
   clean0(h);
   h->Draw();
   f->SetParLimits(4,-1000,0);
   f->SetParLimits(2,0.01,0.05);
   f->SetParLimits(8,0.01,0.05);
   f->SetParLimits(7,0,1);
   f->SetParameter(0,setparam0);
   f->SetParameter(1,setparam1);
   f->SetParameter(2,setparam2);
   f->SetParameter(8,setparam3);
   f->FixParameter(1,fixparam1);
   f->SetParLimits(11,0.,1000);
   f->SetParameter(11,10);
   f->FixParameter(11,0.);
   h->GetEntries();

   if(dopseudo){//pe
      f->SetParameter(0,pseudotf1->GetParameter(0));//pe
      f->SetParameter(1,pseudotf1->GetParameter(1));//pe
      f->SetParameter(2,pseudotf1->GetParameter(2));//pe                                                                                                                                                    
      f->SetParameter(3,pseudotf1->GetParameter(3));//pe
      f->SetParameter(4,pseudotf1->GetParameter(4));//pe
      f->SetParameter(5,pseudotf1->GetParameter(5));//pe
//      f->SetParameter(6,pseudotf1->GetParameter(6));//pe
      f->SetParameter(7,pseudotf1->GetParameter(7));//pe
      f->SetParameter(8,pseudotf1->GetParameter(8));//pe
   }//pe
   
   if(!dopseudo){//pe
     hMC->Fit(Form("f%d",count),"q","",5,6);
     hMC->Fit(Form("f%d",count),"q","",5,6);
     f->ReleaseParameter(1);
     hMC->Fit(Form("f%d",count),"L q","",5,6);
     hMC->Fit(Form("f%d",count),"L q","",5,6);
     hMC->Fit(Form("f%d",count),"L q","",5,6);
     hMC->Fit(Form("f%d",count),"L m","",5,6);
   }//pe   
   fitfn_MC = *f;//pe

   f->FixParameter(1,f->GetParameter(1));
   f->FixParameter(2,f->GetParameter(2));
   f->FixParameter(7,f->GetParameter(7));
   f->FixParameter(8,f->GetParameter(8));
   
   h->Fit(Form("f%d",count),"q","",5,6);
   h->Fit(Form("f%d",count),"q","",5,6);
   f->ReleaseParameter(1);
   h->Fit(Form("f%d",count),"L q","",5,6);
   h->Fit(Form("f%d",count),"L q","",5,6);
   h->Fit(Form("f%d",count),"L q","",5,6);
   h->Fit(Form("f%d",count),"L m","",5,6);
   h->SetMarkerSize(0.8);
   h->SetMarkerStyle(20);
   cout <<h->GetEntries()<<endl;

   TF1 *all_background = new TF1(Form("all_background%d",count),"[3]+[4]*x+[5]*x*x + [11]*("+iNP+")");
   all_background->SetParameter(3,f->GetParameter(3));
   all_background->SetParameter(4,f->GetParameter(4));
   all_background->SetParameter(5,f->GetParameter(5));
   all_background->SetParameter(11,f->GetParameter(11));
   all_background->SetRange(5,6);

   // function for background shape plotting. take the fit result from f
   TF1 *background = new TF1(Form("background%d",count),"[0]+[1]*x+[2]*x*x");
   background->SetParameter(0,f->GetParameter(3));
   background->SetParameter(1,f->GetParameter(4));
   background->SetParameter(2,f->GetParameter(5));
   background->SetParameter(3,f->GetParameter(6));
   background->SetLineColor(4);
   background->SetRange(5,6);
   background->SetLineStyle(2);
   
   // function for signal shape plotting. take the fit result from f
//   TF1 *Bkpi = new TF1(Form("fBkpi",count),"[0]*(1.24e2*Gaus(x,5.107,0.02987)+1.886e2*Gaus(x,5.0116,5.546e-2))");
   TF1 *Bkpi = new TF1(Form("fBkpi",count),"[0]*("+iNP+")");
   Bkpi->SetParameter(0,f->GetParameter(5));
   Bkpi->SetLineColor(kGreen+1);
   Bkpi->SetFillColor(kGreen+1);
   Bkpi->SetRange(5.00,5.28);
   Bkpi->SetLineStyle(1);
   Bkpi->SetFillStyle(3004);

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
   h->SetTitleOffset(1.,"Y");
   h->SetAxisRange(0,h->GetMaximum()*1.2,"Y");
//   Bkpi->Draw("same");
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
//   leg->AddEntry(Bkpi,"Non-prompt J/#psi","f");
   leg->Draw();
   TLegend *leg2 = myLegend(0.44,0.33,0.89,0.50);
   leg2->AddEntry(h,"B meson","");
   leg2->AddEntry(h,Form("M_{B}=%.2f #pm %.2f MeV/c^{2}",mass->GetParameter(1)*1000.,mass->GetParError(1)*1000.),"");
   leg2->AddEntry(h,Form("N_{B}=%.0f #pm %.0f", yield, yieldErr),"");
   leg2->Draw();

   c->SaveAs(Form("ResultsBs/BMass-%d.C",count));
   c->SaveAs(Form("ResultsBs/BMass-%d.gif",count));
   //c->SaveAs(Form("ResultsBs/BMass-%d.eps",count));
   //c->SaveAs(Form("ResultsBs/BMass-%d.pdf",count));

   //return mass;
   return_th1[0] = *mass;//pe
   return_th1[1] = *all_background;//pe
   return_th1[2] = *f;//pe
   return_th1[3] = fitfn_MC;//pe
   return return_th1;//pe
}

//void fitBs(TString infname="",bool doweight = 1)
void fitBs(int iseed = 1)//pe
{
//   if (doweight==0) weight="1";//pe
   TString infname="";//pe
   if (infname=="") infname=inputdata.Data();
   TFile *inf = new TFile(infname.Data());
   TTree *nt = (TTree*) inf->Get("ntphi");

   TFile *infMC = new TFile(inputmc.Data());
   TTree *ntGen = (TTree*)infMC->Get("ntGen");
   TTree *ntMC = (TTree*)infMC->Get("ntphi");

   //nt->SetAlias("LD",LDalias.Data());
   //ntMC->SetAlias("LD",LDalias.Data());

   const int nBins = 1;
   double ptBins[nBins+1] = {10,60};
   TH1D *hPt = new TH1D("hPt","",nBins,ptBins);
   TH1D *hPtMC = new TH1D("hPtMC","",nBins,ptBins);
   TH1D *hPtGen = new TH1D("hPtGen","",nBins,ptBins);
   ntMC->Project("hPtMC","pt",TCut(weight)*(selmc.Data()));
   ntGen->Project("hPtGen","pt",TCut(weight)*(selmcgen.Data()));

  TH1D* dummy = new TH1D("dummy","",50, 5, 6);//pe
  TF1* dummytf1;//pe
  TF1 mass[nBins];//pe
  TF1 bg[nBins];//pe
  TF1 sigAndbg[nBins];//pe
  TF1 sigAndbg_MC[nBins];//pe
  TF1 f;//pe
  TF1* results_tf1;//pe
   for (int i=0;i<nBins;i++)
   {
      //TF1 *f = fit(nt,ntMC,ptBins[i],ptBins[i+1]);
      now = time(0);dt = ctime(&now);cout << "The local date and time is: " << dt << endl;
      results_tf1 = fit(nt,ntMC,ptBins[i],ptBins[i+1], 0, dummy, dummytf1);//pe
      now = time(0);dt = ctime(&now);cout << "The local date and time is: " << dt << endl;
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

   divideBinWidth(hPtMC);
   divideBinWidth(hPtGen);

   TCanvas *c=  new TCanvas("cResult","",600,600);
   hPt->SetXTitle("B_{s} p_{T} (GeV/c)");
   hPt->SetYTitle("Uncorrected B_{s} dN/dp_{T}");
   hPt->Draw();
   
  
  hPtMC->Sumw2();
  TH1D *hEff = (TH1D*)hPtMC->Clone("hEff");
  hPtMC->Sumw2();
  hEff->Divide(hPtGen);
  
  TH1D *hPtCor = (TH1D*)hPt->Clone("hPtCor");
  hPtCor->Divide(hEff);
  TCanvas *cCor=  new TCanvas("cCorResult","",600,600);
  hPtCor->SetYTitle("Correctd B_{s} dN/dp_{T}");
  hPtCor->Draw();

  TH1D *hPtSigma= (TH1D*)hPtCor->Clone("hPtSigma");
  double BRchain=2.89977e-5;
  hPtSigma->Scale(1./(2*luminosity*BRchain));
  hPtSigma->SetYTitle("d#sigma/dp_{T} (B_{s})");

  TCanvas *cSigma=  new TCanvas("cSigma","",600,600);

  hPtSigma->Draw();
  
  TFile *outf = new TFile("ResultsBs/SigmaBs.root","recreate");
  outf->cd();

  //Pseudo//////////////////                                                                                                                                                                                
  TNtuple* pseudo = new TNtuple("pseudo","","yield:bin:true_yield:yield_err:data_obs:obs_err");
  int ngen = ntGen->GetEntries();
  TRandom3 *evtrandom = new TRandom3();
  evtrandom->SetSeed(iseed);
  TH1D* thisevt = new TH1D("thisevt","",50, 5, 6);
  TH1D* GenPt = new TH1D("GenPt","",60, 0, 60);
  TH1D *hyield[nBins];
  bool GenSigFromPDF = 0;
  bool PoissonSig = 1; 
  for (int i=0;i<nBins;i++)
    { 
    TF1* pseudo_results;
    double extracted_sig = hPt->GetBinContent(i+1)*(ptBins[i+1]-ptBins[i]);
    double Nsig = hPtCor->GetBinContent(i+1)*(ptBins[i+1]-ptBins[i]);
    if(PoissonSig) Nsig = evtrandom->Poisson(Nsig);
    cout<<"Got Nsig:"<<Nsig<<"/"<<hPt->GetBinContent(i+1)*(ptBins[i+1]-ptBins[i])<<endl;
    cout<<"Got extracted:"<<extracted_sig<<endl;
    double Nbg = bg[i].Integral(5,6)/0.02;
    hyield[i] = new TH1D(Form("hyield%d", i),"",100,Nsig-7*sqrt(Nsig),Nsig+7*sqrt(Nsig));
    TH1D *hMass = new TH1D("hMass","",50,5,6);
//    for(int ex = 0; ex < 1; ex++){
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
                //myran[s] = getthisevt++;
                if(getthisevt >= ngen) getthisevt=0;
            }
            quickSort(myran, 0, 49999);
  
            for(int s=0; s<50000; s++){
                getthisevt = myran[s];
                _count++;
                //cout<<"Getting evt: "<<getthisevt<<" count: "<<count <<endl;
                ntGen->Project("GenPt","pt",Form("isSignal>0 && abs(y+0.465)<1.93 && pt>%f && pt<%f", ptBins[i],ptBins[i+1]),"",1,getthisevt);
                if(GenPt->GetEntries() != 0) {
                    gotsig = 1;
                }
                //gotsig=1;
                if(gotsig){
                    ntMC->Project("thisevt","mass",Form("%s&&pt>%f&&pt<%f",selmc.Data(),ptBins[i],ptBins[i+1]),"",1,getthisevt);
                    hMass->Add(thisevt);
                    count++;
                    gotsig = 0;
                    //cout<<"this"<<thisevt->Integral()<<endl;
                }
                if(count==(int)Nsig) break;
            }
          //}while(count<500);
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
          double bkg=bg[i].GetRandom();
          hMass->Fill(bkg);
        }
        cout<<"Start pseudo exp fit:"<<ex<<endl;
        now = time(0);dt = ctime(&now);cout << "The local date and time is: " << dt << endl;
        pseudo_results = fit(nt,ntMC,ptBins[i],ptBins[i+1], 1, hMass, &(sigAndbg_MC[i]));
        now = time(0);dt = ctime(&now);cout << "The local date and time is: " << dt << endl;
        pseudo->Fill(pseudo_results[0].Integral(5,6)/0.02,i,true_yield,pseudo_results[2].GetParError(0)/0.02,hPt->GetBinContent(i+1)*(ptBins[i+1]-ptBins[i]),hPt->GetBinError(i+1)*(ptBins[i+1]-ptBins[i]));
        hMass->Reset();
    }//pseudo exp loop
  }//pt bin loop
  pseudo->Write();
  //Pseudo//////////////////                                         

  hPt->Write();
  hEff->Write();
  hPtCor->Write();
  hPtSigma->Write();
  outf->Close();
  delete outf;
}
