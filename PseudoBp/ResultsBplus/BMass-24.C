{
//=========Macro generated from canvas: c24/
//=========  (Fri Apr 25 19:08:18 2014) by ROOT version5.32/00
   TCanvas *c24 = new TCanvas("c24", "",0,0,600,600);
   gStyle->SetOptFit(1);
   gStyle->SetOptStat(0);
   c24->Range(0,0,1,1);
   c24->SetFillColor(0);
   c24->SetBorderMode(0);
   c24->SetBorderSize(2);
   c24->SetTickx(1);
   c24->SetTicky(1);
   c24->SetLeftMargin(0.13);
   c24->SetRightMargin(0.05);
   c24->SetTopMargin(0.05);
   c24->SetBottomMargin(0.13);
   c24->SetFrameFillStyle(0);
   c24->SetFrameBorderMode(0);
   
   TH1D *h24 = new TH1D("h24","",50,5,6);
   h24->SetBinContent(1,19);
   h24->SetBinContent(2,17);
   h24->SetBinContent(3,21);
   h24->SetBinContent(4,22);
   h24->SetBinContent(5,25);
   h24->SetBinContent(6,24);
   h24->SetBinContent(7,26);
   h24->SetBinContent(8,28);
   h24->SetBinContent(9,29);
   h24->SetBinContent(10,29);
   h24->SetBinContent(11,34);
   h24->SetBinContent(12,27);
   h24->SetBinContent(13,65);
   h24->SetBinContent(14,119);
   h24->SetBinContent(15,117);
   h24->SetBinContent(16,43);
   h24->SetBinContent(17,22);
   h24->SetBinContent(18,22);
   h24->SetBinContent(19,24);
   h24->SetBinContent(20,21);
   h24->SetBinContent(21,18);
   h24->SetBinContent(22,21);
   h24->SetBinContent(23,16);
   h24->SetBinContent(24,20);
   h24->SetBinContent(25,19);
   h24->SetBinContent(26,15);
   h24->SetBinContent(27,20);
   h24->SetBinContent(28,15);
   h24->SetBinContent(29,17);
   h24->SetBinContent(30,19);
   h24->SetBinContent(31,15);
   h24->SetBinContent(32,17);
   h24->SetBinContent(33,22);
   h24->SetBinContent(34,7);
   h24->SetBinContent(35,6);
   h24->SetBinContent(36,21);
   h24->SetBinContent(37,12);
   h24->SetBinContent(38,12);
   h24->SetBinContent(39,7);
   h24->SetBinContent(40,18);
   h24->SetBinContent(41,14);
   h24->SetBinContent(42,16);
   h24->SetBinContent(43,15);
   h24->SetBinContent(44,19);
   h24->SetBinContent(45,13);
   h24->SetBinContent(46,15);
   h24->SetBinContent(47,11);
   h24->SetBinContent(48,14);
   h24->SetBinContent(49,17);
   h24->SetBinContent(50,14);
   h24->SetMinimum(0);
   h24->SetMaximum(142.8);
   h24->SetEntries(1199);
   h24->SetStats(0);
   
   TF1 *f24 = new TF1("f24","[0]*([7]*Gaus(x,[1],[2])/(sqrt(2*3.14159)*[2])+(1-[7])*Gaus(x,[1],[8])/(sqrt(2*3.14159)*[8]))+[3]+[4]*x+[5]*(4.18604e+01*Gaus(x,4.97611e+00,8.81611e-02)/(sqrt(2*3.14159)*8.81611e-02)+6.72820e+00*Gaus(x,5.10752e+00,2.66663e-02)/(sqrt(2*3.14159)*2.66663e-02)+1.93889e+00*Gaus(x,5.33740e+00,3.52905e-02)/(sqrt(2*3.14159)*3.52905e-02))",5,6);
   f24->SetFillColor(19);
   f24->SetFillStyle(0);
   f24->SetMarkerStyle(20);
   f24->SetLineColor(2);
   f24->SetLineWidth(1);
   f24->SetChisquare(63.58519);
   f24->SetNDF(44);
   f24->GetXaxis()->SetLabelFont(42);
   f24->GetXaxis()->SetLabelOffset(0.007);
   f24->GetXaxis()->SetLabelSize(0.05);
   f24->GetXaxis()->SetTitleSize(0.06);
   f24->GetXaxis()->SetTitleOffset(0.9);
   f24->GetXaxis()->SetTitleFont(42);
   f24->GetYaxis()->SetLabelFont(42);
   f24->GetYaxis()->SetLabelOffset(0.007);
   f24->GetYaxis()->SetLabelSize(0.05);
   f24->GetYaxis()->SetTitleSize(0.06);
   f24->GetYaxis()->SetTitleOffset(1.05);
   f24->GetYaxis()->SetTitleFont(42);
   f24->SetParameter(0,5.97305);
   f24->SetParError(0,0.4100022);
   f24->SetParLimits(0,0,0);
   f24->SetParameter(1,5.276888);
   f24->SetParError(1,0.001680463);
   f24->SetParLimits(1,0,0);
   f24->SetParameter(2,0.05);
   f24->SetParError(2,0);
   f24->SetParLimits(2,0.05,0.05);
   f24->SetParameter(3,79.0631);
   f24->SetParError(3,0.6122624);
   f24->SetParLimits(3,0,0);
   f24->SetParameter(4,-11.10111);
   f24->SetParError(4,0.1096515);
   f24->SetParLimits(4,-1000,0);
   f24->SetParameter(5,2.401801e-09);
   f24->SetParError(5,0.01828643);
   f24->SetParLimits(5,0,1000);
   f24->SetParameter(6,0);
   f24->SetParError(6,2.4);
   f24->SetParLimits(6,0,0);
   f24->SetParameter(7,0.2397448);
   f24->SetParError(7,0);
   f24->SetParLimits(7,0.2397448,0.2397448);
   f24->SetParameter(8,0.01876219);
   f24->SetParError(8,0);
   f24->SetParLimits(8,0.01876219,0.01876219);
   h24->GetListOfFunctions()->Add(f24);
   h24->SetLineStyle(0);
   h24->SetMarkerStyle(24);
   h24->SetMarkerSize(0.8);
   h24->GetXaxis()->SetTitle("M_{B} (GeV/c^{2})");
   h24->GetXaxis()->CenterTitle(true);
   h24->GetXaxis()->SetLabelFont(42);
   h24->GetXaxis()->SetLabelOffset(0.007);
   h24->GetXaxis()->SetLabelSize(0.05);
   h24->GetXaxis()->SetTitleSize(0.06);
   h24->GetXaxis()->SetTitleOffset(0.9);
   h24->GetXaxis()->SetTitleFont(42);
   h24->GetYaxis()->SetTitle("Entries / (20 MeV/c^{2})");
   h24->GetYaxis()->CenterTitle(true);
   h24->GetYaxis()->SetLabelFont(42);
   h24->GetYaxis()->SetLabelOffset(0.007);
   h24->GetYaxis()->SetLabelSize(0.05);
   h24->GetYaxis()->SetTitleSize(0.06);
   h24->GetYaxis()->SetTitleFont(42);
   h24->GetZaxis()->SetLabelFont(42);
   h24->GetZaxis()->SetLabelOffset(0.007);
   h24->GetZaxis()->SetLabelSize(0.05);
   h24->GetZaxis()->SetTitleSize(0.06);
   h24->GetZaxis()->SetTitleFont(42);
   h24->Draw("e");
   
   TF1 *fBkpi = new TF1("fBkpi","[0]*(4.18604e+01*Gaus(x,4.97611e+00,8.81611e-02)/(sqrt(2*3.14159)*8.81611e-02)+6.72820e+00*Gaus(x,5.10752e+00,2.66663e-02)/(sqrt(2*3.14159)*2.66663e-02)+1.93889e+00*Gaus(x,5.33740e+00,3.52905e-02)/(sqrt(2*3.14159)*3.52905e-02))",5,5.28);

   Int_t ci;   // for color index setting
   ci = TColor::GetColor("#00cc00");
   fBkpi->SetFillColor(ci);
   fBkpi->SetFillStyle(3004);
   fBkpi->SetMarkerStyle(20);

   ci = TColor::GetColor("#00cc00");
   fBkpi->SetLineColor(ci);
   fBkpi->SetLineWidth(1);
   fBkpi->GetXaxis()->SetLabelFont(42);
   fBkpi->GetXaxis()->SetLabelOffset(0.007);
   fBkpi->GetXaxis()->SetLabelSize(0.05);
   fBkpi->GetXaxis()->SetTitleSize(0.06);
   fBkpi->GetXaxis()->SetTitleOffset(0.9);
   fBkpi->GetXaxis()->SetTitleFont(42);
   fBkpi->GetYaxis()->SetLabelFont(42);
   fBkpi->GetYaxis()->SetLabelOffset(0.007);
   fBkpi->GetYaxis()->SetLabelSize(0.05);
   fBkpi->GetYaxis()->SetTitleSize(0.06);
   fBkpi->GetYaxis()->SetTitleOffset(1.05);
   fBkpi->GetYaxis()->SetTitleFont(42);
   fBkpi->SetParameter(0,2.401801e-09);
   fBkpi->SetParError(0,0);
   fBkpi->SetParLimits(0,0,0);
   fBkpi->Draw("same");
   
   TF1 *background24 = new TF1("background24","[0]+[1]*x",5,6);
   background24->SetFillColor(19);
   background24->SetFillStyle(0);
   background24->SetMarkerStyle(20);
   background24->SetLineColor(4);
   background24->SetLineWidth(1);
   background24->SetLineStyle(2);
   background24->GetXaxis()->SetLabelFont(42);
   background24->GetXaxis()->SetLabelOffset(0.007);
   background24->GetXaxis()->SetLabelSize(0.05);
   background24->GetXaxis()->SetTitleSize(0.06);
   background24->GetXaxis()->SetTitleOffset(0.9);
   background24->GetXaxis()->SetTitleFont(42);
   background24->GetYaxis()->SetLabelFont(42);
   background24->GetYaxis()->SetLabelOffset(0.007);
   background24->GetYaxis()->SetLabelSize(0.05);
   background24->GetYaxis()->SetTitleSize(0.06);
   background24->GetYaxis()->SetTitleOffset(1.05);
   background24->GetYaxis()->SetTitleFont(42);
   background24->SetParameter(0,79.0631);
   background24->SetParError(0,0);
   background24->SetParLimits(0,0,0);
   background24->SetParameter(1,-11.10111);
   background24->SetParError(1,0);
   background24->SetParLimits(1,0,0);
   background24->Draw("same");
   
   TF1 *fmass = new TF1("fmass","[0]*([3]*Gaus(x,[1],[2])/(sqrt(2*3.14159)*[2])+(1-[3])*Gaus(x,[1],[4])/(sqrt(2*3.14159)*[4]))",5,6);
   fmass->SetFillColor(2);
   fmass->SetFillStyle(3004);
   fmass->SetMarkerStyle(20);
   fmass->SetLineColor(2);
   fmass->SetLineWidth(1);
   fmass->SetLineStyle(2);
   fmass->GetXaxis()->SetLabelFont(42);
   fmass->GetXaxis()->SetLabelOffset(0.007);
   fmass->GetXaxis()->SetLabelSize(0.05);
   fmass->GetXaxis()->SetTitleSize(0.06);
   fmass->GetXaxis()->SetTitleOffset(0.9);
   fmass->GetXaxis()->SetTitleFont(42);
   fmass->GetYaxis()->SetLabelFont(42);
   fmass->GetYaxis()->SetLabelOffset(0.007);
   fmass->GetYaxis()->SetLabelSize(0.05);
   fmass->GetYaxis()->SetTitleSize(0.06);
   fmass->GetYaxis()->SetTitleOffset(1.05);
   fmass->GetYaxis()->SetTitleFont(42);
   fmass->SetParameter(0,5.97305);
   fmass->SetParError(0,0.4100022);
   fmass->SetParLimits(0,0,0);
   fmass->SetParameter(1,5.276888);
   fmass->SetParError(1,0.001680463);
   fmass->SetParLimits(1,0,0);
   fmass->SetParameter(2,0.05);
   fmass->SetParError(2,0);
   fmass->SetParLimits(2,0,0);
   fmass->SetParameter(3,0.2397448);
   fmass->SetParError(3,0);
   fmass->SetParLimits(3,0,0);
   fmass->SetParameter(4,0.01876219);
   fmass->SetParError(4,0);
   fmass->SetParLimits(4,0,0);
   fmass->Draw("same");
   
   TF1 *f24 = new TF1("f24","[0]*([7]*Gaus(x,[1],[2])/(sqrt(2*3.14159)*[2])+(1-[7])*Gaus(x,[1],[8])/(sqrt(2*3.14159)*[8]))+[3]+[4]*x+[5]*(4.18604e+01*Gaus(x,4.97611e+00,8.81611e-02)/(sqrt(2*3.14159)*8.81611e-02)+6.72820e+00*Gaus(x,5.10752e+00,2.66663e-02)/(sqrt(2*3.14159)*2.66663e-02)+1.93889e+00*Gaus(x,5.33740e+00,3.52905e-02)/(sqrt(2*3.14159)*3.52905e-02))",0,1);
   f24->SetFillColor(19);
   f24->SetFillStyle(0);
   f24->SetMarkerStyle(20);
   f24->SetLineColor(2);
   f24->SetLineWidth(1);
   f24->SetChisquare(63.58519);
   f24->SetNDF(44);
   f24->GetXaxis()->SetLabelFont(42);
   f24->GetXaxis()->SetLabelOffset(0.007);
   f24->GetXaxis()->SetLabelSize(0.05);
   f24->GetXaxis()->SetTitleSize(0.06);
   f24->GetXaxis()->SetTitleOffset(0.9);
   f24->GetXaxis()->SetTitleFont(42);
   f24->GetYaxis()->SetLabelFont(42);
   f24->GetYaxis()->SetLabelOffset(0.007);
   f24->GetYaxis()->SetLabelSize(0.05);
   f24->GetYaxis()->SetTitleSize(0.06);
   f24->GetYaxis()->SetTitleOffset(1.05);
   f24->GetYaxis()->SetTitleFont(42);
   f24->SetParameter(0,5.97305);
   f24->SetParError(0,0.4100022);
   f24->SetParLimits(0,0,0);
   f24->SetParameter(1,5.276888);
   f24->SetParError(1,0.001680463);
   f24->SetParLimits(1,0,0);
   f24->SetParameter(2,0.05);
   f24->SetParError(2,0);
   f24->SetParLimits(2,0.05,0.05);
   f24->SetParameter(3,79.0631);
   f24->SetParError(3,0.6122624);
   f24->SetParLimits(3,0,0);
   f24->SetParameter(4,-11.10111);
   f24->SetParError(4,0.1096515);
   f24->SetParLimits(4,-1000,0);
   f24->SetParameter(5,2.401801e-09);
   f24->SetParError(5,0.01828643);
   f24->SetParLimits(5,0,1000);
   f24->SetParameter(6,0);
   f24->SetParError(6,2.4);
   f24->SetParLimits(6,0,0);
   f24->SetParameter(7,0.2397448);
   f24->SetParError(7,0);
   f24->SetParLimits(7,0.2397448,0.2397448);
   f24->SetParameter(8,0.01876219);
   f24->SetParError(8,0);
   f24->SetParLimits(8,0.01876219,0.01876219);
   f24->Draw("same");
   
   TLegend *leg = new TLegend(2.546395e-313,2.758595e-313,3.182994e-313,7.90505e-323,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextFont(62);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(19);
   leg->SetFillStyle(0);
   TLegendEntry *entry=leg->AddEntry("h24","CMS Preliminary","");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("h24","p+Pb #sqrt{s_{NN}}= 5.02 TeV","");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("h24","10<p_{T}^{B}<15 GeV/c","");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("h24","Data","pl");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("f24","Fit","l");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("fmass","Signal","f");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("background24","Combinatorial Background","l");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("fBkpi","Non-prompt J/#psi","f");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   leg->Draw();
   
   leg = new TLegend(8.454935e-308,8.458568e-308,8.466564e-308,8.475188e-308,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextFont(62);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(19);
   leg->SetFillStyle(0);
   entry=leg->AddEntry("h24","B meson","");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("h24","M_{B}=5276.89 #pm 1.68 MeV/c^{2}","");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("h24","N_{B}=299 #pm 21","");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   leg->Draw();
   
   TPaveText *pt = new TPaveText(0.01,0.9341608,0.71,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(2);
   pt->SetFillColor(10);
   TText *text = pt->AddText("[0]*([7]*Gaus(x,[1],[2])/(sqrt(2*3.14159)*[2])+(1-[7])*Gaus(x,[1],[8])/(sqrt(2*3.14159)*[8]))+[3]+[4]*x+[5]*(4.18604e+01*Gaus(x,4.97611e+00,8.81611e-02)/(sqrt(2*3.14159)*8.81611e-02)+6.72820e+00*Gaus(x,5.10752e+00,2.66663e-02)/(sqrt(2*3.14159)*2.66663e-02)+1.93889e+00*Gaus(x,5.33740e+00,3.52905e-02)/(sqrt(2*3.14159)*3.52905e-02))");
   pt->Draw();
   c24->Modified();
   c24->cd();
   c24->SetSelected(c24);
}
