{
//=========Macro generated from canvas: c14/
//=========  (Fri Apr 25 19:03:10 2014) by ROOT version5.32/00
   TCanvas *c14 = new TCanvas("c14", "",0,0,600,600);
   gStyle->SetOptFit(1);
   gStyle->SetOptStat(0);
   c14->Range(0,0,1,1);
   c14->SetFillColor(0);
   c14->SetBorderMode(0);
   c14->SetBorderSize(2);
   c14->SetTickx(1);
   c14->SetTicky(1);
   c14->SetLeftMargin(0.13);
   c14->SetRightMargin(0.05);
   c14->SetTopMargin(0.05);
   c14->SetBottomMargin(0.13);
   c14->SetFrameFillStyle(0);
   c14->SetFrameBorderMode(0);
   
   TH1D *h14 = new TH1D("h14","",50,5,6);
   h14->SetBinContent(1,29);
   h14->SetBinContent(2,27);
   h14->SetBinContent(3,27);
   h14->SetBinContent(4,20);
   h14->SetBinContent(5,26);
   h14->SetBinContent(6,20);
   h14->SetBinContent(7,16);
   h14->SetBinContent(8,24);
   h14->SetBinContent(9,16);
   h14->SetBinContent(10,23);
   h14->SetBinContent(11,18);
   h14->SetBinContent(12,31);
   h14->SetBinContent(13,67);
   h14->SetBinContent(14,136);
   h14->SetBinContent(15,145);
   h14->SetBinContent(16,57);
   h14->SetBinContent(17,32);
   h14->SetBinContent(18,28);
   h14->SetBinContent(19,26);
   h14->SetBinContent(20,20);
   h14->SetBinContent(21,18);
   h14->SetBinContent(22,21);
   h14->SetBinContent(23,14);
   h14->SetBinContent(24,17);
   h14->SetBinContent(25,24);
   h14->SetBinContent(26,19);
   h14->SetBinContent(27,18);
   h14->SetBinContent(28,11);
   h14->SetBinContent(29,15);
   h14->SetBinContent(30,23);
   h14->SetBinContent(31,24);
   h14->SetBinContent(32,18);
   h14->SetBinContent(33,15);
   h14->SetBinContent(34,19);
   h14->SetBinContent(35,22);
   h14->SetBinContent(36,11);
   h14->SetBinContent(37,14);
   h14->SetBinContent(38,13);
   h14->SetBinContent(39,13);
   h14->SetBinContent(40,18);
   h14->SetBinContent(41,13);
   h14->SetBinContent(42,22);
   h14->SetBinContent(43,10);
   h14->SetBinContent(44,12);
   h14->SetBinContent(45,12);
   h14->SetBinContent(46,12);
   h14->SetBinContent(47,12);
   h14->SetBinContent(48,15);
   h14->SetBinContent(49,8);
   h14->SetBinContent(50,9);
   h14->SetMinimum(0);
   h14->SetMaximum(174);
   h14->SetEntries(1260);
   h14->SetStats(0);
   
   TF1 *f14 = new TF1("f14","[0]*([7]*Gaus(x,[1],[2])/(sqrt(2*3.14159)*[2])+(1-[7])*Gaus(x,[1],[8])/(sqrt(2*3.14159)*[8]))+[3]+[4]*x+[5]*(4.18604e+01*Gaus(x,4.97611e+00,8.81611e-02)/(sqrt(2*3.14159)*8.81611e-02)+6.72820e+00*Gaus(x,5.10752e+00,2.66663e-02)/(sqrt(2*3.14159)*2.66663e-02)+1.93889e+00*Gaus(x,5.33740e+00,3.52905e-02)/(sqrt(2*3.14159)*3.52905e-02))",5,6);
   f14->SetFillColor(19);
   f14->SetFillStyle(0);
   f14->SetMarkerStyle(20);
   f14->SetLineColor(2);
   f14->SetLineWidth(1);
   f14->SetChisquare(40.78335);
   f14->SetNDF(44);
   f14->GetXaxis()->SetLabelFont(42);
   f14->GetXaxis()->SetLabelOffset(0.007);
   f14->GetXaxis()->SetLabelSize(0.05);
   f14->GetXaxis()->SetTitleSize(0.06);
   f14->GetXaxis()->SetTitleOffset(0.9);
   f14->GetXaxis()->SetTitleFont(42);
   f14->GetYaxis()->SetLabelFont(42);
   f14->GetYaxis()->SetLabelOffset(0.007);
   f14->GetYaxis()->SetLabelSize(0.05);
   f14->GetYaxis()->SetTitleSize(0.06);
   f14->GetYaxis()->SetTitleOffset(1.05);
   f14->GetYaxis()->SetTitleFont(42);
   f14->SetParameter(0,7.365841);
   f14->SetParError(0,0.4436137);
   f14->SetParLimits(0,0,0);
   f14->SetParameter(1,5.279821);
   f14->SetParError(1,0.0014924);
   f14->SetParLimits(1,0,0);
   f14->SetParameter(2,0.05);
   f14->SetParError(2,0);
   f14->SetParLimits(2,0.05,0.05);
   f14->SetParameter(3,78.1134);
   f14->SetParError(3,0.610513);
   f14->SetParLimits(3,0,0);
   f14->SetParameter(4,-11.00503);
   f14->SetParError(4,0.1093133);
   f14->SetParLimits(4,-1000,0);
   f14->SetParameter(5,0.009890862);
   f14->SetParError(5,0.01187727);
   f14->SetParLimits(5,0,1000);
   f14->SetParameter(6,0);
   f14->SetParError(6,1);
   f14->SetParLimits(6,0,0);
   f14->SetParameter(7,0.2397448);
   f14->SetParError(7,0);
   f14->SetParLimits(7,0.2397448,0.2397448);
   f14->SetParameter(8,0.01876219);
   f14->SetParError(8,0);
   f14->SetParLimits(8,0.01876219,0.01876219);
   h14->GetListOfFunctions()->Add(f14);
   h14->SetLineStyle(0);
   h14->SetMarkerStyle(24);
   h14->SetMarkerSize(0.8);
   h14->GetXaxis()->SetTitle("M_{B} (GeV/c^{2})");
   h14->GetXaxis()->CenterTitle(true);
   h14->GetXaxis()->SetLabelFont(42);
   h14->GetXaxis()->SetLabelOffset(0.007);
   h14->GetXaxis()->SetLabelSize(0.05);
   h14->GetXaxis()->SetTitleSize(0.06);
   h14->GetXaxis()->SetTitleOffset(0.9);
   h14->GetXaxis()->SetTitleFont(42);
   h14->GetYaxis()->SetTitle("Entries / (20 MeV/c^{2})");
   h14->GetYaxis()->CenterTitle(true);
   h14->GetYaxis()->SetLabelFont(42);
   h14->GetYaxis()->SetLabelOffset(0.007);
   h14->GetYaxis()->SetLabelSize(0.05);
   h14->GetYaxis()->SetTitleSize(0.06);
   h14->GetYaxis()->SetTitleFont(42);
   h14->GetZaxis()->SetLabelFont(42);
   h14->GetZaxis()->SetLabelOffset(0.007);
   h14->GetZaxis()->SetLabelSize(0.05);
   h14->GetZaxis()->SetTitleSize(0.06);
   h14->GetZaxis()->SetTitleFont(42);
   h14->Draw("e");
   
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
   fBkpi->SetParameter(0,0.009890862);
   fBkpi->SetParError(0,0);
   fBkpi->SetParLimits(0,0,0);
   fBkpi->Draw("same");
   
   TF1 *background14 = new TF1("background14","[0]+[1]*x",5,6);
   background14->SetFillColor(19);
   background14->SetFillStyle(0);
   background14->SetMarkerStyle(20);
   background14->SetLineColor(4);
   background14->SetLineWidth(1);
   background14->SetLineStyle(2);
   background14->GetXaxis()->SetLabelFont(42);
   background14->GetXaxis()->SetLabelOffset(0.007);
   background14->GetXaxis()->SetLabelSize(0.05);
   background14->GetXaxis()->SetTitleSize(0.06);
   background14->GetXaxis()->SetTitleOffset(0.9);
   background14->GetXaxis()->SetTitleFont(42);
   background14->GetYaxis()->SetLabelFont(42);
   background14->GetYaxis()->SetLabelOffset(0.007);
   background14->GetYaxis()->SetLabelSize(0.05);
   background14->GetYaxis()->SetTitleSize(0.06);
   background14->GetYaxis()->SetTitleOffset(1.05);
   background14->GetYaxis()->SetTitleFont(42);
   background14->SetParameter(0,78.1134);
   background14->SetParError(0,0);
   background14->SetParLimits(0,0,0);
   background14->SetParameter(1,-11.00503);
   background14->SetParError(1,0);
   background14->SetParLimits(1,0,0);
   background14->Draw("same");
   
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
   fmass->SetParameter(0,7.365841);
   fmass->SetParError(0,0.4436137);
   fmass->SetParLimits(0,0,0);
   fmass->SetParameter(1,5.279821);
   fmass->SetParError(1,0.0014924);
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
   
   TF1 *f14 = new TF1("f14","[0]*([7]*Gaus(x,[1],[2])/(sqrt(2*3.14159)*[2])+(1-[7])*Gaus(x,[1],[8])/(sqrt(2*3.14159)*[8]))+[3]+[4]*x+[5]*(4.18604e+01*Gaus(x,4.97611e+00,8.81611e-02)/(sqrt(2*3.14159)*8.81611e-02)+6.72820e+00*Gaus(x,5.10752e+00,2.66663e-02)/(sqrt(2*3.14159)*2.66663e-02)+1.93889e+00*Gaus(x,5.33740e+00,3.52905e-02)/(sqrt(2*3.14159)*3.52905e-02))",0,1);
   f14->SetFillColor(19);
   f14->SetFillStyle(0);
   f14->SetMarkerStyle(20);
   f14->SetLineColor(2);
   f14->SetLineWidth(1);
   f14->SetChisquare(40.78335);
   f14->SetNDF(44);
   f14->GetXaxis()->SetLabelFont(42);
   f14->GetXaxis()->SetLabelOffset(0.007);
   f14->GetXaxis()->SetLabelSize(0.05);
   f14->GetXaxis()->SetTitleSize(0.06);
   f14->GetXaxis()->SetTitleOffset(0.9);
   f14->GetXaxis()->SetTitleFont(42);
   f14->GetYaxis()->SetLabelFont(42);
   f14->GetYaxis()->SetLabelOffset(0.007);
   f14->GetYaxis()->SetLabelSize(0.05);
   f14->GetYaxis()->SetTitleSize(0.06);
   f14->GetYaxis()->SetTitleOffset(1.05);
   f14->GetYaxis()->SetTitleFont(42);
   f14->SetParameter(0,7.365841);
   f14->SetParError(0,0.4436137);
   f14->SetParLimits(0,0,0);
   f14->SetParameter(1,5.279821);
   f14->SetParError(1,0.0014924);
   f14->SetParLimits(1,0,0);
   f14->SetParameter(2,0.05);
   f14->SetParError(2,0);
   f14->SetParLimits(2,0.05,0.05);
   f14->SetParameter(3,78.1134);
   f14->SetParError(3,0.610513);
   f14->SetParLimits(3,0,0);
   f14->SetParameter(4,-11.00503);
   f14->SetParError(4,0.1093133);
   f14->SetParLimits(4,-1000,0);
   f14->SetParameter(5,0.009890862);
   f14->SetParError(5,0.01187727);
   f14->SetParLimits(5,0,1000);
   f14->SetParameter(6,0);
   f14->SetParError(6,1);
   f14->SetParLimits(6,0,0);
   f14->SetParameter(7,0.2397448);
   f14->SetParError(7,0);
   f14->SetParLimits(7,0.2397448,0.2397448);
   f14->SetParameter(8,0.01876219);
   f14->SetParError(8,0);
   f14->SetParLimits(8,0.01876219,0.01876219);
   f14->Draw("same");
   
   TLegend *leg = new TLegend(2.546395e-313,2.758595e-313,3.182994e-313,5.203447e-90,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextFont(62);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(19);
   leg->SetFillStyle(0);
   TLegendEntry *entry=leg->AddEntry("h14","CMS Preliminary","");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("h14","p+Pb #sqrt{s_{NN}}= 5.02 TeV","");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("h14","10<p_{T}^{B}<15 GeV/c","");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("h14","Data","pl");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("f14","Fit","l");
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
   entry=leg->AddEntry("background14","Combinatorial Background","l");
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
   
   leg = new TLegend(1.120043e-312,5.37353e-316,2.334195e-313,2.758595e-313,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextFont(62);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(19);
   leg->SetFillStyle(0);
   entry=leg->AddEntry("h14","B meson","");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("h14","M_{B}=5279.82 #pm 1.49 MeV/c^{2}","");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("h14","N_{B}=368 #pm 22","");
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
   c14->Modified();
   c14->cd();
   c14->SetSelected(c14);
}
