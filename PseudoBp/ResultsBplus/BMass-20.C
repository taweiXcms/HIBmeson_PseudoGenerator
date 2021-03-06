{
//=========Macro generated from canvas: c20/
//=========  (Fri Apr 25 19:06:15 2014) by ROOT version5.32/00
   TCanvas *c20 = new TCanvas("c20", "",0,0,600,600);
   gStyle->SetOptFit(1);
   gStyle->SetOptStat(0);
   c20->Range(0,0,1,1);
   c20->SetFillColor(0);
   c20->SetBorderMode(0);
   c20->SetBorderSize(2);
   c20->SetTickx(1);
   c20->SetTicky(1);
   c20->SetLeftMargin(0.13);
   c20->SetRightMargin(0.05);
   c20->SetTopMargin(0.05);
   c20->SetBottomMargin(0.13);
   c20->SetFrameFillStyle(0);
   c20->SetFrameBorderMode(0);
   
   TH1D *h20 = new TH1D("h20","",50,5,6);
   h20->SetBinContent(1,32);
   h20->SetBinContent(2,26);
   h20->SetBinContent(3,26);
   h20->SetBinContent(4,22);
   h20->SetBinContent(5,18);
   h20->SetBinContent(6,24);
   h20->SetBinContent(7,20);
   h20->SetBinContent(8,15);
   h20->SetBinContent(9,21);
   h20->SetBinContent(10,18);
   h20->SetBinContent(11,25);
   h20->SetBinContent(12,33);
   h20->SetBinContent(13,62);
   h20->SetBinContent(14,120);
   h20->SetBinContent(15,114);
   h20->SetBinContent(16,45);
   h20->SetBinContent(17,37);
   h20->SetBinContent(18,27);
   h20->SetBinContent(19,20);
   h20->SetBinContent(20,16);
   h20->SetBinContent(21,22);
   h20->SetBinContent(22,20);
   h20->SetBinContent(23,23);
   h20->SetBinContent(24,21);
   h20->SetBinContent(25,13);
   h20->SetBinContent(26,21);
   h20->SetBinContent(27,19);
   h20->SetBinContent(28,17);
   h20->SetBinContent(29,14);
   h20->SetBinContent(30,16);
   h20->SetBinContent(31,13);
   h20->SetBinContent(32,20);
   h20->SetBinContent(33,16);
   h20->SetBinContent(34,19);
   h20->SetBinContent(35,19);
   h20->SetBinContent(36,15);
   h20->SetBinContent(37,11);
   h20->SetBinContent(38,20);
   h20->SetBinContent(39,12);
   h20->SetBinContent(40,20);
   h20->SetBinContent(41,13);
   h20->SetBinContent(42,14);
   h20->SetBinContent(43,12);
   h20->SetBinContent(44,17);
   h20->SetBinContent(45,13);
   h20->SetBinContent(46,16);
   h20->SetBinContent(47,19);
   h20->SetBinContent(48,15);
   h20->SetBinContent(49,15);
   h20->SetBinContent(50,12);
   h20->SetMinimum(0);
   h20->SetMaximum(144);
   h20->SetEntries(1218);
   h20->SetStats(0);
   
   TF1 *f20 = new TF1("f20","[0]*([7]*Gaus(x,[1],[2])/(sqrt(2*3.14159)*[2])+(1-[7])*Gaus(x,[1],[8])/(sqrt(2*3.14159)*[8]))+[3]+[4]*x+[5]*(4.18604e+01*Gaus(x,4.97611e+00,8.81611e-02)/(sqrt(2*3.14159)*8.81611e-02)+6.72820e+00*Gaus(x,5.10752e+00,2.66663e-02)/(sqrt(2*3.14159)*2.66663e-02)+1.93889e+00*Gaus(x,5.33740e+00,3.52905e-02)/(sqrt(2*3.14159)*3.52905e-02))",5,6);
   f20->SetFillColor(19);
   f20->SetFillStyle(0);
   f20->SetMarkerStyle(20);
   f20->SetLineColor(2);
   f20->SetLineWidth(1);
   f20->SetChisquare(30.08777);
   f20->SetNDF(44);
   f20->GetXaxis()->SetLabelFont(42);
   f20->GetXaxis()->SetLabelOffset(0.007);
   f20->GetXaxis()->SetLabelSize(0.05);
   f20->GetXaxis()->SetTitleSize(0.06);
   f20->GetXaxis()->SetTitleOffset(0.9);
   f20->GetXaxis()->SetTitleFont(42);
   f20->GetYaxis()->SetLabelFont(42);
   f20->GetYaxis()->SetLabelOffset(0.007);
   f20->GetYaxis()->SetLabelSize(0.05);
   f20->GetYaxis()->SetTitleSize(0.06);
   f20->GetYaxis()->SetTitleOffset(1.05);
   f20->GetYaxis()->SetTitleFont(42);
   f20->SetParameter(0,6.125927);
   f20->SetParError(0,0.4130116);
   f20->SetParLimits(0,0,0);
   f20->SetParameter(1,5.27776);
   f20->SetParError(1,0.001725009);
   f20->SetParLimits(1,0,0);
   f20->SetParameter(2,0.05);
   f20->SetParError(2,0);
   f20->SetParLimits(2,0.05,0.05);
   f20->SetParameter(3,55.22554);
   f20->SetParError(3,0.6231947);
   f20->SetParLimits(3,0,0);
   f20->SetParameter(4,-6.822182);
   f20->SetParError(4,0.1118586);
   f20->SetParLimits(4,-1000,0);
   f20->SetParameter(5,0.02112239);
   f20->SetParError(5,0.01187967);
   f20->SetParLimits(5,0,1000);
   f20->SetParameter(6,0);
   f20->SetParError(6,2.4);
   f20->SetParLimits(6,0,0);
   f20->SetParameter(7,0.2397448);
   f20->SetParError(7,0);
   f20->SetParLimits(7,0.2397448,0.2397448);
   f20->SetParameter(8,0.01876219);
   f20->SetParError(8,0);
   f20->SetParLimits(8,0.01876219,0.01876219);
   h20->GetListOfFunctions()->Add(f20);
   h20->SetLineStyle(0);
   h20->SetMarkerStyle(24);
   h20->SetMarkerSize(0.8);
   h20->GetXaxis()->SetTitle("M_{B} (GeV/c^{2})");
   h20->GetXaxis()->CenterTitle(true);
   h20->GetXaxis()->SetLabelFont(42);
   h20->GetXaxis()->SetLabelOffset(0.007);
   h20->GetXaxis()->SetLabelSize(0.05);
   h20->GetXaxis()->SetTitleSize(0.06);
   h20->GetXaxis()->SetTitleOffset(0.9);
   h20->GetXaxis()->SetTitleFont(42);
   h20->GetYaxis()->SetTitle("Entries / (20 MeV/c^{2})");
   h20->GetYaxis()->CenterTitle(true);
   h20->GetYaxis()->SetLabelFont(42);
   h20->GetYaxis()->SetLabelOffset(0.007);
   h20->GetYaxis()->SetLabelSize(0.05);
   h20->GetYaxis()->SetTitleSize(0.06);
   h20->GetYaxis()->SetTitleFont(42);
   h20->GetZaxis()->SetLabelFont(42);
   h20->GetZaxis()->SetLabelOffset(0.007);
   h20->GetZaxis()->SetLabelSize(0.05);
   h20->GetZaxis()->SetTitleSize(0.06);
   h20->GetZaxis()->SetTitleFont(42);
   h20->Draw("e");
   
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
   fBkpi->SetParameter(0,0.02112239);
   fBkpi->SetParError(0,0);
   fBkpi->SetParLimits(0,0,0);
   fBkpi->Draw("same");
   
   TF1 *background20 = new TF1("background20","[0]+[1]*x",5,6);
   background20->SetFillColor(19);
   background20->SetFillStyle(0);
   background20->SetMarkerStyle(20);
   background20->SetLineColor(4);
   background20->SetLineWidth(1);
   background20->SetLineStyle(2);
   background20->GetXaxis()->SetLabelFont(42);
   background20->GetXaxis()->SetLabelOffset(0.007);
   background20->GetXaxis()->SetLabelSize(0.05);
   background20->GetXaxis()->SetTitleSize(0.06);
   background20->GetXaxis()->SetTitleOffset(0.9);
   background20->GetXaxis()->SetTitleFont(42);
   background20->GetYaxis()->SetLabelFont(42);
   background20->GetYaxis()->SetLabelOffset(0.007);
   background20->GetYaxis()->SetLabelSize(0.05);
   background20->GetYaxis()->SetTitleSize(0.06);
   background20->GetYaxis()->SetTitleOffset(1.05);
   background20->GetYaxis()->SetTitleFont(42);
   background20->SetParameter(0,55.22554);
   background20->SetParError(0,0);
   background20->SetParLimits(0,0,0);
   background20->SetParameter(1,-6.822182);
   background20->SetParError(1,0);
   background20->SetParLimits(1,0,0);
   background20->Draw("same");
   
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
   fmass->SetParameter(0,6.125927);
   fmass->SetParError(0,0.4130116);
   fmass->SetParLimits(0,0,0);
   fmass->SetParameter(1,5.27776);
   fmass->SetParError(1,0.001725009);
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
   
   TF1 *f20 = new TF1("f20","[0]*([7]*Gaus(x,[1],[2])/(sqrt(2*3.14159)*[2])+(1-[7])*Gaus(x,[1],[8])/(sqrt(2*3.14159)*[8]))+[3]+[4]*x+[5]*(4.18604e+01*Gaus(x,4.97611e+00,8.81611e-02)/(sqrt(2*3.14159)*8.81611e-02)+6.72820e+00*Gaus(x,5.10752e+00,2.66663e-02)/(sqrt(2*3.14159)*2.66663e-02)+1.93889e+00*Gaus(x,5.33740e+00,3.52905e-02)/(sqrt(2*3.14159)*3.52905e-02))",0,1);
   f20->SetFillColor(19);
   f20->SetFillStyle(0);
   f20->SetMarkerStyle(20);
   f20->SetLineColor(2);
   f20->SetLineWidth(1);
   f20->SetChisquare(30.08777);
   f20->SetNDF(44);
   f20->GetXaxis()->SetLabelFont(42);
   f20->GetXaxis()->SetLabelOffset(0.007);
   f20->GetXaxis()->SetLabelSize(0.05);
   f20->GetXaxis()->SetTitleSize(0.06);
   f20->GetXaxis()->SetTitleOffset(0.9);
   f20->GetXaxis()->SetTitleFont(42);
   f20->GetYaxis()->SetLabelFont(42);
   f20->GetYaxis()->SetLabelOffset(0.007);
   f20->GetYaxis()->SetLabelSize(0.05);
   f20->GetYaxis()->SetTitleSize(0.06);
   f20->GetYaxis()->SetTitleOffset(1.05);
   f20->GetYaxis()->SetTitleFont(42);
   f20->SetParameter(0,6.125927);
   f20->SetParError(0,0.4130116);
   f20->SetParLimits(0,0,0);
   f20->SetParameter(1,5.27776);
   f20->SetParError(1,0.001725009);
   f20->SetParLimits(1,0,0);
   f20->SetParameter(2,0.05);
   f20->SetParError(2,0);
   f20->SetParLimits(2,0.05,0.05);
   f20->SetParameter(3,55.22554);
   f20->SetParError(3,0.6231947);
   f20->SetParLimits(3,0,0);
   f20->SetParameter(4,-6.822182);
   f20->SetParError(4,0.1118586);
   f20->SetParLimits(4,-1000,0);
   f20->SetParameter(5,0.02112239);
   f20->SetParError(5,0.01187967);
   f20->SetParLimits(5,0,1000);
   f20->SetParameter(6,0);
   f20->SetParError(6,2.4);
   f20->SetParLimits(6,0,0);
   f20->SetParameter(7,0.2397448);
   f20->SetParError(7,0);
   f20->SetParLimits(7,0.2397448,0.2397448);
   f20->SetParameter(8,0.01876219);
   f20->SetParError(8,0);
   f20->SetParLimits(8,0.01876219,0.01876219);
   f20->Draw("same");
   
   TLegend *leg = new TLegend(2.121996e-314,0,0,0,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextFont(62);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(19);
   leg->SetFillStyle(0);
   TLegendEntry *entry=leg->AddEntry("h20","CMS Preliminary","");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("h20","p+Pb #sqrt{s_{NN}}= 5.02 TeV","");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("h20","10<p_{T}^{B}<15 GeV/c","");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("h20","Data","pl");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("f20","Fit","l");
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
   entry=leg->AddEntry("background20","Combinatorial Background","l");
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
   
   leg = new TLegend(8.766126e-308,8.770081e-308,8.772899e-308,8.783204e-308,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextFont(62);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(19);
   leg->SetFillStyle(0);
   entry=leg->AddEntry("h20","B meson","");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("h20","M_{B}=5277.76 #pm 1.73 MeV/c^{2}","");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("h20","N_{B}=306 #pm 21","");
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
   c20->Modified();
   c20->cd();
   c20->SetSelected(c20);
}
