{
//=========Macro generated from canvas: c31/
//=========  (Fri Apr 25 19:11:02 2014) by ROOT version5.32/00
   TCanvas *c31 = new TCanvas("c31", "",0,0,600,600);
   gStyle->SetOptFit(1);
   gStyle->SetOptStat(0);
   c31->Range(0,0,1,1);
   c31->SetFillColor(0);
   c31->SetBorderMode(0);
   c31->SetBorderSize(2);
   c31->SetTickx(1);
   c31->SetTicky(1);
   c31->SetLeftMargin(0.13);
   c31->SetRightMargin(0.05);
   c31->SetTopMargin(0.05);
   c31->SetBottomMargin(0.13);
   c31->SetFrameFillStyle(0);
   c31->SetFrameBorderMode(0);
   
   TH1D *h31 = new TH1D("h31","",50,5,6);
   h31->SetBinContent(1,7);
   h31->SetBinContent(2,4);
   h31->SetBinContent(3,5);
   h31->SetBinContent(4,3);
   h31->SetBinContent(5,5);
   h31->SetBinContent(6,7);
   h31->SetBinContent(7,5);
   h31->SetBinContent(8,4);
   h31->SetBinContent(9,4);
   h31->SetBinContent(10,5);
   h31->SetBinContent(11,5);
   h31->SetBinContent(12,12);
   h31->SetBinContent(13,22);
   h31->SetBinContent(14,50);
   h31->SetBinContent(15,44);
   h31->SetBinContent(16,15);
   h31->SetBinContent(17,9);
   h31->SetBinContent(18,5);
   h31->SetBinContent(19,4);
   h31->SetBinContent(20,5);
   h31->SetBinContent(21,3);
   h31->SetBinContent(22,4);
   h31->SetBinContent(23,4);
   h31->SetBinContent(24,4);
   h31->SetBinContent(25,2);
   h31->SetBinContent(26,4);
   h31->SetBinContent(27,4);
   h31->SetBinContent(28,5);
   h31->SetBinContent(29,3);
   h31->SetBinContent(30,2);
   h31->SetBinContent(31,1);
   h31->SetBinContent(32,1);
   h31->SetBinContent(34,5);
   h31->SetBinContent(35,1);
   h31->SetBinContent(36,2);
   h31->SetBinContent(37,2);
   h31->SetBinContent(38,1);
   h31->SetBinContent(39,2);
   h31->SetBinContent(40,5);
   h31->SetBinContent(41,4);
   h31->SetBinContent(42,3);
   h31->SetBinContent(43,2);
   h31->SetBinContent(44,1);
   h31->SetBinContent(45,3);
   h31->SetBinContent(46,2);
   h31->SetBinContent(47,2);
   h31->SetBinContent(48,1);
   h31->SetBinContent(49,3);
   h31->SetBinError(1,2.645751);
   h31->SetBinError(2,2);
   h31->SetBinError(3,2.236068);
   h31->SetBinError(4,1.732051);
   h31->SetBinError(5,2.236068);
   h31->SetBinError(6,2.645751);
   h31->SetBinError(7,2.236068);
   h31->SetBinError(8,2);
   h31->SetBinError(9,2);
   h31->SetBinError(10,2.236068);
   h31->SetBinError(11,2.236068);
   h31->SetBinError(12,3.464102);
   h31->SetBinError(13,4.690416);
   h31->SetBinError(14,7.071068);
   h31->SetBinError(15,6.63325);
   h31->SetBinError(16,3.872983);
   h31->SetBinError(17,3);
   h31->SetBinError(18,2.236068);
   h31->SetBinError(19,2);
   h31->SetBinError(20,2.236068);
   h31->SetBinError(21,1.732051);
   h31->SetBinError(22,2);
   h31->SetBinError(23,2);
   h31->SetBinError(24,2);
   h31->SetBinError(25,1.414214);
   h31->SetBinError(26,2);
   h31->SetBinError(27,2);
   h31->SetBinError(28,2.236068);
   h31->SetBinError(29,1.732051);
   h31->SetBinError(30,1.414214);
   h31->SetBinError(31,1);
   h31->SetBinError(32,1);
   h31->SetBinError(33,1);
   h31->SetBinError(34,2.236068);
   h31->SetBinError(35,1);
   h31->SetBinError(36,1.414214);
   h31->SetBinError(37,1.414214);
   h31->SetBinError(38,1);
   h31->SetBinError(39,1.414214);
   h31->SetBinError(40,2.236068);
   h31->SetBinError(41,2);
   h31->SetBinError(42,1.732051);
   h31->SetBinError(43,1.414214);
   h31->SetBinError(44,1);
   h31->SetBinError(45,1.732051);
   h31->SetBinError(46,1.414214);
   h31->SetBinError(47,1.414214);
   h31->SetBinError(48,1);
   h31->SetBinError(49,1.732051);
   h31->SetBinError(50,1);
   h31->SetMinimum(0);
   h31->SetMaximum(60);
   h31->SetEntries(296);
   h31->SetStats(0);
   
   TF1 *f31 = new TF1("f31","[0]*([7]*Gaus(x,[1],[2])/(sqrt(2*3.14159)*[2])+(1-[7])*Gaus(x,[1],[8])/(sqrt(2*3.14159)*[8]))+[3]+[4]*x+[5]*(4.18604e+01*Gaus(x,4.97611e+00,8.81611e-02)/(sqrt(2*3.14159)*8.81611e-02)+6.72820e+00*Gaus(x,5.10752e+00,2.66663e-02)/(sqrt(2*3.14159)*2.66663e-02)+1.93889e+00*Gaus(x,5.33740e+00,3.52905e-02)/(sqrt(2*3.14159)*3.52905e-02))",5,6);
   f31->SetFillColor(19);
   f31->SetFillStyle(0);
   f31->SetMarkerStyle(20);
   f31->SetLineColor(2);
   f31->SetLineWidth(1);
   f31->SetChisquare(34.83398);
   f31->SetNDF(44);
   f31->GetXaxis()->SetLabelFont(42);
   f31->GetXaxis()->SetLabelOffset(0.007);
   f31->GetXaxis()->SetLabelSize(0.05);
   f31->GetXaxis()->SetTitleSize(0.06);
   f31->GetXaxis()->SetTitleOffset(0.9);
   f31->GetXaxis()->SetTitleFont(42);
   f31->GetYaxis()->SetLabelFont(42);
   f31->GetYaxis()->SetLabelOffset(0.007);
   f31->GetYaxis()->SetLabelSize(0.05);
   f31->GetYaxis()->SetTitleSize(0.06);
   f31->GetYaxis()->SetTitleOffset(1.05);
   f31->GetYaxis()->SetTitleFont(42);
   f31->SetParameter(0,2.636113);
   f31->SetParError(0,0.2516804);
   f31->SetParLimits(0,0,0);
   f31->SetParameter(1,5.276907);
   f31->SetParError(1,0.002249287);
   f31->SetParLimits(1,0,0);
   f31->SetParameter(2,0.03959756);
   f31->SetParError(2,0);
   f31->SetParLimits(2,0.03959756,0.03959756);
   f31->SetParameter(3,20.5921);
   f31->SetParError(3,0.2535745);
   f31->SetParLimits(3,0,0);
   f31->SetParameter(4,-3.165488);
   f31->SetParError(4,0.04503655);
   f31->SetParLimits(4,-1000,0);
   f31->SetParameter(5,0.004060204);
   f31->SetParError(5,0.005526025);
   f31->SetParLimits(5,0,1000);
   f31->SetParameter(6,0);
   f31->SetParError(6,2.4);
   f31->SetParLimits(6,0,0);
   f31->SetParameter(7,0.331348);
   f31->SetParError(7,0);
   f31->SetParLimits(7,0.331348,0.331348);
   f31->SetParameter(8,0.01638322);
   f31->SetParError(8,0);
   f31->SetParLimits(8,0.01638322,0.01638322);
   h31->GetListOfFunctions()->Add(f31);
   h31->SetLineStyle(0);
   h31->SetMarkerStyle(24);
   h31->SetMarkerSize(0.8);
   h31->GetXaxis()->SetTitle("M_{B} (GeV/c^{2})");
   h31->GetXaxis()->CenterTitle(true);
   h31->GetXaxis()->SetLabelFont(42);
   h31->GetXaxis()->SetLabelOffset(0.007);
   h31->GetXaxis()->SetLabelSize(0.05);
   h31->GetXaxis()->SetTitleSize(0.06);
   h31->GetXaxis()->SetTitleOffset(0.9);
   h31->GetXaxis()->SetTitleFont(42);
   h31->GetYaxis()->SetTitle("Entries / (20 MeV/c^{2})");
   h31->GetYaxis()->CenterTitle(true);
   h31->GetYaxis()->SetLabelFont(42);
   h31->GetYaxis()->SetLabelOffset(0.007);
   h31->GetYaxis()->SetLabelSize(0.05);
   h31->GetYaxis()->SetTitleSize(0.06);
   h31->GetYaxis()->SetTitleFont(42);
   h31->GetZaxis()->SetLabelFont(42);
   h31->GetZaxis()->SetLabelOffset(0.007);
   h31->GetZaxis()->SetLabelSize(0.05);
   h31->GetZaxis()->SetTitleSize(0.06);
   h31->GetZaxis()->SetTitleFont(42);
   h31->Draw("e");
   
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
   fBkpi->SetParameter(0,0.004060204);
   fBkpi->SetParError(0,0);
   fBkpi->SetParLimits(0,0,0);
   fBkpi->Draw("same");
   
   TF1 *background31 = new TF1("background31","[0]+[1]*x",5,6);
   background31->SetFillColor(19);
   background31->SetFillStyle(0);
   background31->SetMarkerStyle(20);
   background31->SetLineColor(4);
   background31->SetLineWidth(1);
   background31->SetLineStyle(2);
   background31->GetXaxis()->SetLabelFont(42);
   background31->GetXaxis()->SetLabelOffset(0.007);
   background31->GetXaxis()->SetLabelSize(0.05);
   background31->GetXaxis()->SetTitleSize(0.06);
   background31->GetXaxis()->SetTitleOffset(0.9);
   background31->GetXaxis()->SetTitleFont(42);
   background31->GetYaxis()->SetLabelFont(42);
   background31->GetYaxis()->SetLabelOffset(0.007);
   background31->GetYaxis()->SetLabelSize(0.05);
   background31->GetYaxis()->SetTitleSize(0.06);
   background31->GetYaxis()->SetTitleOffset(1.05);
   background31->GetYaxis()->SetTitleFont(42);
   background31->SetParameter(0,20.5921);
   background31->SetParError(0,0);
   background31->SetParLimits(0,0,0);
   background31->SetParameter(1,-3.165488);
   background31->SetParError(1,0);
   background31->SetParLimits(1,0,0);
   background31->Draw("same");
   
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
   fmass->SetParameter(0,2.636113);
   fmass->SetParError(0,0.2516804);
   fmass->SetParLimits(0,0,0);
   fmass->SetParameter(1,5.276907);
   fmass->SetParError(1,0.002249287);
   fmass->SetParLimits(1,0,0);
   fmass->SetParameter(2,0.03959756);
   fmass->SetParError(2,0);
   fmass->SetParLimits(2,0,0);
   fmass->SetParameter(3,0.331348);
   fmass->SetParError(3,0);
   fmass->SetParLimits(3,0,0);
   fmass->SetParameter(4,0.01638322);
   fmass->SetParError(4,0);
   fmass->SetParLimits(4,0,0);
   fmass->Draw("same");
   
   TF1 *f31 = new TF1("f31","[0]*([7]*Gaus(x,[1],[2])/(sqrt(2*3.14159)*[2])+(1-[7])*Gaus(x,[1],[8])/(sqrt(2*3.14159)*[8]))+[3]+[4]*x+[5]*(4.18604e+01*Gaus(x,4.97611e+00,8.81611e-02)/(sqrt(2*3.14159)*8.81611e-02)+6.72820e+00*Gaus(x,5.10752e+00,2.66663e-02)/(sqrt(2*3.14159)*2.66663e-02)+1.93889e+00*Gaus(x,5.33740e+00,3.52905e-02)/(sqrt(2*3.14159)*3.52905e-02))",0,1);
   f31->SetFillColor(19);
   f31->SetFillStyle(0);
   f31->SetMarkerStyle(20);
   f31->SetLineColor(2);
   f31->SetLineWidth(1);
   f31->SetChisquare(34.83398);
   f31->SetNDF(44);
   f31->GetXaxis()->SetLabelFont(42);
   f31->GetXaxis()->SetLabelOffset(0.007);
   f31->GetXaxis()->SetLabelSize(0.05);
   f31->GetXaxis()->SetTitleSize(0.06);
   f31->GetXaxis()->SetTitleOffset(0.9);
   f31->GetXaxis()->SetTitleFont(42);
   f31->GetYaxis()->SetLabelFont(42);
   f31->GetYaxis()->SetLabelOffset(0.007);
   f31->GetYaxis()->SetLabelSize(0.05);
   f31->GetYaxis()->SetTitleSize(0.06);
   f31->GetYaxis()->SetTitleOffset(1.05);
   f31->GetYaxis()->SetTitleFont(42);
   f31->SetParameter(0,2.636113);
   f31->SetParError(0,0.2516804);
   f31->SetParLimits(0,0,0);
   f31->SetParameter(1,5.276907);
   f31->SetParError(1,0.002249287);
   f31->SetParLimits(1,0,0);
   f31->SetParameter(2,0.03959756);
   f31->SetParError(2,0);
   f31->SetParLimits(2,0.03959756,0.03959756);
   f31->SetParameter(3,20.5921);
   f31->SetParError(3,0.2535745);
   f31->SetParLimits(3,0,0);
   f31->SetParameter(4,-3.165488);
   f31->SetParError(4,0.04503655);
   f31->SetParLimits(4,-1000,0);
   f31->SetParameter(5,0.004060204);
   f31->SetParError(5,0.005526025);
   f31->SetParLimits(5,0,1000);
   f31->SetParameter(6,0);
   f31->SetParError(6,2.4);
   f31->SetParLimits(6,0,0);
   f31->SetParameter(7,0.331348);
   f31->SetParError(7,0);
   f31->SetParLimits(7,0.331348,0.331348);
   f31->SetParameter(8,0.01638322);
   f31->SetParError(8,0);
   f31->SetParLimits(8,0.01638322,0.01638322);
   f31->Draw("same");
   
   TLegend *leg = new TLegend(2.546395e-313,2.758595e-313,3.182994e-313,3.853595e-86,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextFont(62);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(19);
   leg->SetFillStyle(0);
   TLegendEntry *entry=leg->AddEntry("h31","CMS Preliminary","");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("h31","p+Pb #sqrt{s_{NN}}= 5.02 TeV","");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("h31","15<p_{T}^{B}<20 GeV/c","");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("h31","Data","pl");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("f31","Fit","l");
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
   entry=leg->AddEntry("background31","Combinatorial Background","l");
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
   
   leg = new TLegend(-6.751917e-93,-3.377998e-192,9.111623e+30,-1.300946e-56,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextFont(62);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(19);
   leg->SetFillStyle(0);
   entry=leg->AddEntry("h31","B meson","");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("h31","M_{B}=5276.91 #pm 2.25 MeV/c^{2}","");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("h31","N_{B}=132 #pm 13","");
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
   c31->Modified();
   c31->cd();
   c31->SetSelected(c31);
}
