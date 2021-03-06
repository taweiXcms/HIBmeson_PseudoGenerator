{
//=========Macro generated from canvas: c26/
//=========  (Fri Apr 25 19:09:10 2014) by ROOT version5.32/00
   TCanvas *c26 = new TCanvas("c26", "",0,0,600,600);
   gStyle->SetOptFit(1);
   gStyle->SetOptStat(0);
   c26->Range(0,0,1,1);
   c26->SetFillColor(0);
   c26->SetBorderMode(0);
   c26->SetBorderSize(2);
   c26->SetTickx(1);
   c26->SetTicky(1);
   c26->SetLeftMargin(0.13);
   c26->SetRightMargin(0.05);
   c26->SetTopMargin(0.05);
   c26->SetBottomMargin(0.13);
   c26->SetFrameFillStyle(0);
   c26->SetFrameBorderMode(0);
   
   TH1D *h26 = new TH1D("h26","",50,5,6);
   h26->SetBinContent(1,5);
   h26->SetBinContent(2,3);
   h26->SetBinContent(3,2);
   h26->SetBinContent(4,5);
   h26->SetBinContent(5,6);
   h26->SetBinContent(6,2);
   h26->SetBinContent(7,3);
   h26->SetBinContent(8,3);
   h26->SetBinContent(9,5);
   h26->SetBinContent(10,3);
   h26->SetBinContent(11,7);
   h26->SetBinContent(12,14);
   h26->SetBinContent(13,14);
   h26->SetBinContent(14,57);
   h26->SetBinContent(15,59);
   h26->SetBinContent(16,18);
   h26->SetBinContent(17,8);
   h26->SetBinContent(18,6);
   h26->SetBinContent(19,9);
   h26->SetBinContent(20,3);
   h26->SetBinContent(21,8);
   h26->SetBinContent(22,2);
   h26->SetBinContent(23,6);
   h26->SetBinContent(24,6);
   h26->SetBinContent(25,2);
   h26->SetBinContent(26,3);
   h26->SetBinContent(27,1);
   h26->SetBinContent(28,2);
   h26->SetBinContent(29,2);
   h26->SetBinContent(30,6);
   h26->SetBinContent(31,2);
   h26->SetBinContent(32,1);
   h26->SetBinContent(33,3);
   h26->SetBinContent(34,3);
   h26->SetBinContent(35,1);
   h26->SetBinContent(36,2);
   h26->SetBinContent(37,1);
   h26->SetBinContent(38,2);
   h26->SetBinContent(39,2);
   h26->SetBinContent(40,2);
   h26->SetBinContent(41,1);
   h26->SetBinContent(42,1);
   h26->SetBinContent(43,3);
   h26->SetBinContent(44,3);
   h26->SetBinContent(45,2);
   h26->SetBinContent(46,1);
   h26->SetBinContent(49,2);
   h26->SetBinContent(50,4);
   h26->SetBinError(1,2.236068);
   h26->SetBinError(2,1.732051);
   h26->SetBinError(3,1.414214);
   h26->SetBinError(4,2.236068);
   h26->SetBinError(5,2.44949);
   h26->SetBinError(6,1.414214);
   h26->SetBinError(7,1.732051);
   h26->SetBinError(8,1.732051);
   h26->SetBinError(9,2.236068);
   h26->SetBinError(10,1.732051);
   h26->SetBinError(11,2.645751);
   h26->SetBinError(12,3.741657);
   h26->SetBinError(13,3.741657);
   h26->SetBinError(14,7.549834);
   h26->SetBinError(15,7.681146);
   h26->SetBinError(16,4.242641);
   h26->SetBinError(17,2.828427);
   h26->SetBinError(18,2.44949);
   h26->SetBinError(19,3);
   h26->SetBinError(20,1.732051);
   h26->SetBinError(21,2.828427);
   h26->SetBinError(22,1.414214);
   h26->SetBinError(23,2.44949);
   h26->SetBinError(24,2.44949);
   h26->SetBinError(25,1.414214);
   h26->SetBinError(26,1.732051);
   h26->SetBinError(27,1);
   h26->SetBinError(28,1.414214);
   h26->SetBinError(29,1.414214);
   h26->SetBinError(30,2.44949);
   h26->SetBinError(31,1.414214);
   h26->SetBinError(32,1);
   h26->SetBinError(33,1.732051);
   h26->SetBinError(34,1.732051);
   h26->SetBinError(35,1);
   h26->SetBinError(36,1.414214);
   h26->SetBinError(37,1);
   h26->SetBinError(38,1.414214);
   h26->SetBinError(39,1.414214);
   h26->SetBinError(40,1.414214);
   h26->SetBinError(41,1);
   h26->SetBinError(42,1);
   h26->SetBinError(43,1.732051);
   h26->SetBinError(44,1.732051);
   h26->SetBinError(45,1.414214);
   h26->SetBinError(46,1);
   h26->SetBinError(47,1);
   h26->SetBinError(48,1);
   h26->SetBinError(49,1.414214);
   h26->SetBinError(50,2);
   h26->SetMinimum(0);
   h26->SetMaximum(70.8);
   h26->SetEntries(306);
   h26->SetStats(0);
   
   TF1 *f26 = new TF1("f26","[0]*([7]*Gaus(x,[1],[2])/(sqrt(2*3.14159)*[2])+(1-[7])*Gaus(x,[1],[8])/(sqrt(2*3.14159)*[8]))+[3]+[4]*x+[5]*(4.18604e+01*Gaus(x,4.97611e+00,8.81611e-02)/(sqrt(2*3.14159)*8.81611e-02)+6.72820e+00*Gaus(x,5.10752e+00,2.66663e-02)/(sqrt(2*3.14159)*2.66663e-02)+1.93889e+00*Gaus(x,5.33740e+00,3.52905e-02)/(sqrt(2*3.14159)*3.52905e-02))",5,6);
   f26->SetFillColor(19);
   f26->SetFillStyle(0);
   f26->SetMarkerStyle(20);
   f26->SetLineColor(2);
   f26->SetLineWidth(1);
   f26->SetChisquare(47.53492);
   f26->SetNDF(44);
   f26->GetXaxis()->SetLabelFont(42);
   f26->GetXaxis()->SetLabelOffset(0.007);
   f26->GetXaxis()->SetLabelSize(0.05);
   f26->GetXaxis()->SetTitleSize(0.06);
   f26->GetXaxis()->SetTitleOffset(0.9);
   f26->GetXaxis()->SetTitleFont(42);
   f26->GetYaxis()->SetLabelFont(42);
   f26->GetYaxis()->SetLabelOffset(0.007);
   f26->GetYaxis()->SetLabelSize(0.05);
   f26->GetYaxis()->SetTitleSize(0.06);
   f26->GetYaxis()->SetTitleOffset(1.05);
   f26->GetYaxis()->SetTitleFont(42);
   f26->SetParameter(0,3.094941);
   f26->SetParError(0,0.2692973);
   f26->SetParLimits(0,0,0);
   f26->SetParameter(1,5.280797);
   f26->SetParError(1,0.001950012);
   f26->SetParLimits(1,0,0);
   f26->SetParameter(2,0.03959756);
   f26->SetParError(2,0);
   f26->SetParLimits(2,0.03959756,0.03959756);
   f26->SetParameter(3,20.56235);
   f26->SetParError(3,0.2442843);
   f26->SetParLimits(3,0,0);
   f26->SetParameter(4,-3.188607);
   f26->SetParError(4,0.04335772);
   f26->SetParLimits(4,-1000,0);
   f26->SetParameter(5,7.127077e-10);
   f26->SetParError(5,0.003331329);
   f26->SetParLimits(5,0,1000);
   f26->SetParameter(6,0);
   f26->SetParError(6,2.4);
   f26->SetParLimits(6,0,0);
   f26->SetParameter(7,0.331348);
   f26->SetParError(7,0);
   f26->SetParLimits(7,0.331348,0.331348);
   f26->SetParameter(8,0.01638322);
   f26->SetParError(8,0);
   f26->SetParLimits(8,0.01638322,0.01638322);
   h26->GetListOfFunctions()->Add(f26);
   h26->SetLineStyle(0);
   h26->SetMarkerStyle(24);
   h26->SetMarkerSize(0.8);
   h26->GetXaxis()->SetTitle("M_{B} (GeV/c^{2})");
   h26->GetXaxis()->CenterTitle(true);
   h26->GetXaxis()->SetLabelFont(42);
   h26->GetXaxis()->SetLabelOffset(0.007);
   h26->GetXaxis()->SetLabelSize(0.05);
   h26->GetXaxis()->SetTitleSize(0.06);
   h26->GetXaxis()->SetTitleOffset(0.9);
   h26->GetXaxis()->SetTitleFont(42);
   h26->GetYaxis()->SetTitle("Entries / (20 MeV/c^{2})");
   h26->GetYaxis()->CenterTitle(true);
   h26->GetYaxis()->SetLabelFont(42);
   h26->GetYaxis()->SetLabelOffset(0.007);
   h26->GetYaxis()->SetLabelSize(0.05);
   h26->GetYaxis()->SetTitleSize(0.06);
   h26->GetYaxis()->SetTitleFont(42);
   h26->GetZaxis()->SetLabelFont(42);
   h26->GetZaxis()->SetLabelOffset(0.007);
   h26->GetZaxis()->SetLabelSize(0.05);
   h26->GetZaxis()->SetTitleSize(0.06);
   h26->GetZaxis()->SetTitleFont(42);
   h26->Draw("e");
   
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
   fBkpi->SetParameter(0,7.127077e-10);
   fBkpi->SetParError(0,0);
   fBkpi->SetParLimits(0,0,0);
   fBkpi->Draw("same");
   
   TF1 *background26 = new TF1("background26","[0]+[1]*x",5,6);
   background26->SetFillColor(19);
   background26->SetFillStyle(0);
   background26->SetMarkerStyle(20);
   background26->SetLineColor(4);
   background26->SetLineWidth(1);
   background26->SetLineStyle(2);
   background26->GetXaxis()->SetLabelFont(42);
   background26->GetXaxis()->SetLabelOffset(0.007);
   background26->GetXaxis()->SetLabelSize(0.05);
   background26->GetXaxis()->SetTitleSize(0.06);
   background26->GetXaxis()->SetTitleOffset(0.9);
   background26->GetXaxis()->SetTitleFont(42);
   background26->GetYaxis()->SetLabelFont(42);
   background26->GetYaxis()->SetLabelOffset(0.007);
   background26->GetYaxis()->SetLabelSize(0.05);
   background26->GetYaxis()->SetTitleSize(0.06);
   background26->GetYaxis()->SetTitleOffset(1.05);
   background26->GetYaxis()->SetTitleFont(42);
   background26->SetParameter(0,20.56235);
   background26->SetParError(0,0);
   background26->SetParLimits(0,0,0);
   background26->SetParameter(1,-3.188607);
   background26->SetParError(1,0);
   background26->SetParLimits(1,0,0);
   background26->Draw("same");
   
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
   fmass->SetParameter(0,3.094941);
   fmass->SetParError(0,0.2692973);
   fmass->SetParLimits(0,0,0);
   fmass->SetParameter(1,5.280797);
   fmass->SetParError(1,0.001950012);
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
   
   TF1 *f26 = new TF1("f26","[0]*([7]*Gaus(x,[1],[2])/(sqrt(2*3.14159)*[2])+(1-[7])*Gaus(x,[1],[8])/(sqrt(2*3.14159)*[8]))+[3]+[4]*x+[5]*(4.18604e+01*Gaus(x,4.97611e+00,8.81611e-02)/(sqrt(2*3.14159)*8.81611e-02)+6.72820e+00*Gaus(x,5.10752e+00,2.66663e-02)/(sqrt(2*3.14159)*2.66663e-02)+1.93889e+00*Gaus(x,5.33740e+00,3.52905e-02)/(sqrt(2*3.14159)*3.52905e-02))",0,1);
   f26->SetFillColor(19);
   f26->SetFillStyle(0);
   f26->SetMarkerStyle(20);
   f26->SetLineColor(2);
   f26->SetLineWidth(1);
   f26->SetChisquare(47.53492);
   f26->SetNDF(44);
   f26->GetXaxis()->SetLabelFont(42);
   f26->GetXaxis()->SetLabelOffset(0.007);
   f26->GetXaxis()->SetLabelSize(0.05);
   f26->GetXaxis()->SetTitleSize(0.06);
   f26->GetXaxis()->SetTitleOffset(0.9);
   f26->GetXaxis()->SetTitleFont(42);
   f26->GetYaxis()->SetLabelFont(42);
   f26->GetYaxis()->SetLabelOffset(0.007);
   f26->GetYaxis()->SetLabelSize(0.05);
   f26->GetYaxis()->SetTitleSize(0.06);
   f26->GetYaxis()->SetTitleOffset(1.05);
   f26->GetYaxis()->SetTitleFont(42);
   f26->SetParameter(0,3.094941);
   f26->SetParError(0,0.2692973);
   f26->SetParLimits(0,0,0);
   f26->SetParameter(1,5.280797);
   f26->SetParError(1,0.001950012);
   f26->SetParLimits(1,0,0);
   f26->SetParameter(2,0.03959756);
   f26->SetParError(2,0);
   f26->SetParLimits(2,0.03959756,0.03959756);
   f26->SetParameter(3,20.56235);
   f26->SetParError(3,0.2442843);
   f26->SetParLimits(3,0,0);
   f26->SetParameter(4,-3.188607);
   f26->SetParError(4,0.04335772);
   f26->SetParLimits(4,-1000,0);
   f26->SetParameter(5,7.127077e-10);
   f26->SetParError(5,0.003331329);
   f26->SetParLimits(5,0,1000);
   f26->SetParameter(6,0);
   f26->SetParError(6,2.4);
   f26->SetParLimits(6,0,0);
   f26->SetParameter(7,0.331348);
   f26->SetParError(7,0);
   f26->SetParLimits(7,0.331348,0.331348);
   f26->SetParameter(8,0.01638322);
   f26->SetParError(8,0);
   f26->SetParLimits(8,0.01638322,0.01638322);
   f26->Draw("same");
   
   TLegend *leg = new TLegend(2.121996e-314,0,0,0,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextFont(62);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(19);
   leg->SetFillStyle(0);
   TLegendEntry *entry=leg->AddEntry("h26","CMS Preliminary","");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("h26","p+Pb #sqrt{s_{NN}}= 5.02 TeV","");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("h26","15<p_{T}^{B}<20 GeV/c","");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("h26","Data","pl");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("f26","Fit","l");
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
   entry=leg->AddEntry("background26","Combinatorial Background","l");
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
   
   leg = new TLegend(2.824376e-311,2.824376e-311,2.824376e-311,2.824376e-311,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextFont(62);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(19);
   leg->SetFillStyle(0);
   entry=leg->AddEntry("h26","B meson","");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("h26","M_{B}=5280.80 #pm 1.95 MeV/c^{2}","");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("h26","N_{B}=155 #pm 13","");
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
   c26->Modified();
   c26->cd();
   c26->SetSelected(c26);
}
