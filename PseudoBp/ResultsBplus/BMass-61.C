{
//=========Macro generated from canvas: c61/
//=========  (Fri Apr 25 19:20:51 2014) by ROOT version5.32/00
   TCanvas *c61 = new TCanvas("c61", "",0,0,600,600);
   gStyle->SetOptFit(1);
   gStyle->SetOptStat(0);
   c61->Range(0,0,1,1);
   c61->SetFillColor(0);
   c61->SetBorderMode(0);
   c61->SetBorderSize(2);
   c61->SetTickx(1);
   c61->SetTicky(1);
   c61->SetLeftMargin(0.13);
   c61->SetRightMargin(0.05);
   c61->SetTopMargin(0.05);
   c61->SetBottomMargin(0.13);
   c61->SetFrameFillStyle(0);
   c61->SetFrameBorderMode(0);
   
   TH1D *h61 = new TH1D("h61","",50,5,6);
   h61->SetBinContent(2,1);
   h61->SetBinContent(3,1);
   h61->SetBinContent(4,1);
   h61->SetBinContent(5,1);
   h61->SetBinContent(6,2);
   h61->SetBinContent(8,3);
   h61->SetBinContent(9,1);
   h61->SetBinContent(11,5);
   h61->SetBinContent(12,4);
   h61->SetBinContent(13,7);
   h61->SetBinContent(14,25);
   h61->SetBinContent(15,14);
   h61->SetBinContent(16,4);
   h61->SetBinContent(17,1);
   h61->SetBinContent(18,1);
   h61->SetBinContent(19,1);
   h61->SetBinContent(20,1);
   h61->SetBinContent(22,2);
   h61->SetBinContent(25,2);
   h61->SetBinContent(27,2);
   h61->SetBinContent(28,1);
   h61->SetBinContent(31,1);
   h61->SetBinContent(33,3);
   h61->SetBinContent(34,1);
   h61->SetBinContent(35,1);
   h61->SetBinContent(41,1);
   h61->SetBinContent(42,2);
   h61->SetBinContent(43,1);
   h61->SetBinContent(45,2);
   h61->SetBinContent(47,1);
   h61->SetBinContent(49,2);
   h61->SetBinError(1,1);
   h61->SetBinError(2,1);
   h61->SetBinError(3,1);
   h61->SetBinError(4,1);
   h61->SetBinError(5,1);
   h61->SetBinError(6,1.414214);
   h61->SetBinError(7,1);
   h61->SetBinError(8,1.732051);
   h61->SetBinError(9,1);
   h61->SetBinError(10,1);
   h61->SetBinError(11,2.236068);
   h61->SetBinError(12,2);
   h61->SetBinError(13,2.645751);
   h61->SetBinError(14,5);
   h61->SetBinError(15,3.741657);
   h61->SetBinError(16,2);
   h61->SetBinError(17,1);
   h61->SetBinError(18,1);
   h61->SetBinError(19,1);
   h61->SetBinError(20,1);
   h61->SetBinError(21,1);
   h61->SetBinError(22,1.414214);
   h61->SetBinError(23,1);
   h61->SetBinError(24,1);
   h61->SetBinError(25,1.414214);
   h61->SetBinError(26,1);
   h61->SetBinError(27,1.414214);
   h61->SetBinError(28,1);
   h61->SetBinError(29,1);
   h61->SetBinError(30,1);
   h61->SetBinError(31,1);
   h61->SetBinError(32,1);
   h61->SetBinError(33,1.732051);
   h61->SetBinError(34,1);
   h61->SetBinError(35,1);
   h61->SetBinError(36,1);
   h61->SetBinError(37,1);
   h61->SetBinError(38,1);
   h61->SetBinError(39,1);
   h61->SetBinError(40,1);
   h61->SetBinError(41,1);
   h61->SetBinError(42,1.414214);
   h61->SetBinError(43,1);
   h61->SetBinError(44,1);
   h61->SetBinError(45,1.414214);
   h61->SetBinError(46,1);
   h61->SetBinError(47,1);
   h61->SetBinError(48,1);
   h61->SetBinError(49,1.414214);
   h61->SetBinError(50,1);
   h61->SetMinimum(0);
   h61->SetMaximum(30);
   h61->SetEntries(95);
   h61->SetStats(0);
   
   TF1 *f61 = new TF1("f61","[0]*([7]*Gaus(x,[1],[2])/(sqrt(2*3.14159)*[2])+(1-[7])*Gaus(x,[1],[8])/(sqrt(2*3.14159)*[8]))+[3]+[4]*x+[5]*(4.18604e+01*Gaus(x,4.97611e+00,8.81611e-02)/(sqrt(2*3.14159)*8.81611e-02)+6.72820e+00*Gaus(x,5.10752e+00,2.66663e-02)/(sqrt(2*3.14159)*2.66663e-02)+1.93889e+00*Gaus(x,5.33740e+00,3.52905e-02)/(sqrt(2*3.14159)*3.52905e-02))",5,6);
   f61->SetFillColor(19);
   f61->SetFillStyle(0);
   f61->SetMarkerStyle(20);
   f61->SetLineColor(2);
   f61->SetLineWidth(1);
   f61->SetChisquare(27.94271);
   f61->SetNDF(44);
   f61->GetXaxis()->SetLabelFont(42);
   f61->GetXaxis()->SetLabelOffset(0.007);
   f61->GetXaxis()->SetLabelSize(0.05);
   f61->GetXaxis()->SetTitleSize(0.06);
   f61->GetXaxis()->SetTitleOffset(0.9);
   f61->GetXaxis()->SetTitleFont(42);
   f61->GetYaxis()->SetLabelFont(42);
   f61->GetYaxis()->SetLabelOffset(0.007);
   f61->GetYaxis()->SetLabelSize(0.05);
   f61->GetYaxis()->SetTitleSize(0.06);
   f61->GetYaxis()->SetTitleOffset(1.05);
   f61->GetYaxis()->SetTitleFont(42);
   f61->SetParameter(0,1.069141);
   f61->SetParError(0,0.1549212);
   f61->SetParLimits(0,0,0);
   f61->SetParameter(1,5.272907);
   f61->SetParError(1,0.003323236);
   f61->SetParLimits(1,0,0);
   f61->SetParameter(2,0.03850397);
   f61->SetParError(2,0);
   f61->SetParLimits(2,0.03850397,0.03850397);
   f61->SetParameter(3,2.901405);
   f61->SetParError(3,0.1355893);
   f61->SetParLimits(3,0,0);
   f61->SetParameter(4,-0.3764747);
   f61->SetParError(4,0.02428815);
   f61->SetParLimits(4,-1000,0);
   f61->SetParameter(5,3.189338e-09);
   f61->SetParError(5,0.01712313);
   f61->SetParLimits(5,0,1000);
   f61->SetParameter(6,0.04242641);
   f61->SetParError(6,0.1018234);
   f61->SetParLimits(6,0,0);
   f61->SetParameter(7,0.3230165);
   f61->SetParError(7,0);
   f61->SetParLimits(7,0.3230165,0.3230165);
   f61->SetParameter(8,0.01722579);
   f61->SetParError(8,0);
   f61->SetParLimits(8,0.01722579,0.01722579);
   h61->GetListOfFunctions()->Add(f61);
   h61->SetLineStyle(0);
   h61->SetMarkerStyle(24);
   h61->SetMarkerSize(0.8);
   h61->GetXaxis()->SetTitle("M_{B} (GeV/c^{2})");
   h61->GetXaxis()->CenterTitle(true);
   h61->GetXaxis()->SetLabelFont(42);
   h61->GetXaxis()->SetLabelOffset(0.007);
   h61->GetXaxis()->SetLabelSize(0.05);
   h61->GetXaxis()->SetTitleSize(0.06);
   h61->GetXaxis()->SetTitleOffset(0.9);
   h61->GetXaxis()->SetTitleFont(42);
   h61->GetYaxis()->SetTitle("Entries / (20 MeV/c^{2})");
   h61->GetYaxis()->CenterTitle(true);
   h61->GetYaxis()->SetLabelFont(42);
   h61->GetYaxis()->SetLabelOffset(0.007);
   h61->GetYaxis()->SetLabelSize(0.05);
   h61->GetYaxis()->SetTitleSize(0.06);
   h61->GetYaxis()->SetTitleFont(42);
   h61->GetZaxis()->SetLabelFont(42);
   h61->GetZaxis()->SetLabelOffset(0.007);
   h61->GetZaxis()->SetLabelSize(0.05);
   h61->GetZaxis()->SetTitleSize(0.06);
   h61->GetZaxis()->SetTitleFont(42);
   h61->Draw("e");
   
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
   fBkpi->SetParameter(0,3.189338e-09);
   fBkpi->SetParError(0,0);
   fBkpi->SetParLimits(0,0,0);
   fBkpi->Draw("same");
   
   TF1 *background61 = new TF1("background61","[0]+[1]*x",5,6);
   background61->SetFillColor(19);
   background61->SetFillStyle(0);
   background61->SetMarkerStyle(20);
   background61->SetLineColor(4);
   background61->SetLineWidth(1);
   background61->SetLineStyle(2);
   background61->GetXaxis()->SetLabelFont(42);
   background61->GetXaxis()->SetLabelOffset(0.007);
   background61->GetXaxis()->SetLabelSize(0.05);
   background61->GetXaxis()->SetTitleSize(0.06);
   background61->GetXaxis()->SetTitleOffset(0.9);
   background61->GetXaxis()->SetTitleFont(42);
   background61->GetYaxis()->SetLabelFont(42);
   background61->GetYaxis()->SetLabelOffset(0.007);
   background61->GetYaxis()->SetLabelSize(0.05);
   background61->GetYaxis()->SetTitleSize(0.06);
   background61->GetYaxis()->SetTitleOffset(1.05);
   background61->GetYaxis()->SetTitleFont(42);
   background61->SetParameter(0,2.901405);
   background61->SetParError(0,0);
   background61->SetParLimits(0,0,0);
   background61->SetParameter(1,-0.3764747);
   background61->SetParError(1,0);
   background61->SetParLimits(1,0,0);
   background61->Draw("same");
   
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
   fmass->SetParameter(0,1.069141);
   fmass->SetParError(0,0.1549212);
   fmass->SetParLimits(0,0,0);
   fmass->SetParameter(1,5.272907);
   fmass->SetParError(1,0.003323236);
   fmass->SetParLimits(1,0,0);
   fmass->SetParameter(2,0.03850397);
   fmass->SetParError(2,0);
   fmass->SetParLimits(2,0,0);
   fmass->SetParameter(3,0.3230165);
   fmass->SetParError(3,0);
   fmass->SetParLimits(3,0,0);
   fmass->SetParameter(4,0.01722579);
   fmass->SetParError(4,0);
   fmass->SetParLimits(4,0,0);
   fmass->Draw("same");
   
   TF1 *f61 = new TF1("f61","[0]*([7]*Gaus(x,[1],[2])/(sqrt(2*3.14159)*[2])+(1-[7])*Gaus(x,[1],[8])/(sqrt(2*3.14159)*[8]))+[3]+[4]*x+[5]*(4.18604e+01*Gaus(x,4.97611e+00,8.81611e-02)/(sqrt(2*3.14159)*8.81611e-02)+6.72820e+00*Gaus(x,5.10752e+00,2.66663e-02)/(sqrt(2*3.14159)*2.66663e-02)+1.93889e+00*Gaus(x,5.33740e+00,3.52905e-02)/(sqrt(2*3.14159)*3.52905e-02))",0,1);
   f61->SetFillColor(19);
   f61->SetFillStyle(0);
   f61->SetMarkerStyle(20);
   f61->SetLineColor(2);
   f61->SetLineWidth(1);
   f61->SetChisquare(27.94271);
   f61->SetNDF(44);
   f61->GetXaxis()->SetLabelFont(42);
   f61->GetXaxis()->SetLabelOffset(0.007);
   f61->GetXaxis()->SetLabelSize(0.05);
   f61->GetXaxis()->SetTitleSize(0.06);
   f61->GetXaxis()->SetTitleOffset(0.9);
   f61->GetXaxis()->SetTitleFont(42);
   f61->GetYaxis()->SetLabelFont(42);
   f61->GetYaxis()->SetLabelOffset(0.007);
   f61->GetYaxis()->SetLabelSize(0.05);
   f61->GetYaxis()->SetTitleSize(0.06);
   f61->GetYaxis()->SetTitleOffset(1.05);
   f61->GetYaxis()->SetTitleFont(42);
   f61->SetParameter(0,1.069141);
   f61->SetParError(0,0.1549212);
   f61->SetParLimits(0,0,0);
   f61->SetParameter(1,5.272907);
   f61->SetParError(1,0.003323236);
   f61->SetParLimits(1,0,0);
   f61->SetParameter(2,0.03850397);
   f61->SetParError(2,0);
   f61->SetParLimits(2,0.03850397,0.03850397);
   f61->SetParameter(3,2.901405);
   f61->SetParError(3,0.1355893);
   f61->SetParLimits(3,0,0);
   f61->SetParameter(4,-0.3764747);
   f61->SetParError(4,0.02428815);
   f61->SetParLimits(4,-1000,0);
   f61->SetParameter(5,3.189338e-09);
   f61->SetParError(5,0.01712313);
   f61->SetParLimits(5,0,1000);
   f61->SetParameter(6,0.04242641);
   f61->SetParError(6,0.1018234);
   f61->SetParLimits(6,0,0);
   f61->SetParameter(7,0.3230165);
   f61->SetParError(7,0);
   f61->SetParLimits(7,0.3230165,0.3230165);
   f61->SetParameter(8,0.01722579);
   f61->SetParError(8,0);
   f61->SetParLimits(8,0.01722579,0.01722579);
   f61->Draw("same");
   
   TLegend *leg = new TLegend(2.121996e-314,0,0,8.487983e-314,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextFont(62);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(19);
   leg->SetFillStyle(0);
   TLegendEntry *entry=leg->AddEntry("h61","CMS Preliminary","");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("h61","p+Pb #sqrt{s_{NN}}= 5.02 TeV","");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("h61","20<p_{T}^{B}<25 GeV/c","");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("h61","Data","pl");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("f61","Fit","l");
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
   entry=leg->AddEntry("background61","Combinatorial Background","l");
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
   
   leg = new TLegend(2.3705e-310,4.08147e-33,3.398581e-315,2.3705e-310,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextFont(62);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(19);
   leg->SetFillStyle(0);
   entry=leg->AddEntry("h61","B meson","");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("h61","M_{B}=5272.91 #pm 3.32 MeV/c^{2}","");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("h61","N_{B}=53 #pm 8","");
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
   c61->Modified();
   c61->cd();
   c61->SetSelected(c61);
}
