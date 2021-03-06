{
//=========Macro generated from canvas: c62/
//=========  (Fri Apr 25 19:21:07 2014) by ROOT version5.32/00
   TCanvas *c62 = new TCanvas("c62", "",0,0,600,600);
   gStyle->SetOptFit(1);
   gStyle->SetOptStat(0);
   c62->Range(0,0,1,1);
   c62->SetFillColor(0);
   c62->SetBorderMode(0);
   c62->SetBorderSize(2);
   c62->SetTickx(1);
   c62->SetTicky(1);
   c62->SetLeftMargin(0.13);
   c62->SetRightMargin(0.05);
   c62->SetTopMargin(0.05);
   c62->SetBottomMargin(0.13);
   c62->SetFrameFillStyle(0);
   c62->SetFrameBorderMode(0);
   
   TH1D *h62 = new TH1D("h62","",50,5,6);
   h62->SetBinContent(1,1);
   h62->SetBinContent(4,1);
   h62->SetBinContent(5,2);
   h62->SetBinContent(8,4);
   h62->SetBinContent(10,1);
   h62->SetBinContent(11,2);
   h62->SetBinContent(13,3);
   h62->SetBinContent(14,23);
   h62->SetBinContent(15,14);
   h62->SetBinContent(16,9);
   h62->SetBinContent(17,3);
   h62->SetBinContent(18,2);
   h62->SetBinContent(20,2);
   h62->SetBinContent(21,4);
   h62->SetBinContent(22,2);
   h62->SetBinContent(24,2);
   h62->SetBinContent(25,2);
   h62->SetBinContent(27,1);
   h62->SetBinContent(28,1);
   h62->SetBinContent(29,1);
   h62->SetBinContent(30,2);
   h62->SetBinContent(32,1);
   h62->SetBinContent(35,2);
   h62->SetBinContent(40,1);
   h62->SetBinContent(41,1);
   h62->SetBinContent(42,2);
   h62->SetBinContent(43,1);
   h62->SetBinContent(45,1);
   h62->SetBinError(1,1);
   h62->SetBinError(2,1);
   h62->SetBinError(3,1);
   h62->SetBinError(4,1);
   h62->SetBinError(5,1.414214);
   h62->SetBinError(6,1);
   h62->SetBinError(7,1);
   h62->SetBinError(8,2);
   h62->SetBinError(9,1);
   h62->SetBinError(10,1);
   h62->SetBinError(11,1.414214);
   h62->SetBinError(12,1);
   h62->SetBinError(13,1.732051);
   h62->SetBinError(14,4.795832);
   h62->SetBinError(15,3.741657);
   h62->SetBinError(16,3);
   h62->SetBinError(17,1.732051);
   h62->SetBinError(18,1.414214);
   h62->SetBinError(19,1);
   h62->SetBinError(20,1.414214);
   h62->SetBinError(21,2);
   h62->SetBinError(22,1.414214);
   h62->SetBinError(23,1);
   h62->SetBinError(24,1.414214);
   h62->SetBinError(25,1.414214);
   h62->SetBinError(26,1);
   h62->SetBinError(27,1);
   h62->SetBinError(28,1);
   h62->SetBinError(29,1);
   h62->SetBinError(30,1.414214);
   h62->SetBinError(31,1);
   h62->SetBinError(32,1);
   h62->SetBinError(33,1);
   h62->SetBinError(34,1);
   h62->SetBinError(35,1.414214);
   h62->SetBinError(36,1);
   h62->SetBinError(37,1);
   h62->SetBinError(38,1);
   h62->SetBinError(39,1);
   h62->SetBinError(40,1);
   h62->SetBinError(41,1);
   h62->SetBinError(42,1.414214);
   h62->SetBinError(43,1);
   h62->SetBinError(44,1);
   h62->SetBinError(45,1);
   h62->SetBinError(46,1);
   h62->SetBinError(47,1);
   h62->SetBinError(48,1);
   h62->SetBinError(49,1);
   h62->SetBinError(50,1);
   h62->SetMinimum(0);
   h62->SetMaximum(27.6);
   h62->SetEntries(91);
   h62->SetStats(0);
   
   TF1 *f62 = new TF1("f62","[0]*([7]*Gaus(x,[1],[2])/(sqrt(2*3.14159)*[2])+(1-[7])*Gaus(x,[1],[8])/(sqrt(2*3.14159)*[8]))+[3]+[4]*x+[5]*(4.18604e+01*Gaus(x,4.97611e+00,8.81611e-02)/(sqrt(2*3.14159)*8.81611e-02)+6.72820e+00*Gaus(x,5.10752e+00,2.66663e-02)/(sqrt(2*3.14159)*2.66663e-02)+1.93889e+00*Gaus(x,5.33740e+00,3.52905e-02)/(sqrt(2*3.14159)*3.52905e-02))",5,6);
   f62->SetFillColor(19);
   f62->SetFillStyle(0);
   f62->SetMarkerStyle(20);
   f62->SetLineColor(2);
   f62->SetLineWidth(1);
   f62->SetChisquare(39.18342);
   f62->SetNDF(44);
   f62->GetXaxis()->SetLabelFont(42);
   f62->GetXaxis()->SetLabelOffset(0.007);
   f62->GetXaxis()->SetLabelSize(0.05);
   f62->GetXaxis()->SetTitleSize(0.06);
   f62->GetXaxis()->SetTitleOffset(0.9);
   f62->GetXaxis()->SetTitleFont(42);
   f62->GetYaxis()->SetLabelFont(42);
   f62->GetYaxis()->SetLabelOffset(0.007);
   f62->GetYaxis()->SetLabelSize(0.05);
   f62->GetYaxis()->SetTitleSize(0.06);
   f62->GetYaxis()->SetTitleOffset(1.05);
   f62->GetYaxis()->SetTitleFont(42);
   f62->SetParameter(0,0.9731622);
   f62->SetParError(0,0.1483875);
   f62->SetParLimits(0,0,0);
   f62->SetParameter(1,5.282383);
   f62->SetParError(1,0.003413806);
   f62->SetParLimits(1,0,0);
   f62->SetParameter(2,0.03850397);
   f62->SetParError(2,0);
   f62->SetParLimits(2,0.03850397,0.03850397);
   f62->SetParameter(3,6.346529);
   f62->SetParError(3,0.129387);
   f62->SetParLimits(3,0,0);
   f62->SetParameter(4,-0.9999478);
   f62->SetParError(4,0.02308689);
   f62->SetParLimits(4,-1000,0);
   f62->SetParameter(5,2.403633e-11);
   f62->SetParError(5,0.00112315);
   f62->SetParLimits(5,0,1000);
   f62->SetParameter(6,0.04242641);
   f62->SetParError(6,0.1018234);
   f62->SetParLimits(6,0,0);
   f62->SetParameter(7,0.3230165);
   f62->SetParError(7,0);
   f62->SetParLimits(7,0.3230165,0.3230165);
   f62->SetParameter(8,0.01722579);
   f62->SetParError(8,0);
   f62->SetParLimits(8,0.01722579,0.01722579);
   h62->GetListOfFunctions()->Add(f62);
   h62->SetLineStyle(0);
   h62->SetMarkerStyle(24);
   h62->SetMarkerSize(0.8);
   h62->GetXaxis()->SetTitle("M_{B} (GeV/c^{2})");
   h62->GetXaxis()->CenterTitle(true);
   h62->GetXaxis()->SetLabelFont(42);
   h62->GetXaxis()->SetLabelOffset(0.007);
   h62->GetXaxis()->SetLabelSize(0.05);
   h62->GetXaxis()->SetTitleSize(0.06);
   h62->GetXaxis()->SetTitleOffset(0.9);
   h62->GetXaxis()->SetTitleFont(42);
   h62->GetYaxis()->SetTitle("Entries / (20 MeV/c^{2})");
   h62->GetYaxis()->CenterTitle(true);
   h62->GetYaxis()->SetLabelFont(42);
   h62->GetYaxis()->SetLabelOffset(0.007);
   h62->GetYaxis()->SetLabelSize(0.05);
   h62->GetYaxis()->SetTitleSize(0.06);
   h62->GetYaxis()->SetTitleFont(42);
   h62->GetZaxis()->SetLabelFont(42);
   h62->GetZaxis()->SetLabelOffset(0.007);
   h62->GetZaxis()->SetLabelSize(0.05);
   h62->GetZaxis()->SetTitleSize(0.06);
   h62->GetZaxis()->SetTitleFont(42);
   h62->Draw("e");
   
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
   fBkpi->SetParameter(0,2.403633e-11);
   fBkpi->SetParError(0,0);
   fBkpi->SetParLimits(0,0,0);
   fBkpi->Draw("same");
   
   TF1 *background62 = new TF1("background62","[0]+[1]*x",5,6);
   background62->SetFillColor(19);
   background62->SetFillStyle(0);
   background62->SetMarkerStyle(20);
   background62->SetLineColor(4);
   background62->SetLineWidth(1);
   background62->SetLineStyle(2);
   background62->GetXaxis()->SetLabelFont(42);
   background62->GetXaxis()->SetLabelOffset(0.007);
   background62->GetXaxis()->SetLabelSize(0.05);
   background62->GetXaxis()->SetTitleSize(0.06);
   background62->GetXaxis()->SetTitleOffset(0.9);
   background62->GetXaxis()->SetTitleFont(42);
   background62->GetYaxis()->SetLabelFont(42);
   background62->GetYaxis()->SetLabelOffset(0.007);
   background62->GetYaxis()->SetLabelSize(0.05);
   background62->GetYaxis()->SetTitleSize(0.06);
   background62->GetYaxis()->SetTitleOffset(1.05);
   background62->GetYaxis()->SetTitleFont(42);
   background62->SetParameter(0,6.346529);
   background62->SetParError(0,0);
   background62->SetParLimits(0,0,0);
   background62->SetParameter(1,-0.9999478);
   background62->SetParError(1,0);
   background62->SetParLimits(1,0,0);
   background62->Draw("same");
   
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
   fmass->SetParameter(0,0.9731622);
   fmass->SetParError(0,0.1483875);
   fmass->SetParLimits(0,0,0);
   fmass->SetParameter(1,5.282383);
   fmass->SetParError(1,0.003413806);
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
   
   TF1 *f62 = new TF1("f62","[0]*([7]*Gaus(x,[1],[2])/(sqrt(2*3.14159)*[2])+(1-[7])*Gaus(x,[1],[8])/(sqrt(2*3.14159)*[8]))+[3]+[4]*x+[5]*(4.18604e+01*Gaus(x,4.97611e+00,8.81611e-02)/(sqrt(2*3.14159)*8.81611e-02)+6.72820e+00*Gaus(x,5.10752e+00,2.66663e-02)/(sqrt(2*3.14159)*2.66663e-02)+1.93889e+00*Gaus(x,5.33740e+00,3.52905e-02)/(sqrt(2*3.14159)*3.52905e-02))",0,1);
   f62->SetFillColor(19);
   f62->SetFillStyle(0);
   f62->SetMarkerStyle(20);
   f62->SetLineColor(2);
   f62->SetLineWidth(1);
   f62->SetChisquare(39.18342);
   f62->SetNDF(44);
   f62->GetXaxis()->SetLabelFont(42);
   f62->GetXaxis()->SetLabelOffset(0.007);
   f62->GetXaxis()->SetLabelSize(0.05);
   f62->GetXaxis()->SetTitleSize(0.06);
   f62->GetXaxis()->SetTitleOffset(0.9);
   f62->GetXaxis()->SetTitleFont(42);
   f62->GetYaxis()->SetLabelFont(42);
   f62->GetYaxis()->SetLabelOffset(0.007);
   f62->GetYaxis()->SetLabelSize(0.05);
   f62->GetYaxis()->SetTitleSize(0.06);
   f62->GetYaxis()->SetTitleOffset(1.05);
   f62->GetYaxis()->SetTitleFont(42);
   f62->SetParameter(0,0.9731622);
   f62->SetParError(0,0.1483875);
   f62->SetParLimits(0,0,0);
   f62->SetParameter(1,5.282383);
   f62->SetParError(1,0.003413806);
   f62->SetParLimits(1,0,0);
   f62->SetParameter(2,0.03850397);
   f62->SetParError(2,0);
   f62->SetParLimits(2,0.03850397,0.03850397);
   f62->SetParameter(3,6.346529);
   f62->SetParError(3,0.129387);
   f62->SetParLimits(3,0,0);
   f62->SetParameter(4,-0.9999478);
   f62->SetParError(4,0.02308689);
   f62->SetParLimits(4,-1000,0);
   f62->SetParameter(5,2.403633e-11);
   f62->SetParError(5,0.00112315);
   f62->SetParLimits(5,0,1000);
   f62->SetParameter(6,0.04242641);
   f62->SetParError(6,0.1018234);
   f62->SetParLimits(6,0,0);
   f62->SetParameter(7,0.3230165);
   f62->SetParError(7,0);
   f62->SetParLimits(7,0.3230165,0.3230165);
   f62->SetParameter(8,0.01722579);
   f62->SetParError(8,0);
   f62->SetParLimits(8,0.01722579,0.01722579);
   f62->Draw("same");
   
   TLegend *leg = new TLegend(2.837108e-311,2.837108e-311,2.845596e-311,2.845596e-311,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextFont(62);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(19);
   leg->SetFillStyle(0);
   TLegendEntry *entry=leg->AddEntry("h62","CMS Preliminary","");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("h62","p+Pb #sqrt{s_{NN}}= 5.02 TeV","");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("h62","20<p_{T}^{B}<25 GeV/c","");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("h62","Data","pl");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("f62","Fit","l");
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
   entry=leg->AddEntry("background62","Combinatorial Background","l");
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
   
   leg = new TLegend(2.637641e-311,2.680081e-311,2.680081e-311,2.680081e-311,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextFont(62);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(19);
   leg->SetFillStyle(0);
   entry=leg->AddEntry("h62","B meson","");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("h62","M_{B}=5282.38 #pm 3.41 MeV/c^{2}","");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("h62","N_{B}=49 #pm 7","");
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
   c62->Modified();
   c62->cd();
   c62->SetSelected(c62);
}
