{
//=========Macro generated from canvas: c3/
//=========  (Sun Apr 27 17:54:32 2014) by ROOT version5.32/00
   TCanvas *c3 = new TCanvas("c3", "",0,0,600,600);
   gStyle->SetOptFit(1);
   gStyle->SetOptStat(0);
   c3->Range(0,0,1,1);
   c3->SetFillColor(0);
   c3->SetBorderMode(0);
   c3->SetBorderSize(2);
   c3->SetTickx(1);
   c3->SetTicky(1);
   c3->SetLeftMargin(0.13);
   c3->SetRightMargin(0.05);
   c3->SetTopMargin(0.05);
   c3->SetBottomMargin(0.13);
   c3->SetFrameFillStyle(0);
   c3->SetFrameBorderMode(0);
   
   TH1D *h3 = new TH1D("h3","",50,5,6);
   h3->SetBinContent(1,14);
   h3->SetBinContent(2,14);
   h3->SetBinContent(3,10);
   h3->SetBinContent(4,9);
   h3->SetBinContent(5,14);
   h3->SetBinContent(6,11);
   h3->SetBinContent(7,12);
   h3->SetBinContent(8,3);
   h3->SetBinContent(9,8);
   h3->SetBinContent(10,5);
   h3->SetBinContent(11,6);
   h3->SetBinContent(12,18);
   h3->SetBinContent(13,20);
   h3->SetBinContent(14,47);
   h3->SetBinContent(15,63);
   h3->SetBinContent(16,17);
   h3->SetBinContent(17,18);
   h3->SetBinContent(18,5);
   h3->SetBinContent(19,6);
   h3->SetBinContent(20,1);
   h3->SetBinContent(21,3);
   h3->SetBinContent(22,4);
   h3->SetBinContent(23,3);
   h3->SetBinContent(24,5);
   h3->SetBinContent(25,3);
   h3->SetBinContent(26,11);
   h3->SetBinContent(27,6);
   h3->SetBinContent(28,3);
   h3->SetBinContent(29,6);
   h3->SetBinContent(30,4);
   h3->SetBinContent(31,3);
   h3->SetBinContent(32,6);
   h3->SetBinContent(33,1);
   h3->SetBinContent(34,2);
   h3->SetBinContent(35,2);
   h3->SetBinContent(36,4);
   h3->SetBinContent(37,4);
   h3->SetBinContent(39,1);
   h3->SetBinContent(40,7);
   h3->SetBinContent(41,3);
   h3->SetBinContent(42,5);
   h3->SetBinContent(43,7);
   h3->SetBinContent(44,5);
   h3->SetBinContent(45,1);
   h3->SetBinContent(46,6);
   h3->SetBinContent(47,7);
   h3->SetBinContent(48,2);
   h3->SetBinContent(49,6);
   h3->SetBinContent(50,4);
   h3->SetBinError(1,3.741657);
   h3->SetBinError(2,3.741657);
   h3->SetBinError(3,3.162278);
   h3->SetBinError(4,3);
   h3->SetBinError(5,3.741657);
   h3->SetBinError(6,3.316625);
   h3->SetBinError(7,3.464102);
   h3->SetBinError(8,1.732051);
   h3->SetBinError(9,2.828427);
   h3->SetBinError(10,2.236068);
   h3->SetBinError(11,2.44949);
   h3->SetBinError(12,4.242641);
   h3->SetBinError(13,4.472136);
   h3->SetBinError(14,6.855655);
   h3->SetBinError(15,7.937254);
   h3->SetBinError(16,4.123106);
   h3->SetBinError(17,4.242641);
   h3->SetBinError(18,2.236068);
   h3->SetBinError(19,2.44949);
   h3->SetBinError(20,1);
   h3->SetBinError(21,1.732051);
   h3->SetBinError(22,2);
   h3->SetBinError(23,1.732051);
   h3->SetBinError(24,2.236068);
   h3->SetBinError(25,1.732051);
   h3->SetBinError(26,3.316625);
   h3->SetBinError(27,2.44949);
   h3->SetBinError(28,1.732051);
   h3->SetBinError(29,2.44949);
   h3->SetBinError(30,2);
   h3->SetBinError(31,1.732051);
   h3->SetBinError(32,2.44949);
   h3->SetBinError(33,1);
   h3->SetBinError(34,1.414214);
   h3->SetBinError(35,1.414214);
   h3->SetBinError(36,2);
   h3->SetBinError(37,2);
   h3->SetBinError(38,1);
   h3->SetBinError(39,1);
   h3->SetBinError(40,2.645751);
   h3->SetBinError(41,1.732051);
   h3->SetBinError(42,2.236068);
   h3->SetBinError(43,2.645751);
   h3->SetBinError(44,2.236068);
   h3->SetBinError(45,1);
   h3->SetBinError(46,2.44949);
   h3->SetBinError(47,2.645751);
   h3->SetBinError(48,1.414214);
   h3->SetBinError(49,2.44949);
   h3->SetBinError(50,2);
   h3->SetMinimum(0);
   h3->SetMaximum(75.6);
   h3->SetEntries(425);
   h3->SetStats(0);
   
   TF1 *f3 = new TF1("f3","[0]*([7]*Gaus(x,[1],[2])/(sqrt(2*3.14159)*[2])+(1-[7])*Gaus(x,[1],[8])/(sqrt(2*3.14159)*[8]))+[3]+[4]*x+[5]*(4.18604e+01*Gaus(x,4.97611e+00,8.81611e-02)/(sqrt(2*3.14159)*8.81611e-02)+6.72820e+00*Gaus(x,5.10752e+00,2.66663e-02)/(sqrt(2*3.14159)*2.66663e-02)+1.93889e+00*Gaus(x,5.33740e+00,3.52905e-02)/(sqrt(2*3.14159)*3.52905e-02))",5,6);
   f3->SetFillColor(19);
   f3->SetFillStyle(0);
   f3->SetMarkerStyle(20);
   f3->SetLineColor(2);
   f3->SetLineWidth(1);
   f3->SetChisquare(122.2302);
   f3->SetNDF(44);
   f3->GetXaxis()->SetLabelFont(42);
   f3->GetXaxis()->SetLabelOffset(0.007);
   f3->GetXaxis()->SetLabelSize(0.05);
   f3->GetXaxis()->SetTitleSize(0.06);
   f3->GetXaxis()->SetTitleOffset(0.9);
   f3->GetXaxis()->SetTitleFont(42);
   f3->GetYaxis()->SetLabelFont(42);
   f3->GetYaxis()->SetLabelOffset(0.007);
   f3->GetYaxis()->SetLabelSize(0.05);
   f3->GetYaxis()->SetTitleSize(0.06);
   f3->GetYaxis()->SetTitleOffset(1.05);
   f3->GetYaxis()->SetTitleFont(42);
   f3->SetParameter(0,2.620801);
   f3->SetParError(0,0.2543975);
   f3->SetParLimits(0,0,0);
   f3->SetParameter(1,5.281129);
   f3->SetParError(1,0.001946285);
   f3->SetParLimits(1,0,0);
   f3->SetParameter(2,0.01573478);
   f3->SetParError(2,0);
   f3->SetParLimits(2,0.01573478,0.01573478);
   f3->SetParameter(3,20.05992);
   f3->SetParError(3,0.333841);
   f3->SetParLimits(3,0,0);
   f3->SetParameter(4,-2.751972);
   f3->SetParError(4,0.05943311);
   f3->SetParLimits(4,-1000,0);
   f3->SetParameter(5,0.03802641);
   f3->SetParError(5,0.008391221);
   f3->SetParLimits(5,0,1000);
   f3->SetParameter(6,0);
   f3->SetParError(6,19.2);
   f3->SetParLimits(6,0,0);
   f3->SetParameter(7,0.9331789);
   f3->SetParError(7,0);
   f3->SetParLimits(7,0.9331789,0.9331789);
   f3->SetParameter(8,0.05);
   f3->SetParError(8,0);
   f3->SetParLimits(8,0.05,0.05);
   h3->GetListOfFunctions()->Add(f3);
   h3->SetLineStyle(0);
   h3->SetMarkerStyle(24);
   h3->SetMarkerSize(0.8);
   h3->GetXaxis()->SetTitle("M_{B} (GeV/c^{2})");
   h3->GetXaxis()->CenterTitle(true);
   h3->GetXaxis()->SetLabelFont(42);
   h3->GetXaxis()->SetLabelOffset(0.007);
   h3->GetXaxis()->SetLabelSize(0.05);
   h3->GetXaxis()->SetTitleSize(0.06);
   h3->GetXaxis()->SetTitleOffset(0.9);
   h3->GetXaxis()->SetTitleFont(42);
   h3->GetYaxis()->SetTitle("Entries / (20 MeV/c^{2})");
   h3->GetYaxis()->CenterTitle(true);
   h3->GetYaxis()->SetLabelFont(42);
   h3->GetYaxis()->SetLabelOffset(0.007);
   h3->GetYaxis()->SetLabelSize(0.05);
   h3->GetYaxis()->SetTitleSize(0.06);
   h3->GetYaxis()->SetTitleFont(42);
   h3->GetZaxis()->SetLabelFont(42);
   h3->GetZaxis()->SetLabelOffset(0.007);
   h3->GetZaxis()->SetLabelSize(0.05);
   h3->GetZaxis()->SetTitleSize(0.06);
   h3->GetZaxis()->SetTitleFont(42);
   h3->Draw("e");
   
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
   fBkpi->SetParameter(0,0.03802641);
   fBkpi->SetParError(0,0);
   fBkpi->SetParLimits(0,0,0);
   fBkpi->Draw("same");
   
   TF1 *background3 = new TF1("background3","[0]+[1]*x",5,6);
   background3->SetFillColor(19);
   background3->SetFillStyle(0);
   background3->SetMarkerStyle(20);
   background3->SetLineColor(4);
   background3->SetLineWidth(1);
   background3->SetLineStyle(2);
   background3->GetXaxis()->SetLabelFont(42);
   background3->GetXaxis()->SetLabelOffset(0.007);
   background3->GetXaxis()->SetLabelSize(0.05);
   background3->GetXaxis()->SetTitleSize(0.06);
   background3->GetXaxis()->SetTitleOffset(0.9);
   background3->GetXaxis()->SetTitleFont(42);
   background3->GetYaxis()->SetLabelFont(42);
   background3->GetYaxis()->SetLabelOffset(0.007);
   background3->GetYaxis()->SetLabelSize(0.05);
   background3->GetYaxis()->SetTitleSize(0.06);
   background3->GetYaxis()->SetTitleOffset(1.05);
   background3->GetYaxis()->SetTitleFont(42);
   background3->SetParameter(0,20.05992);
   background3->SetParError(0,0);
   background3->SetParLimits(0,0,0);
   background3->SetParameter(1,-2.751972);
   background3->SetParError(1,0);
   background3->SetParLimits(1,0,0);
   background3->Draw("same");
   
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
   fmass->SetParameter(0,2.620801);
   fmass->SetParError(0,0.2543975);
   fmass->SetParLimits(0,0,0);
   fmass->SetParameter(1,5.281129);
   fmass->SetParError(1,0.001946285);
   fmass->SetParLimits(1,0,0);
   fmass->SetParameter(2,0.01573478);
   fmass->SetParError(2,0);
   fmass->SetParLimits(2,0,0);
   fmass->SetParameter(3,0.9331789);
   fmass->SetParError(3,0);
   fmass->SetParLimits(3,0,0);
   fmass->SetParameter(4,0.05);
   fmass->SetParError(4,0);
   fmass->SetParLimits(4,0,0);
   fmass->Draw("same");
   
   TF1 *f3 = new TF1("f3","[0]*([7]*Gaus(x,[1],[2])/(sqrt(2*3.14159)*[2])+(1-[7])*Gaus(x,[1],[8])/(sqrt(2*3.14159)*[8]))+[3]+[4]*x+[5]*(4.18604e+01*Gaus(x,4.97611e+00,8.81611e-02)/(sqrt(2*3.14159)*8.81611e-02)+6.72820e+00*Gaus(x,5.10752e+00,2.66663e-02)/(sqrt(2*3.14159)*2.66663e-02)+1.93889e+00*Gaus(x,5.33740e+00,3.52905e-02)/(sqrt(2*3.14159)*3.52905e-02))",0,1);
   f3->SetFillColor(19);
   f3->SetFillStyle(0);
   f3->SetMarkerStyle(20);
   f3->SetLineColor(2);
   f3->SetLineWidth(1);
   f3->SetChisquare(122.2302);
   f3->SetNDF(44);
   f3->GetXaxis()->SetLabelFont(42);
   f3->GetXaxis()->SetLabelOffset(0.007);
   f3->GetXaxis()->SetLabelSize(0.05);
   f3->GetXaxis()->SetTitleSize(0.06);
   f3->GetXaxis()->SetTitleOffset(0.9);
   f3->GetXaxis()->SetTitleFont(42);
   f3->GetYaxis()->SetLabelFont(42);
   f3->GetYaxis()->SetLabelOffset(0.007);
   f3->GetYaxis()->SetLabelSize(0.05);
   f3->GetYaxis()->SetTitleSize(0.06);
   f3->GetYaxis()->SetTitleOffset(1.05);
   f3->GetYaxis()->SetTitleFont(42);
   f3->SetParameter(0,2.620801);
   f3->SetParError(0,0.2543975);
   f3->SetParLimits(0,0,0);
   f3->SetParameter(1,5.281129);
   f3->SetParError(1,0.001946285);
   f3->SetParLimits(1,0,0);
   f3->SetParameter(2,0.01573478);
   f3->SetParError(2,0);
   f3->SetParLimits(2,0.01573478,0.01573478);
   f3->SetParameter(3,20.05992);
   f3->SetParError(3,0.333841);
   f3->SetParLimits(3,0,0);
   f3->SetParameter(4,-2.751972);
   f3->SetParError(4,0.05943311);
   f3->SetParLimits(4,-1000,0);
   f3->SetParameter(5,0.03802641);
   f3->SetParError(5,0.008391221);
   f3->SetParLimits(5,0,1000);
   f3->SetParameter(6,0);
   f3->SetParError(6,19.2);
   f3->SetParLimits(6,0,0);
   f3->SetParameter(7,0.9331789);
   f3->SetParError(7,0);
   f3->SetParLimits(7,0.9331789,0.9331789);
   f3->SetParameter(8,0.05);
   f3->SetParError(8,0);
   f3->SetParLimits(8,0.05,0.05);
   f3->Draw("same");
   
   TLegend *leg = new TLegend(2.353513e-310,0,0,8.611237e-314,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextFont(62);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(19);
   leg->SetFillStyle(0);
   TLegendEntry *entry=leg->AddEntry("h3","CMS Preliminary","");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("h3","p+Pb #sqrt{s_{NN}}= 5.02 TeV","");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("h3","0<y_{CM}^{B}<1","");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("h3","Data","pl");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("f3","Fit","l");
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
   entry=leg->AddEntry("background3","Combinatorial Background","l");
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
   
   leg = new TLegend(0,2.353513e-310,0,0,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextFont(62);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(19);
   leg->SetFillStyle(0);
   entry=leg->AddEntry("h3","B meson","");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("h3","M_{B}=5281.13 #pm 1.95 MeV/c^{2}","");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("h3","N_{B}=131 #pm 13","");
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
   c3->Modified();
   c3->cd();
   c3->SetSelected(c3);
}
