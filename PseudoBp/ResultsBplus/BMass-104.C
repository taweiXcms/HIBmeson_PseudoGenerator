{
//=========Macro generated from canvas: c104/
//=========  (Fri Apr 25 19:35:12 2014) by ROOT version5.32/00
   TCanvas *c104 = new TCanvas("c104", "",0,0,600,600);
   gStyle->SetOptFit(1);
   gStyle->SetOptStat(0);
   c104->Range(0,0,1,1);
   c104->SetFillColor(0);
   c104->SetBorderMode(0);
   c104->SetBorderSize(2);
   c104->SetTickx(1);
   c104->SetTicky(1);
   c104->SetLeftMargin(0.13);
   c104->SetRightMargin(0.05);
   c104->SetTopMargin(0.05);
   c104->SetBottomMargin(0.13);
   c104->SetFrameFillStyle(0);
   c104->SetFrameBorderMode(0);
   
   TH1D *h104 = new TH1D("h104","",50,5,6);
   h104->SetBinContent(9,2);
   h104->SetBinContent(13,6);
   h104->SetBinContent(14,11);
   h104->SetBinContent(15,10);
   h104->SetBinContent(16,4);
   h104->SetBinContent(27,2);
   h104->SetBinContent(30,1);
   h104->SetBinContent(34,1);
   h104->SetBinContent(37,1);
   h104->SetBinContent(40,1);
   h104->SetBinContent(42,1);
   h104->SetBinContent(49,1);
   h104->SetBinError(1,1);
   h104->SetBinError(2,1);
   h104->SetBinError(3,1);
   h104->SetBinError(4,1);
   h104->SetBinError(5,1);
   h104->SetBinError(6,1);
   h104->SetBinError(7,1);
   h104->SetBinError(8,1);
   h104->SetBinError(9,1.414214);
   h104->SetBinError(10,1);
   h104->SetBinError(11,1);
   h104->SetBinError(12,1);
   h104->SetBinError(13,2.44949);
   h104->SetBinError(14,3.316625);
   h104->SetBinError(15,3.162278);
   h104->SetBinError(16,2);
   h104->SetBinError(17,1);
   h104->SetBinError(18,1);
   h104->SetBinError(19,1);
   h104->SetBinError(20,1);
   h104->SetBinError(21,1);
   h104->SetBinError(22,1);
   h104->SetBinError(23,1);
   h104->SetBinError(24,1);
   h104->SetBinError(25,1);
   h104->SetBinError(26,1);
   h104->SetBinError(27,1.414214);
   h104->SetBinError(28,1);
   h104->SetBinError(29,1);
   h104->SetBinError(30,1);
   h104->SetBinError(31,1);
   h104->SetBinError(32,1);
   h104->SetBinError(33,1);
   h104->SetBinError(34,1);
   h104->SetBinError(35,1);
   h104->SetBinError(36,1);
   h104->SetBinError(37,1);
   h104->SetBinError(38,1);
   h104->SetBinError(39,1);
   h104->SetBinError(40,1);
   h104->SetBinError(41,1);
   h104->SetBinError(42,1);
   h104->SetBinError(43,1);
   h104->SetBinError(44,1);
   h104->SetBinError(45,1);
   h104->SetBinError(46,1);
   h104->SetBinError(47,1);
   h104->SetBinError(48,1);
   h104->SetBinError(49,1);
   h104->SetBinError(50,1);
   h104->SetMinimum(0);
   h104->SetMaximum(13.2);
   h104->SetEntries(41);
   h104->SetStats(0);
   
   TF1 *f104 = new TF1("f104","[0]*([7]*Gaus(x,[1],[2])/(sqrt(2*3.14159)*[2])+(1-[7])*Gaus(x,[1],[8])/(sqrt(2*3.14159)*[8]))+[3]+[4]*x+[5]*(4.18604e+01*Gaus(x,4.97611e+00,8.81611e-02)/(sqrt(2*3.14159)*8.81611e-02)+6.72820e+00*Gaus(x,5.10752e+00,2.66663e-02)/(sqrt(2*3.14159)*2.66663e-02)+1.93889e+00*Gaus(x,5.33740e+00,3.52905e-02)/(sqrt(2*3.14159)*3.52905e-02))",5,6);
   f104->SetFillColor(19);
   f104->SetFillStyle(0);
   f104->SetMarkerStyle(20);
   f104->SetLineColor(2);
   f104->SetLineWidth(1);
   f104->SetChisquare(11.94051);
   f104->SetNDF(44);
   f104->GetXaxis()->SetLabelFont(42);
   f104->GetXaxis()->SetLabelOffset(0.007);
   f104->GetXaxis()->SetLabelSize(0.05);
   f104->GetXaxis()->SetTitleSize(0.06);
   f104->GetXaxis()->SetTitleOffset(0.9);
   f104->GetXaxis()->SetTitleFont(42);
   f104->GetYaxis()->SetLabelFont(42);
   f104->GetYaxis()->SetLabelOffset(0.007);
   f104->GetYaxis()->SetLabelSize(0.05);
   f104->GetYaxis()->SetTitleSize(0.06);
   f104->GetYaxis()->SetTitleOffset(1.05);
   f104->GetYaxis()->SetTitleFont(42);
   f104->SetParameter(0,0.611595);
   f104->SetParError(0,0.1129674);
   f104->SetParLimits(0,0,0);
   f104->SetParameter(1,5.277254);
   f104->SetParError(1,0.004058184);
   f104->SetParLimits(1,0,0);
   f104->SetParameter(2,0.05);
   f104->SetParError(2,0);
   f104->SetParLimits(2,0.05,0.05);
   f104->SetParameter(3,0.208353);
   f104->SetParError(3,0.06887499);
   f104->SetParLimits(3,0,0);
   f104->SetParameter(4,-3.777018e-09);
   f104->SetParError(4,0.1086861);
   f104->SetParLimits(4,-1000,0);
   f104->SetParameter(5,4.252154e-11);
   f104->SetParError(5,0.0004750378);
   f104->SetParLimits(5,0,1000);
   f104->SetParameter(6,0);
   f104->SetParError(6,8);
   f104->SetParLimits(6,0,0);
   f104->SetParameter(7,0.1439248);
   f104->SetParError(7,0);
   f104->SetParLimits(7,0.1439248,0.1439248);
   f104->SetParameter(8,0.02011116);
   f104->SetParError(8,0);
   f104->SetParLimits(8,0.02011116,0.02011116);
   h104->GetListOfFunctions()->Add(f104);
   h104->SetLineStyle(0);
   h104->SetMarkerStyle(24);
   h104->SetMarkerSize(0.8);
   h104->GetXaxis()->SetTitle("M_{B} (GeV/c^{2})");
   h104->GetXaxis()->CenterTitle(true);
   h104->GetXaxis()->SetLabelFont(42);
   h104->GetXaxis()->SetLabelOffset(0.007);
   h104->GetXaxis()->SetLabelSize(0.05);
   h104->GetXaxis()->SetTitleSize(0.06);
   h104->GetXaxis()->SetTitleOffset(0.9);
   h104->GetXaxis()->SetTitleFont(42);
   h104->GetYaxis()->SetTitle("Entries / (20 MeV/c^{2})");
   h104->GetYaxis()->CenterTitle(true);
   h104->GetYaxis()->SetLabelFont(42);
   h104->GetYaxis()->SetLabelOffset(0.007);
   h104->GetYaxis()->SetLabelSize(0.05);
   h104->GetYaxis()->SetTitleSize(0.06);
   h104->GetYaxis()->SetTitleFont(42);
   h104->GetZaxis()->SetLabelFont(42);
   h104->GetZaxis()->SetLabelOffset(0.007);
   h104->GetZaxis()->SetLabelSize(0.05);
   h104->GetZaxis()->SetTitleSize(0.06);
   h104->GetZaxis()->SetTitleFont(42);
   h104->Draw("e");
   
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
   fBkpi->SetParameter(0,4.252154e-11);
   fBkpi->SetParError(0,0);
   fBkpi->SetParLimits(0,0,0);
   fBkpi->Draw("same");
   
   TF1 *background104 = new TF1("background104","[0]+[1]*x",5,6);
   background104->SetFillColor(19);
   background104->SetFillStyle(0);
   background104->SetMarkerStyle(20);
   background104->SetLineColor(4);
   background104->SetLineWidth(1);
   background104->SetLineStyle(2);
   background104->GetXaxis()->SetLabelFont(42);
   background104->GetXaxis()->SetLabelOffset(0.007);
   background104->GetXaxis()->SetLabelSize(0.05);
   background104->GetXaxis()->SetTitleSize(0.06);
   background104->GetXaxis()->SetTitleOffset(0.9);
   background104->GetXaxis()->SetTitleFont(42);
   background104->GetYaxis()->SetLabelFont(42);
   background104->GetYaxis()->SetLabelOffset(0.007);
   background104->GetYaxis()->SetLabelSize(0.05);
   background104->GetYaxis()->SetTitleSize(0.06);
   background104->GetYaxis()->SetTitleOffset(1.05);
   background104->GetYaxis()->SetTitleFont(42);
   background104->SetParameter(0,0.208353);
   background104->SetParError(0,0);
   background104->SetParLimits(0,0,0);
   background104->SetParameter(1,-3.777018e-09);
   background104->SetParError(1,0);
   background104->SetParLimits(1,0,0);
   background104->Draw("same");
   
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
   fmass->SetParameter(0,0.611595);
   fmass->SetParError(0,0.1129674);
   fmass->SetParLimits(0,0,0);
   fmass->SetParameter(1,5.277254);
   fmass->SetParError(1,0.004058184);
   fmass->SetParLimits(1,0,0);
   fmass->SetParameter(2,0.05);
   fmass->SetParError(2,0);
   fmass->SetParLimits(2,0,0);
   fmass->SetParameter(3,0.1439248);
   fmass->SetParError(3,0);
   fmass->SetParLimits(3,0,0);
   fmass->SetParameter(4,0.02011116);
   fmass->SetParError(4,0);
   fmass->SetParLimits(4,0,0);
   fmass->Draw("same");
   
   TF1 *f104 = new TF1("f104","[0]*([7]*Gaus(x,[1],[2])/(sqrt(2*3.14159)*[2])+(1-[7])*Gaus(x,[1],[8])/(sqrt(2*3.14159)*[8]))+[3]+[4]*x+[5]*(4.18604e+01*Gaus(x,4.97611e+00,8.81611e-02)/(sqrt(2*3.14159)*8.81611e-02)+6.72820e+00*Gaus(x,5.10752e+00,2.66663e-02)/(sqrt(2*3.14159)*2.66663e-02)+1.93889e+00*Gaus(x,5.33740e+00,3.52905e-02)/(sqrt(2*3.14159)*3.52905e-02))",0,1);
   f104->SetFillColor(19);
   f104->SetFillStyle(0);
   f104->SetMarkerStyle(20);
   f104->SetLineColor(2);
   f104->SetLineWidth(1);
   f104->SetChisquare(11.94051);
   f104->SetNDF(44);
   f104->GetXaxis()->SetLabelFont(42);
   f104->GetXaxis()->SetLabelOffset(0.007);
   f104->GetXaxis()->SetLabelSize(0.05);
   f104->GetXaxis()->SetTitleSize(0.06);
   f104->GetXaxis()->SetTitleOffset(0.9);
   f104->GetXaxis()->SetTitleFont(42);
   f104->GetYaxis()->SetLabelFont(42);
   f104->GetYaxis()->SetLabelOffset(0.007);
   f104->GetYaxis()->SetLabelSize(0.05);
   f104->GetYaxis()->SetTitleSize(0.06);
   f104->GetYaxis()->SetTitleOffset(1.05);
   f104->GetYaxis()->SetTitleFont(42);
   f104->SetParameter(0,0.611595);
   f104->SetParError(0,0.1129674);
   f104->SetParLimits(0,0,0);
   f104->SetParameter(1,5.277254);
   f104->SetParError(1,0.004058184);
   f104->SetParLimits(1,0,0);
   f104->SetParameter(2,0.05);
   f104->SetParError(2,0);
   f104->SetParLimits(2,0.05,0.05);
   f104->SetParameter(3,0.208353);
   f104->SetParError(3,0.06887499);
   f104->SetParLimits(3,0,0);
   f104->SetParameter(4,-3.777018e-09);
   f104->SetParError(4,0.1086861);
   f104->SetParLimits(4,-1000,0);
   f104->SetParameter(5,4.252154e-11);
   f104->SetParError(5,0.0004750378);
   f104->SetParLimits(5,0,1000);
   f104->SetParameter(6,0);
   f104->SetParError(6,8);
   f104->SetParLimits(6,0,0);
   f104->SetParameter(7,0.1439248);
   f104->SetParError(7,0);
   f104->SetParLimits(7,0.1439248,0.1439248);
   f104->SetParameter(8,0.02011116);
   f104->SetParError(8,0);
   f104->SetParLimits(8,0.02011116,0.02011116);
   f104->Draw("same");
   
   TLegend *leg = new TLegend(2.546395e-313,2.758595e-313,3.182994e-313,8.487983e-314,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextFont(62);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(19);
   leg->SetFillStyle(0);
   TLegendEntry *entry=leg->AddEntry("h104","CMS Preliminary","");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("h104","p+Pb #sqrt{s_{NN}}= 5.02 TeV","");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("h104","30<p_{T}^{B}<60 GeV/c","");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("h104","Data","pl");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("f104","Fit","l");
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
   entry=leg->AddEntry("background104","Combinatorial Background","l");
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
   
   leg = new TLegend(0,0,0,0,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextFont(62);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(19);
   leg->SetFillStyle(0);
   entry=leg->AddEntry("h104","B meson","");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("h104","M_{B}=5277.25 #pm 4.06 MeV/c^{2}","");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("h104","N_{B}=31 #pm 6","");
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
   c104->Modified();
   c104->cd();
   c104->SetSelected(c104);
}
