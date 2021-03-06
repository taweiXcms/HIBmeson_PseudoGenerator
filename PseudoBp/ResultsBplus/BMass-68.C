{
//=========Macro generated from canvas: c68/
//=========  (Fri Apr 25 19:23:04 2014) by ROOT version5.32/00
   TCanvas *c68 = new TCanvas("c68", "",0,0,600,600);
   gStyle->SetOptFit(1);
   gStyle->SetOptStat(0);
   c68->Range(0,0,1,1);
   c68->SetFillColor(0);
   c68->SetBorderMode(0);
   c68->SetBorderSize(2);
   c68->SetTickx(1);
   c68->SetTicky(1);
   c68->SetLeftMargin(0.13);
   c68->SetRightMargin(0.05);
   c68->SetTopMargin(0.05);
   c68->SetBottomMargin(0.13);
   c68->SetFrameFillStyle(0);
   c68->SetFrameBorderMode(0);
   
   TH1D *h68 = new TH1D("h68","",50,5,6);
   h68->SetBinContent(4,1);
   h68->SetBinContent(8,1);
   h68->SetBinContent(9,1);
   h68->SetBinContent(10,1);
   h68->SetBinContent(11,1);
   h68->SetBinContent(13,4);
   h68->SetBinContent(14,4);
   h68->SetBinContent(15,5);
   h68->SetBinContent(18,1);
   h68->SetBinContent(19,2);
   h68->SetBinContent(20,1);
   h68->SetBinContent(24,1);
   h68->SetBinContent(26,1);
   h68->SetBinContent(29,1);
   h68->SetBinContent(31,1);
   h68->SetBinContent(49,1);
   h68->SetBinError(1,1);
   h68->SetBinError(2,1);
   h68->SetBinError(3,1);
   h68->SetBinError(4,1);
   h68->SetBinError(5,1);
   h68->SetBinError(6,1);
   h68->SetBinError(7,1);
   h68->SetBinError(8,1);
   h68->SetBinError(9,1);
   h68->SetBinError(10,1);
   h68->SetBinError(11,1);
   h68->SetBinError(12,1);
   h68->SetBinError(13,2);
   h68->SetBinError(14,2);
   h68->SetBinError(15,2.236068);
   h68->SetBinError(16,1);
   h68->SetBinError(17,1);
   h68->SetBinError(18,1);
   h68->SetBinError(19,1.414214);
   h68->SetBinError(20,1);
   h68->SetBinError(21,1);
   h68->SetBinError(22,1);
   h68->SetBinError(23,1);
   h68->SetBinError(24,1);
   h68->SetBinError(25,1);
   h68->SetBinError(26,1);
   h68->SetBinError(27,1);
   h68->SetBinError(28,1);
   h68->SetBinError(29,1);
   h68->SetBinError(30,1);
   h68->SetBinError(31,1);
   h68->SetBinError(32,1);
   h68->SetBinError(33,1);
   h68->SetBinError(34,1);
   h68->SetBinError(35,1);
   h68->SetBinError(36,1);
   h68->SetBinError(37,1);
   h68->SetBinError(38,1);
   h68->SetBinError(39,1);
   h68->SetBinError(40,1);
   h68->SetBinError(41,1);
   h68->SetBinError(42,1);
   h68->SetBinError(43,1);
   h68->SetBinError(44,1);
   h68->SetBinError(45,1);
   h68->SetBinError(46,1);
   h68->SetBinError(47,1);
   h68->SetBinError(48,1);
   h68->SetBinError(49,1);
   h68->SetBinError(50,1);
   h68->SetMinimum(0);
   h68->SetMaximum(6);
   h68->SetEntries(27);
   h68->SetStats(0);
   
   TF1 *f68 = new TF1("f68","[0]*([7]*Gaus(x,[1],[2])/(sqrt(2*3.14159)*[2])+(1-[7])*Gaus(x,[1],[8])/(sqrt(2*3.14159)*[8]))+[3]+[4]*x+[5]*(4.18604e+01*Gaus(x,4.97611e+00,8.81611e-02)/(sqrt(2*3.14159)*8.81611e-02)+6.72820e+00*Gaus(x,5.10752e+00,2.66663e-02)/(sqrt(2*3.14159)*2.66663e-02)+1.93889e+00*Gaus(x,5.33740e+00,3.52905e-02)/(sqrt(2*3.14159)*3.52905e-02))",5,6);
   f68->SetFillColor(19);
   f68->SetFillStyle(0);
   f68->SetMarkerStyle(20);
   f68->SetLineColor(2);
   f68->SetLineWidth(1);
   f68->SetChisquare(13.05473);
   f68->SetNDF(44);
   f68->GetXaxis()->SetLabelFont(42);
   f68->GetXaxis()->SetLabelOffset(0.007);
   f68->GetXaxis()->SetLabelSize(0.05);
   f68->GetXaxis()->SetTitleSize(0.06);
   f68->GetXaxis()->SetTitleOffset(0.9);
   f68->GetXaxis()->SetTitleFont(42);
   f68->GetYaxis()->SetLabelFont(42);
   f68->GetYaxis()->SetLabelOffset(0.007);
   f68->GetYaxis()->SetLabelSize(0.05);
   f68->GetYaxis()->SetTitleSize(0.06);
   f68->GetYaxis()->SetTitleOffset(1.05);
   f68->GetYaxis()->SetTitleFont(42);
   f68->SetParameter(0,0.2536244);
   f68->SetParError(0,0.07773569);
   f68->SetParLimits(0,0,0);
   f68->SetParameter(1,5.271461);
   f68->SetParError(1,0.006960089);
   f68->SetParLimits(1,0,0);
   f68->SetParameter(2,0.04043018);
   f68->SetParError(2,0);
   f68->SetParLimits(2,0.04043018,0.04043018);
   f68->SetParameter(3,2.546597);
   f68->SetParError(3,0.06720193);
   f68->SetParLimits(3,0,0);
   f68->SetParameter(4,-0.4109588);
   f68->SetParError(4,0.01176099);
   f68->SetParLimits(4,-1000,0);
   f68->SetParameter(5,1.876277e-10);
   f68->SetParError(5,0.0008772055);
   f68->SetParLimits(5,0,1000);
   f68->SetParameter(6,0);
   f68->SetParError(6,2.4);
   f68->SetParLimits(6,0,0);
   f68->SetParameter(7,0.2650151);
   f68->SetParError(7,0);
   f68->SetParLimits(7,0.2650151,0.2650151);
   f68->SetParameter(8,0.01908153);
   f68->SetParError(8,0);
   f68->SetParLimits(8,0.01908153,0.01908153);
   h68->GetListOfFunctions()->Add(f68);
   h68->SetLineStyle(0);
   h68->SetMarkerStyle(24);
   h68->SetMarkerSize(0.8);
   h68->GetXaxis()->SetTitle("M_{B} (GeV/c^{2})");
   h68->GetXaxis()->CenterTitle(true);
   h68->GetXaxis()->SetLabelFont(42);
   h68->GetXaxis()->SetLabelOffset(0.007);
   h68->GetXaxis()->SetLabelSize(0.05);
   h68->GetXaxis()->SetTitleSize(0.06);
   h68->GetXaxis()->SetTitleOffset(0.9);
   h68->GetXaxis()->SetTitleFont(42);
   h68->GetYaxis()->SetTitle("Entries / (20 MeV/c^{2})");
   h68->GetYaxis()->CenterTitle(true);
   h68->GetYaxis()->SetLabelFont(42);
   h68->GetYaxis()->SetLabelOffset(0.007);
   h68->GetYaxis()->SetLabelSize(0.05);
   h68->GetYaxis()->SetTitleSize(0.06);
   h68->GetYaxis()->SetTitleFont(42);
   h68->GetZaxis()->SetLabelFont(42);
   h68->GetZaxis()->SetLabelOffset(0.007);
   h68->GetZaxis()->SetLabelSize(0.05);
   h68->GetZaxis()->SetTitleSize(0.06);
   h68->GetZaxis()->SetTitleFont(42);
   h68->Draw("e");
   
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
   fBkpi->SetParameter(0,1.876277e-10);
   fBkpi->SetParError(0,0);
   fBkpi->SetParLimits(0,0,0);
   fBkpi->Draw("same");
   
   TF1 *background68 = new TF1("background68","[0]+[1]*x",5,6);
   background68->SetFillColor(19);
   background68->SetFillStyle(0);
   background68->SetMarkerStyle(20);
   background68->SetLineColor(4);
   background68->SetLineWidth(1);
   background68->SetLineStyle(2);
   background68->GetXaxis()->SetLabelFont(42);
   background68->GetXaxis()->SetLabelOffset(0.007);
   background68->GetXaxis()->SetLabelSize(0.05);
   background68->GetXaxis()->SetTitleSize(0.06);
   background68->GetXaxis()->SetTitleOffset(0.9);
   background68->GetXaxis()->SetTitleFont(42);
   background68->GetYaxis()->SetLabelFont(42);
   background68->GetYaxis()->SetLabelOffset(0.007);
   background68->GetYaxis()->SetLabelSize(0.05);
   background68->GetYaxis()->SetTitleSize(0.06);
   background68->GetYaxis()->SetTitleOffset(1.05);
   background68->GetYaxis()->SetTitleFont(42);
   background68->SetParameter(0,2.546597);
   background68->SetParError(0,0);
   background68->SetParLimits(0,0,0);
   background68->SetParameter(1,-0.4109588);
   background68->SetParError(1,0);
   background68->SetParLimits(1,0,0);
   background68->Draw("same");
   
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
   fmass->SetParameter(0,0.2536244);
   fmass->SetParError(0,0.07773569);
   fmass->SetParLimits(0,0,0);
   fmass->SetParameter(1,5.271461);
   fmass->SetParError(1,0.006960089);
   fmass->SetParLimits(1,0,0);
   fmass->SetParameter(2,0.04043018);
   fmass->SetParError(2,0);
   fmass->SetParLimits(2,0,0);
   fmass->SetParameter(3,0.2650151);
   fmass->SetParError(3,0);
   fmass->SetParLimits(3,0,0);
   fmass->SetParameter(4,0.01908153);
   fmass->SetParError(4,0);
   fmass->SetParLimits(4,0,0);
   fmass->Draw("same");
   
   TF1 *f68 = new TF1("f68","[0]*([7]*Gaus(x,[1],[2])/(sqrt(2*3.14159)*[2])+(1-[7])*Gaus(x,[1],[8])/(sqrt(2*3.14159)*[8]))+[3]+[4]*x+[5]*(4.18604e+01*Gaus(x,4.97611e+00,8.81611e-02)/(sqrt(2*3.14159)*8.81611e-02)+6.72820e+00*Gaus(x,5.10752e+00,2.66663e-02)/(sqrt(2*3.14159)*2.66663e-02)+1.93889e+00*Gaus(x,5.33740e+00,3.52905e-02)/(sqrt(2*3.14159)*3.52905e-02))",0,1);
   f68->SetFillColor(19);
   f68->SetFillStyle(0);
   f68->SetMarkerStyle(20);
   f68->SetLineColor(2);
   f68->SetLineWidth(1);
   f68->SetChisquare(13.05473);
   f68->SetNDF(44);
   f68->GetXaxis()->SetLabelFont(42);
   f68->GetXaxis()->SetLabelOffset(0.007);
   f68->GetXaxis()->SetLabelSize(0.05);
   f68->GetXaxis()->SetTitleSize(0.06);
   f68->GetXaxis()->SetTitleOffset(0.9);
   f68->GetXaxis()->SetTitleFont(42);
   f68->GetYaxis()->SetLabelFont(42);
   f68->GetYaxis()->SetLabelOffset(0.007);
   f68->GetYaxis()->SetLabelSize(0.05);
   f68->GetYaxis()->SetTitleSize(0.06);
   f68->GetYaxis()->SetTitleOffset(1.05);
   f68->GetYaxis()->SetTitleFont(42);
   f68->SetParameter(0,0.2536244);
   f68->SetParError(0,0.07773569);
   f68->SetParLimits(0,0,0);
   f68->SetParameter(1,5.271461);
   f68->SetParError(1,0.006960089);
   f68->SetParLimits(1,0,0);
   f68->SetParameter(2,0.04043018);
   f68->SetParError(2,0);
   f68->SetParLimits(2,0.04043018,0.04043018);
   f68->SetParameter(3,2.546597);
   f68->SetParError(3,0.06720193);
   f68->SetParLimits(3,0,0);
   f68->SetParameter(4,-0.4109588);
   f68->SetParError(4,0.01176099);
   f68->SetParLimits(4,-1000,0);
   f68->SetParameter(5,1.876277e-10);
   f68->SetParError(5,0.0008772055);
   f68->SetParLimits(5,0,1000);
   f68->SetParameter(6,0);
   f68->SetParError(6,2.4);
   f68->SetParLimits(6,0,0);
   f68->SetParameter(7,0.2650151);
   f68->SetParError(7,0);
   f68->SetParLimits(7,0.2650151,0.2650151);
   f68->SetParameter(8,0.01908153);
   f68->SetParError(8,0);
   f68->SetParLimits(8,0.01908153,0.01908153);
   f68->Draw("same");
   
   TLegend *leg = new TLegend(2.121996e-314,0,0,0,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextFont(62);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(19);
   leg->SetFillStyle(0);
   TLegendEntry *entry=leg->AddEntry("h68","CMS Preliminary","");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("h68","p+Pb #sqrt{s_{NN}}= 5.02 TeV","");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("h68","25<p_{T}^{B}<30 GeV/c","");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("h68","Data","pl");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("f68","Fit","l");
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
   entry=leg->AddEntry("background68","Combinatorial Background","l");
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
   
   leg = new TLegend(2.215364e-311,2.215364e-311,2.215364e-311,2.215364e-311,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextFont(62);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(19);
   leg->SetFillStyle(0);
   entry=leg->AddEntry("h68","B meson","");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("h68","M_{B}=5271.46 #pm 6.96 MeV/c^{2}","");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("h68","N_{B}=13 #pm 4","");
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
   c68->Modified();
   c68->cd();
   c68->SetSelected(c68);
}
