{
//=========Macro generated from canvas: c4/
//=========  (Mon Apr  7 19:53:43 2014) by ROOT version5.32/00
   TCanvas *c4 = new TCanvas("c4", "",0,0,600,600);
   gStyle->SetOptFit(1);
   gStyle->SetOptStat(0);
   c4->Range(0,0,1,1);
   c4->SetFillColor(0);
   c4->SetBorderMode(0);
   c4->SetBorderSize(2);
   c4->SetTickx(1);
   c4->SetTicky(1);
   c4->SetLeftMargin(0.13);
   c4->SetRightMargin(0.05);
   c4->SetTopMargin(0.05);
   c4->SetBottomMargin(0.13);
   c4->SetFrameFillStyle(0);
   c4->SetFrameBorderMode(0);
   
   TH1D *h4 = new TH1D("h4","",50,5,6);
   h4->SetBinContent(1,66);
   h4->SetBinContent(2,66);
   h4->SetBinContent(3,63);
   h4->SetBinContent(4,76);
   h4->SetBinContent(5,76);
   h4->SetBinContent(6,61);
   h4->SetBinContent(7,50);
   h4->SetBinContent(8,49);
   h4->SetBinContent(9,57);
   h4->SetBinContent(10,60);
   h4->SetBinContent(11,59);
   h4->SetBinContent(12,65);
   h4->SetBinContent(13,80);
   h4->SetBinContent(14,123);
   h4->SetBinContent(15,115);
   h4->SetBinContent(16,71);
   h4->SetBinContent(17,56);
   h4->SetBinContent(18,41);
   h4->SetBinContent(19,36);
   h4->SetBinContent(20,39);
   h4->SetBinContent(21,32);
   h4->SetBinContent(22,27);
   h4->SetBinContent(23,33);
   h4->SetBinContent(24,35);
   h4->SetBinContent(25,30);
   h4->SetBinContent(26,20);
   h4->SetBinContent(27,25);
   h4->SetBinContent(28,24);
   h4->SetBinContent(29,26);
   h4->SetBinContent(30,20);
   h4->SetBinContent(31,32);
   h4->SetBinContent(32,17);
   h4->SetBinContent(33,22);
   h4->SetBinContent(34,27);
   h4->SetBinContent(35,12);
   h4->SetBinContent(36,23);
   h4->SetBinContent(37,33);
   h4->SetBinContent(38,19);
   h4->SetBinContent(39,19);
   h4->SetBinContent(40,16);
   h4->SetBinContent(41,21);
   h4->SetBinContent(42,12);
   h4->SetBinContent(43,17);
   h4->SetBinContent(44,16);
   h4->SetBinContent(45,16);
   h4->SetBinContent(46,19);
   h4->SetBinContent(47,11);
   h4->SetBinContent(48,16);
   h4->SetBinContent(49,12);
   h4->SetBinContent(50,18);
   h4->SetMinimum(0);
   h4->SetMaximum(147.6);
   h4->SetEntries(1959);
   h4->SetStats(0);
   
   TF1 *f4 = new TF1("f4","[0]*([7]*Gaus(x,[1],[2])/(sqrt(2*3.14159)*[2])+(1-[7])*Gaus(x,[1],[8])/(sqrt(2*3.14159)*[8]))+[3]+[4]*x+[6]*(38.42*Gaus(x,5.25,0.03473)+15.04*Gaus(x,5.25,0.1121)+104.3*Gaus(x,5.026,0.0935))",5,6);
   f4->SetFillColor(19);
   f4->SetFillStyle(0);
   f4->SetMarkerStyle(20);
   f4->SetLineColor(2);
   f4->SetLineWidth(1);
   f4->SetChisquare(40.98578);
   f4->SetNDF(44);
   f4->GetXaxis()->SetLabelFont(42);
   f4->GetXaxis()->SetLabelOffset(0.007);
   f4->GetXaxis()->SetLabelSize(0.05);
   f4->GetXaxis()->SetTitleSize(0.06);
   f4->GetXaxis()->SetTitleOffset(0.9);
   f4->GetXaxis()->SetTitleFont(42);
   f4->GetYaxis()->SetLabelFont(42);
   f4->GetYaxis()->SetLabelOffset(0.007);
   f4->GetYaxis()->SetLabelSize(0.05);
   f4->GetYaxis()->SetTitleSize(0.06);
   f4->GetYaxis()->SetTitleOffset(1.05);
   f4->GetYaxis()->SetTitleFont(42);
   f4->SetParameter(0,5.159699);
   f4->SetParError(0,0.4573041);
   f4->SetParLimits(0,0,0);
   f4->SetParameter(1,5.280027);
   f4->SetParError(1,0.002439908);
   f4->SetParLimits(1,0,0);
   f4->SetParameter(2,0.05);
   f4->SetParError(2,0);
   f4->SetParLimits(2,0.05,0.05);
   f4->SetParameter(3,188.1797);
   f4->SetParError(3,0.7307479);
   f4->SetParLimits(3,0,0);
   f4->SetParameter(4,-29.36589);
   f4->SetParError(4,0.1292885);
   f4->SetParLimits(4,-1000,0);
   f4->SetParameter(5,0.0684761);
   f4->SetParError(5,8);
   f4->SetParLimits(5,0,0);
   f4->SetParameter(6,0.2913956);
   f4->SetParError(6,0.02850488);
   f4->SetParLimits(6,0,0);
   f4->SetParameter(7,0.4721161);
   f4->SetParError(7,0);
   f4->SetParLimits(7,0.4721161,0.4721161);
   f4->SetParameter(8,0.01445891);
   f4->SetParError(8,0);
   f4->SetParLimits(8,0.01445891,0.01445891);
   h4->GetListOfFunctions()->Add(f4);
   h4->SetLineStyle(0);
   h4->SetMarkerStyle(24);
   h4->SetMarkerSize(0.8);
   h4->GetXaxis()->SetTitle("M_{B} (GeV/c^{2})");
   h4->GetXaxis()->CenterTitle(true);
   h4->GetXaxis()->SetLabelFont(42);
   h4->GetXaxis()->SetLabelOffset(0.007);
   h4->GetXaxis()->SetLabelSize(0.05);
   h4->GetXaxis()->SetTitleSize(0.06);
   h4->GetXaxis()->SetTitleOffset(0.9);
   h4->GetXaxis()->SetTitleFont(42);
   h4->GetYaxis()->SetTitle("Entries / (20 MeV/c^{2})");
   h4->GetYaxis()->CenterTitle(true);
   h4->GetYaxis()->SetLabelFont(42);
   h4->GetYaxis()->SetLabelOffset(0.007);
   h4->GetYaxis()->SetLabelSize(0.05);
   h4->GetYaxis()->SetTitleSize(0.06);
   h4->GetYaxis()->SetTitleOffset(1.4);
   h4->GetYaxis()->SetTitleFont(42);
   h4->GetZaxis()->SetLabelFont(42);
   h4->GetZaxis()->SetLabelOffset(0.007);
   h4->GetZaxis()->SetLabelSize(0.05);
   h4->GetZaxis()->SetTitleSize(0.06);
   h4->GetZaxis()->SetTitleFont(42);
   h4->Draw("e");
   
   TF1 *fBkpi = new TF1("fBkpi","[0]*(38.42*Gaus(x,5.25,0.03473)+15.04*Gaus(x,5.25,0.1121)+104.3*Gaus(x,5.026,0.0935))",5,5.45);

   Int_t ci;   // for color index setting
   ci = TColor::GetColor("#00cc00");
   fBkpi->SetFillColor(ci);
   fBkpi->SetFillStyle(3005);
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
   fBkpi->SetParameter(0,0.2913956);
   fBkpi->SetParError(0,0);
   fBkpi->SetParLimits(0,0,0);
   fBkpi->Draw("same");
   
   TF1 *background4 = new TF1("background4","[0]+[1]*x+[3]*(38.42*Gaus(x,5.25,0.03473)+15.04*Gaus(x,5.25,0.1121)+104.3*Gaus(x,5.026,0.0935))",5,6);
   background4->SetFillColor(19);
   background4->SetFillStyle(0);
   background4->SetMarkerStyle(20);
   background4->SetLineColor(4);
   background4->SetLineWidth(1);
   background4->SetLineStyle(2);
   background4->GetXaxis()->SetLabelFont(42);
   background4->GetXaxis()->SetLabelOffset(0.007);
   background4->GetXaxis()->SetLabelSize(0.05);
   background4->GetXaxis()->SetTitleSize(0.06);
   background4->GetXaxis()->SetTitleOffset(0.9);
   background4->GetXaxis()->SetTitleFont(42);
   background4->GetYaxis()->SetLabelFont(42);
   background4->GetYaxis()->SetLabelOffset(0.007);
   background4->GetYaxis()->SetLabelSize(0.05);
   background4->GetYaxis()->SetTitleSize(0.06);
   background4->GetYaxis()->SetTitleOffset(1.05);
   background4->GetYaxis()->SetTitleFont(42);
   background4->SetParameter(0,188.1797);
   background4->SetParError(0,0);
   background4->SetParLimits(0,0,0);
   background4->SetParameter(1,-29.36589);
   background4->SetParError(1,0);
   background4->SetParLimits(1,0,0);
   background4->SetParameter(2,0.0684761);
   background4->SetParError(2,0);
   background4->SetParLimits(2,0,0);
   background4->SetParameter(3,0.2913956);
   background4->SetParError(3,0);
   background4->SetParLimits(3,0,0);
   background4->Draw("same");
   
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
   fmass->SetParameter(0,5.159699);
   fmass->SetParError(0,0.4573041);
   fmass->SetParLimits(0,0,0);
   fmass->SetParameter(1,5.280027);
   fmass->SetParError(1,0.002439908);
   fmass->SetParLimits(1,0,0);
   fmass->SetParameter(2,0.05);
   fmass->SetParError(2,0);
   fmass->SetParLimits(2,0,0);
   fmass->SetParameter(3,0.4721161);
   fmass->SetParError(3,0);
   fmass->SetParLimits(3,0,0);
   fmass->SetParameter(4,0.01445891);
   fmass->SetParError(4,0);
   fmass->SetParLimits(4,0,0);
   fmass->Draw("same");
   
   TF1 *f4 = new TF1("f4","[0]*([7]*Gaus(x,[1],[2])/(sqrt(2*3.14159)*[2])+(1-[7])*Gaus(x,[1],[8])/(sqrt(2*3.14159)*[8]))+[3]+[4]*x+[6]*(38.42*Gaus(x,5.25,0.03473)+15.04*Gaus(x,5.25,0.1121)+104.3*Gaus(x,5.026,0.0935))",0,1);
   f4->SetFillColor(19);
   f4->SetFillStyle(0);
   f4->SetMarkerStyle(20);
   f4->SetLineColor(2);
   f4->SetLineWidth(1);
   f4->SetChisquare(40.98578);
   f4->SetNDF(44);
   f4->GetXaxis()->SetLabelFont(42);
   f4->GetXaxis()->SetLabelOffset(0.007);
   f4->GetXaxis()->SetLabelSize(0.05);
   f4->GetXaxis()->SetTitleSize(0.06);
   f4->GetXaxis()->SetTitleOffset(0.9);
   f4->GetXaxis()->SetTitleFont(42);
   f4->GetYaxis()->SetLabelFont(42);
   f4->GetYaxis()->SetLabelOffset(0.007);
   f4->GetYaxis()->SetLabelSize(0.05);
   f4->GetYaxis()->SetTitleSize(0.06);
   f4->GetYaxis()->SetTitleOffset(1.05);
   f4->GetYaxis()->SetTitleFont(42);
   f4->SetParameter(0,5.159699);
   f4->SetParError(0,0.4573041);
   f4->SetParLimits(0,0,0);
   f4->SetParameter(1,5.280027);
   f4->SetParError(1,0.002439908);
   f4->SetParLimits(1,0,0);
   f4->SetParameter(2,0.05);
   f4->SetParError(2,0);
   f4->SetParLimits(2,0.05,0.05);
   f4->SetParameter(3,188.1797);
   f4->SetParError(3,0.7307479);
   f4->SetParLimits(3,0,0);
   f4->SetParameter(4,-29.36589);
   f4->SetParError(4,0.1292885);
   f4->SetParLimits(4,-1000,0);
   f4->SetParameter(5,0.0684761);
   f4->SetParError(5,8);
   f4->SetParLimits(5,0,0);
   f4->SetParameter(6,0.2913956);
   f4->SetParError(6,0.02850488);
   f4->SetParLimits(6,0,0);
   f4->SetParameter(7,0.4721161);
   f4->SetParError(7,0);
   f4->SetParLimits(7,0.4721161,0.4721161);
   f4->SetParameter(8,0.01445891);
   f4->SetParError(8,0);
   f4->SetParLimits(8,0.01445891,0.01445891);
   f4->Draw("same");
   
   TLegend *leg = new TLegend(nan,nan,nan,nan,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextFont(62);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(19);
   leg->SetFillStyle(0);
   TLegendEntry *entry=leg->AddEntry("h4","CMS Preliminary","");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("h4","p+Pb #sqrt{s_{NN}}= 5.02 TeV","");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("h4","10<p_{T}^{B}<60 GeV/c","");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("h4","Data","pl");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("f4","Fit","l");
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
   entry=leg->AddEntry("background4","Combinatorial Background","l");
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
   
   leg = new TLegend(4.243992e-313,5.09279e-313,5.729389e-313,6.153788e-313,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextFont(62);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(19);
   leg->SetFillStyle(0);
   entry=leg->AddEntry("h4","B meson","");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("h4","M_{B}=5280.03 #pm 2.44 MeV/c^{2}","");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("h4","N_{B}=258 #pm 23","");
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
   TText *text = pt->AddText("[0]*([7]*Gaus(x,[1],[2])/(sqrt(2*3.14159)*[2])+(1-[7])*Gaus(x,[1],[8])/(sqrt(2*3.14159)*[8]))+[3]+[4]*x+[6]*(38.42*Gaus(x,5.25,0.03473)+15.04*Gaus(x,5.25,0.1121)+104.3*Gaus(x,5.026,0.0935))");
   pt->Draw();
   c4->Modified();
   c4->cd();
   c4->SetSelected(c4);
}
