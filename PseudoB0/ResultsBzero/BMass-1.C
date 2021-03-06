{
//=========Macro generated from canvas: c1/
//=========  (Sat Apr 26 01:55:31 2014) by ROOT version5.32/00
   TCanvas *c1 = new TCanvas("c1", "",100,122,600,600);
   gStyle->SetOptFit(1);
   gStyle->SetOptStat(0);
   c1->Range(0,0,1,1);
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetTickx(1);
   c1->SetTicky(1);
   c1->SetLeftMargin(0.13);
   c1->SetRightMargin(0.05);
   c1->SetTopMargin(0.05);
   c1->SetBottomMargin(0.13);
   c1->SetFrameFillStyle(0);
   c1->SetFrameBorderMode(0);
   
   TH1D *h1 = new TH1D("h1","",50,5,6);
   h1->SetBinContent(1,34);
   h1->SetBinContent(2,27);
   h1->SetBinContent(3,27);
   h1->SetBinContent(4,31);
   h1->SetBinContent(5,29);
   h1->SetBinContent(6,26);
   h1->SetBinContent(7,27);
   h1->SetBinContent(8,18);
   h1->SetBinContent(9,16);
   h1->SetBinContent(10,25);
   h1->SetBinContent(11,22);
   h1->SetBinContent(12,27);
   h1->SetBinContent(13,37);
   h1->SetBinContent(14,47);
   h1->SetBinContent(15,62);
   h1->SetBinContent(16,30);
   h1->SetBinContent(17,18);
   h1->SetBinContent(18,15);
   h1->SetBinContent(19,15);
   h1->SetBinContent(20,13);
   h1->SetBinContent(21,9);
   h1->SetBinContent(22,14);
   h1->SetBinContent(23,3);
   h1->SetBinContent(24,8);
   h1->SetBinContent(25,4);
   h1->SetBinContent(26,8);
   h1->SetBinContent(27,8);
   h1->SetBinContent(28,4);
   h1->SetBinContent(29,4);
   h1->SetBinContent(30,6);
   h1->SetBinContent(31,5);
   h1->SetBinContent(32,6);
   h1->SetBinContent(33,6);
   h1->SetBinContent(34,5);
   h1->SetBinContent(35,6);
   h1->SetBinContent(36,7);
   h1->SetBinContent(37,4);
   h1->SetBinContent(38,2);
   h1->SetBinContent(40,2);
   h1->SetBinContent(41,3);
   h1->SetBinContent(42,3);
   h1->SetBinContent(43,3);
   h1->SetBinContent(44,4);
   h1->SetBinContent(45,7);
   h1->SetBinContent(46,4);
   h1->SetBinContent(47,4);
   h1->SetBinContent(48,2);
   h1->SetBinContent(49,4);
   h1->SetBinContent(50,2);
   h1->SetBinError(1,5.830952);
   h1->SetBinError(2,5.196152);
   h1->SetBinError(3,5.196152);
   h1->SetBinError(4,5.567764);
   h1->SetBinError(5,5.385165);
   h1->SetBinError(6,5.09902);
   h1->SetBinError(7,5.196152);
   h1->SetBinError(8,4.242641);
   h1->SetBinError(9,4);
   h1->SetBinError(10,5);
   h1->SetBinError(11,4.690416);
   h1->SetBinError(12,5.196152);
   h1->SetBinError(13,6.082763);
   h1->SetBinError(14,6.855655);
   h1->SetBinError(15,7.874008);
   h1->SetBinError(16,5.477226);
   h1->SetBinError(17,4.242641);
   h1->SetBinError(18,3.872983);
   h1->SetBinError(19,3.872983);
   h1->SetBinError(20,3.605551);
   h1->SetBinError(21,3);
   h1->SetBinError(22,3.741657);
   h1->SetBinError(23,1.732051);
   h1->SetBinError(24,2.828427);
   h1->SetBinError(25,2);
   h1->SetBinError(26,2.828427);
   h1->SetBinError(27,2.828427);
   h1->SetBinError(28,2);
   h1->SetBinError(29,2);
   h1->SetBinError(30,2.44949);
   h1->SetBinError(31,2.236068);
   h1->SetBinError(32,2.44949);
   h1->SetBinError(33,2.44949);
   h1->SetBinError(34,2.236068);
   h1->SetBinError(35,2.44949);
   h1->SetBinError(36,2.645751);
   h1->SetBinError(37,2);
   h1->SetBinError(38,1.414214);
   h1->SetBinError(39,1);
   h1->SetBinError(40,1.414214);
   h1->SetBinError(41,1.732051);
   h1->SetBinError(42,1.732051);
   h1->SetBinError(43,1.732051);
   h1->SetBinError(44,2);
   h1->SetBinError(45,2.645751);
   h1->SetBinError(46,2);
   h1->SetBinError(47,2);
   h1->SetBinError(48,1.414214);
   h1->SetBinError(49,2);
   h1->SetBinError(50,1.414214);
   h1->SetMinimum(0);
   h1->SetMaximum(74.4);
   h1->SetEntries(693);
   h1->SetStats(0);
   
   TF1 *f1 = new TF1("f1","[0]*([7]*Gaus(x,[1],[2])/(sqrt(2*3.14159)*[2])+(1-[7])*Gaus(x,[1],[8])/(sqrt(2*3.14159)*[8]))+[3]+[4]*x+[6]*(6.71675e+00*Gaus(x,5.30142e+00,8.42680e-02)/(sqrt(2*3.14159)*8.42680e-02)+4.06744e+01*Gaus(x,5.00954e+00,8.11305e-02)/(sqrt(2*3.14159)*8.11305e-02)+5.99974e-01*(2.376716*Gaus(x,5.640619,0.095530)/(sqrt(2*3.14159)*0.095530)+3.702342*Gaus(x,5.501706,0.046222)/(sqrt(2*3.14159)*0.046222))+1.31767e-01*(61.195688*Gaus(x,5.127566,0.087439)/(sqrt(2*3.14159)*0.087439)+58.943919*Gaus(x,5.246471,0.041983)/(sqrt(2*3.14159)*0.041983)))",5,6);
   f1->SetFillColor(19);
   f1->SetFillStyle(0);
   f1->SetMarkerStyle(20);
   f1->SetLineColor(2);
   f1->SetLineWidth(1);
   f1->SetChisquare(53.31804);
   f1->SetNDF(44);
   f1->GetXaxis()->SetLabelFont(42);
   f1->GetXaxis()->SetLabelOffset(0.007);
   f1->GetXaxis()->SetLabelSize(0.05);
   f1->GetXaxis()->SetTitleSize(0.06);
   f1->GetXaxis()->SetTitleOffset(0.9);
   f1->GetXaxis()->SetTitleFont(42);
   f1->GetYaxis()->SetLabelFont(42);
   f1->GetYaxis()->SetLabelOffset(0.007);
   f1->GetYaxis()->SetLabelSize(0.05);
   f1->GetYaxis()->SetTitleSize(0.06);
   f1->GetYaxis()->SetTitleOffset(1.05);
   f1->GetYaxis()->SetTitleFont(42);
   f1->SetParameter(0,2.802961);
   f1->SetParError(0,0.3223323);
   f1->SetParLimits(0,0,0);
   f1->SetParameter(1,5.285265);
   f1->SetParError(1,0.003214787);
   f1->SetParLimits(1,0,0);
   f1->SetParameter(2,0.05);
   f1->SetParError(2,0);
   f1->SetParLimits(2,0.05,0.05);
   f1->SetParameter(3,34.49159);
   f1->SetParError(3,0.3617767);
   f1->SetParLimits(3,0,0);
   f1->SetParameter(4,-5.315307);
   f1->SetParError(4,0.06335287);
   f1->SetParLimits(4,-1000,0);
   f1->SetParameter(5,0);
   f1->SetParError(5,78.38379);
   f1->SetParLimits(5,0,0);
   f1->SetParameter(6,0.1211804);
   f1->SetParError(6,0.009513301);
   f1->SetParLimits(6,0,1000);
   f1->SetParameter(7,0.5827073);
   f1->SetParError(7,0);
   f1->SetParLimits(7,0.5827073,0.5827073);
   f1->SetParameter(8,0.01404381);
   f1->SetParError(8,0);
   f1->SetParLimits(8,0.01404381,0.01404381);
   h1->GetListOfFunctions()->Add(f1);
   h1->SetLineStyle(0);
   h1->SetMarkerStyle(24);
   h1->SetMarkerSize(0.8);
   h1->GetXaxis()->SetTitle("M_{B} (GeV/c^{2})");
   h1->GetXaxis()->CenterTitle(true);
   h1->GetXaxis()->SetLabelFont(42);
   h1->GetXaxis()->SetLabelOffset(0.007);
   h1->GetXaxis()->SetLabelSize(0.05);
   h1->GetXaxis()->SetTitleSize(0.06);
   h1->GetXaxis()->SetTitleOffset(0.9);
   h1->GetXaxis()->SetTitleFont(42);
   h1->GetYaxis()->SetTitle("Entries / (20 MeV/c^{2})");
   h1->GetYaxis()->CenterTitle(true);
   h1->GetYaxis()->SetLabelFont(42);
   h1->GetYaxis()->SetLabelOffset(0.007);
   h1->GetYaxis()->SetLabelSize(0.05);
   h1->GetYaxis()->SetTitleSize(0.06);
   h1->GetYaxis()->SetTitleOffset(1.4);
   h1->GetYaxis()->SetTitleFont(42);
   h1->GetZaxis()->SetLabelFont(42);
   h1->GetZaxis()->SetLabelOffset(0.007);
   h1->GetZaxis()->SetLabelSize(0.05);
   h1->GetZaxis()->SetTitleSize(0.06);
   h1->GetZaxis()->SetTitleFont(42);
   h1->Draw("e");
   
   TF1 *fBkpi = new TF1("fBkpi","[0]*(6.71675e+00*Gaus(x,5.30142e+00,8.42680e-02)/(sqrt(2*3.14159)*8.42680e-02)+4.06744e+01*Gaus(x,5.00954e+00,8.11305e-02)/(sqrt(2*3.14159)*8.11305e-02)+5.99974e-01*(2.376716*Gaus(x,5.640619,0.095530)/(sqrt(2*3.14159)*0.095530)+3.702342*Gaus(x,5.501706,0.046222)/(sqrt(2*3.14159)*0.046222))+1.31767e-01*(61.195688*Gaus(x,5.127566,0.087439)/(sqrt(2*3.14159)*0.087439)+58.943919*Gaus(x,5.246471,0.041983)/(sqrt(2*3.14159)*0.041983)))",5,5.45);

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
   fBkpi->SetParameter(0,0.1211804);
   fBkpi->SetParError(0,0);
   fBkpi->SetParLimits(0,0,0);
   fBkpi->Draw("same");
   
   TF1 *background1 = new TF1("background1","[0]+[1]*x",5,6);
   background1->SetFillColor(19);
   background1->SetFillStyle(0);
   background1->SetMarkerStyle(20);
   background1->SetLineColor(4);
   background1->SetLineWidth(1);
   background1->SetLineStyle(2);
   background1->GetXaxis()->SetLabelFont(42);
   background1->GetXaxis()->SetLabelOffset(0.007);
   background1->GetXaxis()->SetLabelSize(0.05);
   background1->GetXaxis()->SetTitleSize(0.06);
   background1->GetXaxis()->SetTitleOffset(0.9);
   background1->GetXaxis()->SetTitleFont(42);
   background1->GetYaxis()->SetLabelFont(42);
   background1->GetYaxis()->SetLabelOffset(0.007);
   background1->GetYaxis()->SetLabelSize(0.05);
   background1->GetYaxis()->SetTitleSize(0.06);
   background1->GetYaxis()->SetTitleOffset(1.05);
   background1->GetYaxis()->SetTitleFont(42);
   background1->SetParameter(0,34.49159);
   background1->SetParError(0,0);
   background1->SetParLimits(0,0,0);
   background1->SetParameter(1,-5.315307);
   background1->SetParError(1,0);
   background1->SetParLimits(1,0,0);
   background1->Draw("same");
   
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
   fmass->SetParameter(0,2.802961);
   fmass->SetParError(0,0.3223323);
   fmass->SetParLimits(0,0,0);
   fmass->SetParameter(1,5.285265);
   fmass->SetParError(1,0.003214787);
   fmass->SetParLimits(1,0,0);
   fmass->SetParameter(2,0.05);
   fmass->SetParError(2,0);
   fmass->SetParLimits(2,0,0);
   fmass->SetParameter(3,0.5827073);
   fmass->SetParError(3,0);
   fmass->SetParLimits(3,0,0);
   fmass->SetParameter(4,0.01404381);
   fmass->SetParError(4,0);
   fmass->SetParLimits(4,0,0);
   fmass->Draw("same");
   
   TF1 *f1 = new TF1("f1","[0]*([7]*Gaus(x,[1],[2])/(sqrt(2*3.14159)*[2])+(1-[7])*Gaus(x,[1],[8])/(sqrt(2*3.14159)*[8]))+[3]+[4]*x+[6]*(6.71675e+00*Gaus(x,5.30142e+00,8.42680e-02)/(sqrt(2*3.14159)*8.42680e-02)+4.06744e+01*Gaus(x,5.00954e+00,8.11305e-02)/(sqrt(2*3.14159)*8.11305e-02)+5.99974e-01*(2.376716*Gaus(x,5.640619,0.095530)/(sqrt(2*3.14159)*0.095530)+3.702342*Gaus(x,5.501706,0.046222)/(sqrt(2*3.14159)*0.046222))+1.31767e-01*(61.195688*Gaus(x,5.127566,0.087439)/(sqrt(2*3.14159)*0.087439)+58.943919*Gaus(x,5.246471,0.041983)/(sqrt(2*3.14159)*0.041983)))",0,1);
   f1->SetFillColor(19);
   f1->SetFillStyle(0);
   f1->SetMarkerStyle(20);
   f1->SetLineColor(2);
   f1->SetLineWidth(1);
   f1->SetChisquare(53.31804);
   f1->SetNDF(44);
   f1->GetXaxis()->SetLabelFont(42);
   f1->GetXaxis()->SetLabelOffset(0.007);
   f1->GetXaxis()->SetLabelSize(0.05);
   f1->GetXaxis()->SetTitleSize(0.06);
   f1->GetXaxis()->SetTitleOffset(0.9);
   f1->GetXaxis()->SetTitleFont(42);
   f1->GetYaxis()->SetLabelFont(42);
   f1->GetYaxis()->SetLabelOffset(0.007);
   f1->GetYaxis()->SetLabelSize(0.05);
   f1->GetYaxis()->SetTitleSize(0.06);
   f1->GetYaxis()->SetTitleOffset(1.05);
   f1->GetYaxis()->SetTitleFont(42);
   f1->SetParameter(0,2.802961);
   f1->SetParError(0,0.3223323);
   f1->SetParLimits(0,0,0);
   f1->SetParameter(1,5.285265);
   f1->SetParError(1,0.003214787);
   f1->SetParLimits(1,0,0);
   f1->SetParameter(2,0.05);
   f1->SetParError(2,0);
   f1->SetParLimits(2,0.05,0.05);
   f1->SetParameter(3,34.49159);
   f1->SetParError(3,0.3617767);
   f1->SetParLimits(3,0,0);
   f1->SetParameter(4,-5.315307);
   f1->SetParError(4,0.06335287);
   f1->SetParLimits(4,-1000,0);
   f1->SetParameter(5,0);
   f1->SetParError(5,78.38379);
   f1->SetParLimits(5,0,0);
   f1->SetParameter(6,0.1211804);
   f1->SetParError(6,0.009513301);
   f1->SetParLimits(6,0,1000);
   f1->SetParameter(7,0.5827073);
   f1->SetParError(7,0);
   f1->SetParLimits(7,0.5827073,0.5827073);
   f1->SetParameter(8,0.01404381);
   f1->SetParError(8,0);
   f1->SetParLimits(8,0.01404381,0.01404381);
   f1->Draw("same");
   
   TLegend *leg = new TLegend(0,0,0,0,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextFont(62);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(19);
   leg->SetFillStyle(0);
   TLegendEntry *entry=leg->AddEntry("h1","CMS Preliminary","");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("h1","p+Pb #sqrt{s_{NN}}= 5.02 TeV","");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("h1","10<p_{T}^{B}<60 GeV/c","");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("h1","Data","pl");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("f1","Fit","l");
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
   entry=leg->AddEntry("background1","Combinatorial Background","l");
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
   entry=leg->AddEntry("h1","B meson","");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("h1","M_{B}=5285.27 #pm 3.21 MeV/c^{2}","");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("h1","N_{B}=140 #pm 16","");
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
   TText *text = pt->AddText("[0]*([7]*Gaus(x,[1],[2])/(sqrt(2*3.14159)*[2])+(1-[7])*Gaus(x,[1],[8])/(sqrt(2*3.14159)*[8]))+[3]+[4]*x+[6]*(6.71675e+00*Gaus(x,5.30142e+00,8.42680e-02)/(sqrt(2*3.14159)*8.42680e-02)+4.06744e+01*Gaus(x,5.00954e+00,8.11305e-02)/(sqrt(2*3.14159)*8.11305e-02)+5.99974e-01*(2.376716*Gaus(x,5.640619,0.095530)/(sqrt(2*3.14159)*0.095530)+3.702342*Gaus(x,5.501706,0.046222)/(sqrt(2*3.14159)*0.046222))+1.31767e-01*(61.195688*Gaus(x,5.127566,0.087439)/(sqrt(2*3.14159)*0.087439)+58.943919*Gaus(x,5.246471,0.041983)/(sqrt(2*3.14159)*0.041983)))");
   pt->Draw();
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
}
