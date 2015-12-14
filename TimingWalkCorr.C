#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"

void SetCorrAlias(TTree *t, const char* paramToCorr, const char* det);

void TimingWalkCorr(const char* filename) {
	gStyle->SetOptFit(1111);


	TFile *file = new TFile(filename);
	TTree *t = (TTree*)file->Get("labr3");
	TCanvas *c = new TCanvas("c","");
	c->Divide(3,2);

	c->cd(1)->SetGridx();
	t->Draw("tdiff>>hRaw(500,0,100)","en_start>5000");
	TH1F* hRaw = (TH1F*)gROOT->FindObject("hRaw");
	t->Draw(Form("tdiff>>hRaw(500,%f,%f)",hRaw->GetMean() - 2 * hRaw->GetStdDev(), hRaw->GetMean() + 2 * hRaw->GetStdDev()),"en_start>5000");
	TH1F* hRaw = (TH1F*)gROOT->FindObject("hRaw");
	hRaw->Fit("gaus");
	c->Update();

	c->cd(2)->SetGridy();
	SetCorrAlias(t,"tdiff", "start");
	c->Update();
	c->cd(3)->SetGridy();
	t->Draw("tdiff_startCorr:en_start>>hStart(200,0,30000,200,0,100)","en_start>1000","COLZ");
	TH2F *h = (TH2F*) gROOT->FindObject("hStart");
	t->Draw(Form("tdiff_startCorr:en_start>>hStart(200,0,30000,200,%f,%f)",h->GetMean(2) - 3 * h->GetStdDev(2),h->GetMean(2) + 3 * h->GetStdDev(2)),"en_start>1000","COLZ");
	c->Update();

	c->cd(5)->SetGridy();
	SetCorrAlias(t,"tdiff_startCorr", "labr");
	c->Update();
	c->cd(6)->SetGridy();
	t->Draw("tdiff_startCorr_labrCorr:en_labr>>hLabr(200,0,30000,200,0,100)","en_start>1000","COLZ");
	h = (TH2F*) gROOT->FindObject("hLabr");
	t->Draw(Form("tdiff_startCorr_labrCorr:en_labr>>hLabr(200,0,30000,200,%f,%f)",h->GetMean(2) - 3 * h->GetStdDev(2),h->GetMean(2) + 3 * h->GetStdDev(2)),"en_start>1000","COLZ");
	c->Update();

	c->cd(4)->SetGridx();
	t->Draw("tdiff_startCorr_labrCorr>>hCorr(500,0,100)","en_start>5000");
	TH1F* hCorr = (TH1F*)gROOT->FindObject("hCorr");
	t->Draw(Form("tdiff_startCorr_labrCorr>>hCorr(500,%f,%f)",hCorr->GetMean() - 3 * hCorr->GetStdDev(), hCorr->GetMean() + 3 * hCorr->GetStdDev()),"en_start>5000 && en_labr>5000");
	TH1F* hCorr = (TH1F*)gROOT->FindObject("hCorr");
	TF1 *f = new TF1("gausFit","gaus",hCorr->GetMean() - hCorr->GetStdDev(), hCorr->GetMean() + hCorr->GetStdDev());
	hCorr->Fit(f,"RME");
	c->Update();
	

}

TProfile* SetCorrAlias(TTree *t, const char* paramToCorr, const char* det) {
	//t->Draw(Form("%s:en_%s>>(200,0,1,200,0,1)",paramToCorr,det),"en_start>1000","GOFF");
	t->Draw(Form("%s:en_%s>>hWalkCorr_%s(100,0,30000,200,0,100)",paramToCorr,det,det),"en_start>1000","COLZ");
	TH2F *h = (TH2F*) gROOT->FindObject(Form("hWalkCorr_%s",det));
	t->Draw(Form("%s:en_%s>>hWalkCorr_%s(100,0,30000,200,%f,%f)",paramToCorr,det,det,h->GetMean(2) - 3 * h->GetStdDev(2),h->GetMean(2) + 3 * h->GetStdDev(2)),"en_start>1000","COLZ");

	TH2F *h = (TH2F*) gROOT->FindObject(Form("hWalkCorr_%s",det));
	TProfile *prof = h->ProfileX(Form("p_%s",det));
	float min = h->GetXaxis()->GetXmin();
	float max = h->GetXaxis()->GetXmax();
	float range = max - min;

	TF1 *fit = new TF1("fit","pol3",min+0.1*range, max - 0.3 * range);
	prof->Fit(fit,"RME");

	std::string corr = paramToCorr;
	corr.append("-(");
	for (int i=1;i<4;i++) {
		corr.append(Form("+%e * pow(en_%s,%d)",fit->GetParameter(i),det,i));
	}
	corr.append(")");
	t->SetAlias(Form("%s_%sCorr",paramToCorr,det),corr.c_str());

	h->SetMarkerColor(kGreen+1);

	h->Draw();
	h->GetYaxis()->SetRangeUser(h->GetMean(2) - 2 * h->GetStdDev(2),h->GetMean(2) + 2 * h->GetStdDev(2));
	prof->Draw("SAME");

	delete fit;

	return prof;
}
