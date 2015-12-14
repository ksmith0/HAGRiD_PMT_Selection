#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#ifndef STDLIB
	#include <vector>
	#include <utility>
	#define STDLIB
#endif
#include "TSpectrum.h"

void SetCorrAlias(TTree *t, const char* paramToCorr, const char* det, const char* cut);
std::vector<std::pair<float,float>> FindPeaks(TTree* t, const char* var);

void TimingWalkCorr(const char* filename, int chStart, int chLabr) {
	gStyle->SetOptFit(1111);


	TFile *file = new TFile(filename);
	TTree *t = (TTree*)file->Get("labr3");
	TCanvas *c = new TCanvas("c","");
	c->Divide(3,3);

	std::string locCut = Form("loc_start==%d && loc_labr==%d",chStart,chLabr);

	c->cd(1)->SetGridx();
	t->Draw("tdiff>>hRaw(500,0,100)",Form("en_start>5000 && %s",locCut.c_str()));
	TH1F* hRaw;
	for (int i=0;i<3;i++) {
		hRaw = (TH1F*)gROOT->FindObject("hRaw");
		t->Draw(Form("tdiff>>hRaw(500,%f,%f)",hRaw->GetMean() - 2 * hRaw->GetStdDev(), hRaw->GetMean() + 2 * hRaw->GetStdDev()),Form("en_start>5000 && %s",locCut.c_str()));
	}
	hRaw = (TH1F*)gROOT->FindObject("hRaw");
	hRaw->Fit("gaus");
	c->Update();

	c->cd(2)->SetGridy();
	SetCorrAlias(t,"tdiff", "start", locCut.c_str());
	c->Update();
	c->cd(3)->SetGridy();
	t->Draw("tdiff_startCorr:en_start>>hStart(200,0,30000,200,0,100)",Form("en_start>1000 && %s",locCut.c_str()),"GOFF");
	TH2F* h;
	for (int i=0;i<4;i++) {
		h = (TH2F*) gROOT->FindObject("hStart");
		t->Draw(Form("tdiff_startCorr:en_start>>hStart(200,0,30000,200,%f,%f)",h->GetMean(2) - 3 * h->GetStdDev(2),h->GetMean(2) + 3 * h->GetStdDev(2)),Form("en_start>1000 && %s",locCut.c_str()),"GOFF");
	}
	h->SetMinimum(0);
	h->Draw("COLZ");
	c->Update();

	c->cd(5)->SetGridy();
	SetCorrAlias(t,"tdiff_startCorr", "labr",locCut.c_str());
	c->Update();
	c->cd(6)->SetGridy();
	t->Draw("tdiff_startCorr_labrCorr:en_labr>>hLabr(200,0,30000,200,0,100)","en_start>1000","COLZ");
	h = (TH2F*) gROOT->FindObject("hLabr");
	t->Draw(Form("tdiff_startCorr_labrCorr:en_labr>>hLabr(200,0,30000,200,%f,%f)",h->GetMean(2) - 3 * h->GetStdDev(2),h->GetMean(2) + 3 * h->GetStdDev(2)),"en_start>1000","COLZ");
	c->Update();

	c->cd(4)->SetGridx();
	t->Draw("tdiff_startCorr_labrCorr>>hCorr(500,0,100)","en_start>5000");
	TH1F* hCorr;
	for (int i=0;i<3;i++) {
		hCorr = (TH1F*)gROOT->FindObject("hCorr");
		t->Draw(Form("tdiff_startCorr_labrCorr>>hCorr(500,%f,%f)",hCorr->GetMean() - 3 * hCorr->GetStdDev(), hCorr->GetMean() + 3 * hCorr->GetStdDev()),"en_start>5000 && en_labr>5000");
	}
	hCorr = (TH1F*)gROOT->FindObject("hCorr");
	TF1 *f = new TF1("gausFit","gaus",hCorr->GetMean() - hCorr->GetStdDev(), hCorr->GetMean() + hCorr->GetStdDev());
	hCorr->Fit(f,"RME");


	c->cd(8);
	std::vector<std::pair<float,float>> peaks = FindPeaks(t,"en_start");
	

	c->cd(9);

	c->Update();
	

}

std::vector<std::pair<float,float>> FindPeaks(TTree *t, const char* var) {
	std::vector<std::pair<float,float>> peaks;
	t->Draw(Form("%s>>hPeakFind(500)",var),"","");
	TH1F* hist = (TH1F*) gROOT->FindObject("hPeakFind");
	
	TSpectrum *s = new TSpectrum(100);

	s->Search(hist,10);


	return peaks; 
}

TProfile* SetCorrAlias(TTree *t, const char* paramToCorr, const char* det, const char* cut) {
	//t->Draw(Form("%s:en_%s>>(200,0,1,200,0,1)",paramToCorr,det),"en_start>1000","GOFF");
	t->Draw(Form("%s:en_%s>>hWalkCorr_%s(100,0,30000,200,0,100)",paramToCorr,det,det),Form("en_start>1000 && %s",cut),"COLZ");
	TH2F *h;
	for (int i=0;i<3;i++) {
		h = (TH2F*) gROOT->FindObject(Form("hWalkCorr_%s",det));
		t->Draw(Form("%s:en_%s>>hWalkCorr_%s(100,0,30000,200,%f,%f)",paramToCorr,det,det,h->GetMean(2) - 3 * h->GetStdDev(2),h->GetMean(2) + 3 * h->GetStdDev(2)),Form("en_start>1000 && %s",cut),"COLZ");
	}

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
