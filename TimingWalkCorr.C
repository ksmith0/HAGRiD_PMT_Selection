#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TSpectrum.h"

//#define WALKCORR
#define PEAKGATE

TProfile* SetCorrAlias(TTree *t, const char* paramToCorr, const char* det, const char* cut);
void FindPeaks(TTree* t, const char* det, const char* cut, double *&peak, double *&stddev);

const int NumStdDevScans = 2;

void TimingWalkCorr(const char* filename, int chStart, int chLabr) {
	gStyle->SetOptFit(1111);


	TFile *file = new TFile(filename);
	TTree *t = (TTree*)file->Get("labr3");
	TCanvas *c = new TCanvas("c","");
	int numRows = 0;
#ifdef WALKCORR
	numRows += 2;
#endif
#ifdef PEAKGATE
	numRows += 1;
#endif

	if (numRows == 1) c->Divide(numRows,3);
	else c->Divide(3,numRows);


	std::string locCut = Form("loc_start==%d && loc_labr==%d",chStart,chLabr);

	c->cd(1)->SetGridx();
	t->Draw("tdiff>>hRaw(1000,0,100)",Form("en_start>5000 && %s",locCut.c_str()));
	TH1F* hRaw;
	for (int i=0;i<NumStdDevScans;i++) {
		hRaw = (TH1F*)gROOT->FindObject("hRaw");
		t->Draw(Form("tdiff>>hRaw(700,%f,%f)",hRaw->GetMean() - 3 * hRaw->GetStdDev(), hRaw->GetMean() + 3 * hRaw->GetStdDev()),Form("en_start>5000 && %s",locCut.c_str()));
	}
	hRaw = (TH1F*)gROOT->FindObject("hRaw");
	hRaw->Fit("gaus");
	c->Update();

#ifdef WALKCORR
	
	c->cd(2)->SetGridy();
	SetCorrAlias(t,"tdiff", "start", locCut.c_str());
	c->Update();
	c->cd(3)->SetGridy();
	t->Draw("tdiff_startCorr:en_start>>hStart(200,0,30000,200,0,100)",Form("en_start>1000 && %s",locCut.c_str()),"GOFF");
	TH2F* h;
	for (int i=0;i<NumStdDevScans;i++) {
		h = (TH2F*) gROOT->FindObject("hStart");
		t->Draw(Form("tdiff_startCorr:en_start>>hStart(200,0,30000,200,%f,%f)",h->GetMean(2) - 3 * h->GetStdDev(2),h->GetMean(2) + 3 * h->GetStdDev(2)),Form("en_start>1000 && %s",locCut.c_str()),"GOFF");
	}
	h = (TH2F*) gROOT->FindObject("hStart");
	h->SetMinimum(1);
	h->Draw("COLZ");
	c->Update();


	c->cd(5)->SetGridy();
	SetCorrAlias(t,"tdiff_startCorr", "labr",locCut.c_str());
	c->Update();
	c->cd(6)->SetGridy();
	t->Draw("tdiff_startCorr_labrCorr:en_labr>>hLabr(200,0,30000,200,0,100)",Form("en_start>1000 && %s",locCut.c_str()),"GOFF");
	for (int i=0;i<NumStdDevScans;i++) {
		h = (TH2F*) gROOT->FindObject("hLabr");
		t->Draw(Form("tdiff_startCorr_labrCorr:en_labr>>hLabr(200,0,30000,200,%f,%f)",h->GetMean(2) - 3 * h->GetStdDev(2),h->GetMean(2) + 3 * h->GetStdDev(2)),Form("en_start>1000 && %s",locCut.c_str()),"GOFF");
	}
	h = (TH2F*) gROOT->FindObject("hLabr");
	h->SetMinimum(0);
	h->Draw("COLZ");
	c->Update();

	c->cd(4)->SetGridx();
	TH1F* hCorr;
	/*
	t->Draw("tdiff_startCorr_labrCorr>>hCorr(500,0,100)",Form("en_start>5000 && en_labr>5000 && %s",locCut.c_str()));
   for (int i=0;i<NumStdDevScans;i++) {
		hCorr = (TH1F*)gROOT->FindObject("hCorr");
		t->Draw(Form("tdiff_startCorr_labrCorr>>hCorr(500,%f,%f)",hCorr->GetMean() - 3 * hCorr->GetStdDev(), hCorr->GetMean() + 3 * hCorr->GetStdDev()),Form("en_start>5000 && en_labr>5000 && %s",locCut.c_str()));
	}
	*/
	t->Draw(Form("tdiff_startCorr_labrCorr>>hCorr(700,%f,%f)",hRaw->GetXaxis()->GetXmin(), hRaw->GetXaxis()->GetXmax()),Form("en_start>5000 && en_labr>5000 && %s",locCut.c_str()));
	hCorr = (TH1F*)gROOT->FindObject("hCorr");
	TF1 *f = new TF1("gausFit","gaus",hCorr->GetMean() - hCorr->GetStdDev(), hCorr->GetMean() + hCorr->GetStdDev());
	hCorr->Fit(f,"ME");
	c->Update();
#endif

#ifdef PEAKGATE
	c->cd(3*numRows - 1);
	double *peaks, *stddev;
	FindPeaks(t,"start",Form("loc_start==%d",chStart),peaks,stddev);
	c->Update();
	

	c->cd(3*numRows);
	double *peaks_labr, *stddev_labr;
	FindPeaks(t,"labr",Form("loc_labr==%d",chLabr),peaks_labr,stddev_labr);
	c->Update();

	c->cd(3*numRows - 2);
	const char *cut1 = Form("en_start > %f && en_start < %f && en_labr > %f && en_labr < %f",
		peaks[0] - 2 * stddev[0], peaks[0] + 2 * stddev[0],
		peaks_labr[1] - 2 * stddev_labr[1], peaks_labr[1] + 2 * stddev_labr[1]);
	const char *cut2 = Form("en_start > %f && en_start < %f && en_labr > %f && en_labr < %f",
		peaks[1] - 2 * stddev[1], peaks[1] + 2 * stddev[1],
		peaks_labr[0] - 2 * stddev_labr[0], peaks_labr[0] + 2 * stddev_labr[0]);
	TH1F *hCorr2_1, *hCorr2_2;
	/*
	t->Draw("tdiff>>hCorr2_1(50,-100,100)",Form("%s && %s",locCut.c_str(),cut1));
	t->Draw("tdiff>>hCorr2_2(50,-100,100)",Form("%s && %s",locCut.c_str(),cut2));
	for (int i=0;i<NumStdDevScans;i++) {
		hCorr2_1 = (TH1F*)gROOT->FindObject("hCorr2_1");
		t->Draw(Form("tdiff>>hCorr2_1(200,%f,%f)",hRaw->GetMean() - 5 * hRaw->GetStdDev(), hRaw->GetMean() + 5 * hRaw->GetStdDev()),Form("%s && %s",locCut.c_str(),cut1),"GOFF");

		hCorr2_2 = (TH1F*)gROOT->FindObject("hCorr2_2");
		t->Draw(Form("tdiff>>hCorr2_2(200,%f,%f)",hRaw->GetMean() - 5 * hRaw->GetStdDev(), hRaw->GetMean() + 5 * hRaw->GetStdDev()),Form("%s && %s",locCut.c_str(),cut2),"GOFF");
	}
	*/
#ifdef WALKCORR
	const char *tdiffString = "tdiff_startCorr_labrCorr";
#else
	const char *tdiffString = "tdiff";
#endif
//	t->Draw(Form("%s>>hCorr2_1(1000,-20,20)",tdiffString),Form("%s && %s",locCut.c_str(),cut1),"GOFF");
//	t->Draw(Form("%s>>hCorr2_2(1000,-20,20)",tdiffString),Form("%s && %s",locCut.c_str(),cut2),"GOFF");
	t->Draw(Form("%s>>hCorr2_1(200,%f,%f)",tdiffString,hRaw->GetMean()-2,hRaw->GetMean()+2),Form("%s && %s",locCut.c_str(),cut1),"GOFF");
	t->Draw(Form("%s>>hCorr2_2(200,%f,%f)",tdiffString,hRaw->GetMean()-2,hRaw->GetMean()+2),Form("%s && %s",locCut.c_str(),cut2),"GOFF");
	hCorr2_1 = (TH1F*)gROOT->FindObject("hCorr2_1");
	hCorr2_2 = (TH1F*)gROOT->FindObject("hCorr2_2");
	TF1 *f2_1 = new TF1("gausFit","gaus",hCorr2_1->GetMean() - hCorr2_1->GetStdDev(), hCorr2_1->GetMean() + hCorr2_1->GetStdDev());
	f2_1->SetLineColor(kBlue);
	hCorr2_1->Fit(f2_1,"ME");
	TF1 *f2_2 = new TF1("gausFit","gaus",hCorr2_2->GetMean() - hCorr2_2->GetStdDev(), hCorr2_2->GetMean() + hCorr2_2->GetStdDev());
	hCorr2_2->Fit(f2_2,"ME");
	c->Update();
	TPaveStats* stats = (TPaveStats*) hCorr2_2->GetListOfFunctions()->FindObject("stats");
	if (stats) {
		float boxHeight = stats->GetY2NDC() - stats->GetY1NDC();
		stats->SetY2NDC(stats->GetY1NDC() - 0.01);
		stats->SetY1NDC(stats->GetY2NDC() - boxHeight);
	}

	hCorr2_1->Draw("SAMES");
	hCorr2_2->SetLineColor(kRed);

	
	c->Update();
#endif
}

void FindPeaks(TTree *t, const char* det,const char* cut, double *&peaks, double *&stddev) {
	t->Draw(Form("en_%s>>hPeakFind_%s(500)",det,det),cut,"");
	TH1F* hist = (TH1F*) gROOT->FindObject(Form("hPeakFind_%s",det));
	
	TSpectrum *s = new TSpectrum(100);

	s->Search(hist,2);

	double highEn, lowEn;
	double highEnHeight, lowEnHeight;
	bool found = false;
	for (int i=0;i<s->GetNPeaks();i++) {
		for (int j=0;j<s->GetNPeaks();j++) {
			printf("%d %d %f %f %f %f %f\n",i,j,s->GetPositionX()[i],s->GetPositionY()[i],s->GetPositionX()[j], s->GetPositionY()[j], (s->GetPositionX()[i]/s->GetPositionX()[j]/1.136) -1);
			if (i==j) continue;
			if (fabs((s->GetPositionX()[i] / s->GetPositionX()[j] / 1.136 ) - 1) < 0.03 && s->GetPositionY()[i] < s->GetPositionY()[j]) {
				printf("Found Correct Ratio\n");
				if (found) {
					if (s->GetPositionY()[i] >= highEnHeight || s->GetPositionY()[j] >= lowEnHeight) {
						printf("Resetting Values!\n");
						highEn = s->GetPositionX()[i];
						highEnHeight = s->GetPositionY()[i];
						lowEn = s->GetPositionX()[j];
						lowEnHeight = s->GetPositionY()[j];
					}
				}
				else {
					printf("Setting Values!\n");
					found = true;
					highEn = s->GetPositionX()[i];
					highEnHeight = s->GetPositionY()[i];
					lowEn = s->GetPositionX()[j];
					lowEnHeight = s->GetPositionY()[j];
				}
			}
		}
	}

	TF1 *fHigh = new TF1("PeakFitHigh","gaus");
	fHigh->SetParameter(1,highEn);
	fHigh->SetRange(highEn - 0.026 * highEn, highEn + 0.026 * highEn);
	hist->Fit(fHigh,"QMER+");
	TF1 *fLow = new TF1("PeakFitLow","gaus");
	fLow->SetParameter(1,lowEn);
	fLow->SetRange(lowEn - 0.026 * lowEn, lowEn + 0.026 * lowEn);
	hist->Fit(fLow,"QMER+");

	peaks = new Double_t[2];
	stddev = new Double_t[2];

	peaks[0] = fLow->GetParameter(1);
	stddev[0] = fLow->GetParameter(2);
	peaks[1] = fHigh->GetParameter(1);
	stddev[1] = fHigh->GetParameter(2);
}

TProfile* SetCorrAlias(TTree *t, const char* paramToCorr, const char* det, const char* cut) {
	//t->Draw(Form("%s:en_%s>>(200,0,1,200,0,1)",paramToCorr,det),"en_start>1000","GOFF");
	t->Draw(Form("%s:en_%s>>hWalkCorr_%s(100,0,30000,400,-100,100)",paramToCorr,det,det),Form("en_start>1000 && %s",cut),"COLZ");
	TH2F *h;
	for (int i=0;i<NumStdDevScans;i++) {
		h = (TH2F*) gROOT->FindObject(Form("hWalkCorr_%s",det));
		t->Draw(Form("%s:en_%s>>hWalkCorr_%s(100,0,30000,200,%f,%f)",paramToCorr,det,det,h->GetMean(2) - 3 * h->GetStdDev(2),h->GetMean(2) + 3 * h->GetStdDev(2)),Form("en_start>1000 && %s",cut),"COLZ");
	}

	h = (TH2F*) gROOT->FindObject(Form("hWalkCorr_%s",det));
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
