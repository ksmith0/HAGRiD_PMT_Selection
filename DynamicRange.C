void DynamicRange(const char *filename) {
	gStyle->SetOptFit(11111);
	TFile *f = new TFile(filename);
	TTree *t = (TTree*) f->Get("labr3");

	TCanvas *c = new TCanvas("c","Dynamic Range");
	c->DivideSquare(6);

	for (int det=0;det<7;det++) {
		int pad = det + 1;
		//Position 5 was a teeny vandle.
		if (det==5) continue;
		else if (det > 5) pad = det;
		c->cd(pad)->Divide(1,2);
		c->cd(pad)->cd(1);
		t->Draw(Form("peaken_labr+baseline_labr>>hDet%d(2048,0,4095)",det),Form("loc_labr==%i",det));
		TH1F *h = (TH1F*) gROOT->FindObject(Form("hDet%d",det));

		TSpectrum *s = new TSpectrum(10);
		s->Search(h);
		if (s->GetNPeaks() < 3) continue;

		//Sort the peaks
		long long ind[s->GetNPeaks()];
		TMath::Sort((long long)s->GetNPeaks(),s->GetPositionX(),ind);
		long long indTallest[s->GetNPeaks()];
		TMath::Sort((long long)s->GetNPeaks(),s->GetPositionY(),indTallest);

		//The first two peaks are the high energy Co60 identified by the alrgest x psoitions.
		int finalIndex[3];
		finalIndex[2] = ind[0];
		finalIndex[1] = ind[1];
		//The final peak is the Cs137 identified by the highest x position of the two highest peaks.
		finalIndex[0] = indTallest[0];
		for (int i=1;i<3;i++) {
			if (s->GetPositionX()[indTallest[i]] > s->GetPositionX()[finalIndex[0]] &&
				indTallest[i] != finalIndex[1]) 
				finalIndex[0] = indTallest[i];
		}

		

		//Fit these peaks
		TGraphErrors * gr = new TGraphErrors(3);
		Double_t correctEn[3] =    {0.661657,1.173228,1.332492};
		Double_t correctEnUnc[3] = {0.000003,0.000003,0.000004};
		float res, resUnc;
		for (int peak=2;peak >=0; peak--) {
			float pos = s->GetPositionX()[finalIndex[peak]];
			TF1 *fit = new TF1(Form("fit%d_%d",det,peak),"gaus");
			fit->SetRange(pos - 0.015 * pos, pos + 0.015 * pos);
			fit->SetParameter(1,pos);
			h->Fit(fit,"QMER+");

			gr->SetPoint(peak,fit->GetParameter(1),correctEn[peak]);
			gr->SetPointError(peak,fit->GetParameter(2),correctEnUnc[peak]);

			if (peak == 0) {
				res = fit->GetParameter(2) / fit->GetParameter(1) * 2 *sqrt(2 * log(2));
				resUnc = res * sqrt(pow(fit->GetParError(2)/fit->GetParameter(2),2) + pow(fit->GetParError(1) / fit->GetParameter(1),22)) ;
			}

		}
		c->Update();
		c->cd(pad)->cd(2);
		TF1 *fit = new TF1(Form("calibFit%d",det),"pol1");
		gr->Draw("AP");
		gr->Fit(fit,"QME");
		c->Update();
	
		printf("Det %d Max Gamma En %f +/- %f MeV Res %f +/- %f FWHM %%\n",	
			det, 
			fit->GetParameter(0) + fit->GetParameter(1) * 4095, 
			sqrt(pow(fit->GetParError(0),2) + pow(4095 * fit->GetParError(1),2)), 
			res*100, resUnc*100
		);
		
	}


		
}
