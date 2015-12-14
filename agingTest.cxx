#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <math.h>
#include <string.h>

//Program for fiting 137Cs peak at 661.637 keV

#define NUM_FILES 44
#define RES_MAX 50

#define RES_HI 5.0
#define RES_LO 1.25

#define FIT_WIDTH_PERCENT 0.08
void stonehenge(){

int i, j;
int maxCount[7] = {0,0,0,0,0,0,0};
int maxBin[7] = {-1,-1,-1,-1,-1,-1,-1};

double centroid[7][NUM_FILES];
double sigma[7][NUM_FILES];
double amplitude[7][NUM_FILES];
double resolution[7][NUM_FILES];


TGraph *res_v_time0 = new TGraph(NUM_FILES);
TGraph *chan_v_time0 = new TGraph(NUM_FILES);
TGraph *res_v_time1 = new TGraph(NUM_FILES);
TGraph *chan_v_time1 = new TGraph(NUM_FILES);
TGraph *res_v_time2 = new TGraph(NUM_FILES);
TGraph *chan_v_time2 = new TGraph(NUM_FILES);
TGraph *res_v_time3 = new TGraph(NUM_FILES);
TGraph *chan_v_time3 = new TGraph(NUM_FILES);
TGraph *res_v_time4 = new TGraph(NUM_FILES);
TGraph *chan_v_time4 = new TGraph(NUM_FILES);
TGraph *res_v_time5 = new TGraph(NUM_FILES);
TGraph *chan_v_time5 = new TGraph(NUM_FILES);
TGraph *res_v_time6 = new TGraph(NUM_FILES);
TGraph *chan_v_time6 = new TGraph(NUM_FILES);

TMultiGraph *chan_v_time = new TMultiGraph();
TMultiGraph *res_v_time = new TMultiGraph();

TCanvas *c1 = new TCanvas("c1","c1",200,10,800,800);


//List of times
int time[NUM_FILES]{
			0,    3600, 7200, 10800,14400,18000,21600,25200,28800,32400,
			36000,39600,43200,46800,50400,54000,57600,61200,64800,68400,
			72000,75600,79200,82800,86400,90000,93600,97200,100800,104400,
			108000,111600,115200,118800,122400,126000,129600,133200,136800,140400,
			144000,147600,151200,154800
/*
			,158400,162000,165600,169200,172800,176400,
			180000,183600,187200,190800,194400,198000
*/
		};

//List of corresponding file names
char* filename[NUM_FILES]{	
				"Res_0s.root.root",
				"Res_3600s.root.root",
				"Res_7200s.root.root",
				"Res_10800s.root.root",
				"Res_14400s.root.root",
				"Res_18000s.root.root",
				"Res_21600s.root.root",
				"Res_25200s.root.root",
				"Res_28800s.root.root",
				"Res_32400s.root.root",

				"Res_36000s.root.root",
				"Res_39600s.root.root",
				"Res_43200s.root.root",
				"Res_46800s.root.root",
				"Res_50400s.root.root",
				"Res_54000s.root.root",
				"Res_57600s.root.root",
				"Res_61200s.root.root",
				"Res_64800s.root.root",
				"Res_68400s.root.root",

				"Res_72000s.root.root",
				"Res_75600s.root.root",
				"Res_79200s.root.root",
				"Res_82800s.root.root",
				"Res_86400s.root.root",
				"Res_90000s.root.root",
				"Res_93600s.root.root",
				"Res_97200s.root.root",
				"Res_100800s.root.root",
				"Res_104400s.root.root",

				"Res_108000s.root.root",
				"Res_111600s.root.root",
				"Res_115200s.root.root",
				"Res_118800s.root.root",
				"Res_122400s.root.root",
				"Res_126000s.root.root",
				"Res_129600s.root.root",
				"Res_133200s.root.root",
				"Res_136800s.root.root",
				"Res_140400s.root.root",

				"Res_144000s.root.root",
				"Res_147600s.root.root",
				"Res_151200s.root.root",
				"Res_154800s.root.root"
/*
				"Res_158400s.root.root",
				"Res_162000s.root.root",
				"Res_165600s.root.root",
				"Res_169200s.root.root",
				"Res_172800s.root.root",
				"Res_176400s.root.root",

				"Res_180000s.root.root",
				"Res_183600s.root.root",
				"Res_187200s.root.root",
				"Res_190800s.root.root",
				"Res_194400s.root.root",
				"Res_198000s.root.root"
*/
			};

printf("\n\n");

c1->Divide(2,4,0,0);

for(j=0; j<NUM_FILES; j++){		//Loop over files and fit peak


	printf("File: %d \n",j);
	//Load Root File with 137Cs spectrum
	TFile *file = new TFile(filename[j]);
	//Get Historgram from file (h000 or h001) 
	TH1F * h0 = new TH1F("h0","h000 title");
	TH1F * h1 = new TH1F("h1","h001 title");
	TH1F * h2 = new TH1F("h2","h002 title");
	TH1F * h3 = new TH1F("h3","h003 title");
	TH1F * h4 = new TH1F("h4","h004 title");
	TH1F * h5 = new TH1F("h5","h005 title");
	TH1F * h6 = new TH1F("h6","h006 title");

	h0 = (TH1F*)file->Get("h000");
	h1 = (TH1F*)file->Get("h001");
	h2 = (TH1F*)file->Get("h002");
	h3 = (TH1F*)file->Get("h003");
	h4 = (TH1F*)file->Get("h004");
	h5 = (TH1F*)file->Get("h005");
	h6 = (TH1F*)file->Get("h006");

	//Find Peak (Bin with Max counts)
	//only loop over region where peak is expected 
	for(i=7500; i<20000; i++){
		if(h0->GetBinContent(i)>maxCount[0] && i > 150){
			maxCount[0] = h0->GetBinContent(i);
			maxBin[0] = i;
		}

		if(h1->GetBinContent(i)>maxCount[1] && i > 150){
			maxCount[1] = h1->GetBinContent(i);
			maxBin[1] = i;
		}
		if(h2->GetBinContent(i)>maxCount[2] && i > 150){
			maxCount[2] = h2->GetBinContent(i);
			maxBin[2] = i;
		}
		if(h3->GetBinContent(i)>maxCount[3] && i > 150){
			maxCount[3] = h3->GetBinContent(i);
			maxBin[3] = i;
		}
		if(h4->GetBinContent(i)>maxCount[4] && i > 150){
			maxCount[4] = h4->GetBinContent(i);
			maxBin[4] = i;
		}
		if(h5->GetBinContent(i)>maxCount[5] && i > 150){
			maxCount[5] = h5->GetBinContent(i);
			maxBin[5] = i;
		}
		if(h6->GetBinContent(i)>maxCount[6] && i > 150){
			maxCount[6] = h6->GetBinContent(i);
			maxBin[6] = i;
		}

	}
	//Creeate Peak and fit 
	//Fit region is maxBin +/- 1000 bins

	TF1 *pk0 = new TF1("pk0","gaus",maxBin[0]-(FIT_WIDTH_PERCENT*maxBin[0]),maxBin[0]+(FIT_WIDTH_PERCENT*maxBin[0]));
	TF1 *pk1 = new TF1("pk1","gaus",maxBin[1]-(FIT_WIDTH_PERCENT*maxBin[1]),maxBin[1]+(FIT_WIDTH_PERCENT*maxBin[1]));
	TF1 *pk2 = new TF1("pk2","gaus",maxBin[2]-(FIT_WIDTH_PERCENT*maxBin[2]),maxBin[2]+(FIT_WIDTH_PERCENT*maxBin[2]));
	TF1 *pk3 = new TF1("pk3","gaus",maxBin[3]-(FIT_WIDTH_PERCENT*maxBin[3]),maxBin[3]+(FIT_WIDTH_PERCENT*maxBin[3]));
	TF1 *pk4 = new TF1("pk4","gaus",maxBin[4]-(FIT_WIDTH_PERCENT*maxBin[4]),maxBin[4]+(FIT_WIDTH_PERCENT*maxBin[4]));
	TF1 *pk5 = new TF1("pk5","gaus",maxBin[5]-(FIT_WIDTH_PERCENT*maxBin[5]),maxBin[5]+(FIT_WIDTH_PERCENT*maxBin[5]));
	TF1 *pk6 = new TF1("pk6","gaus",maxBin[6]-(FIT_WIDTH_PERCENT*maxBin[6]),maxBin[6]+(FIT_WIDTH_PERCENT*maxBin[6]));


	h0->Fit(pk0,"MQR","same",maxBin[0]-(FIT_WIDTH_PERCENT*maxBin[0]),maxBin[0]+(FIT_WIDTH_PERCENT*maxBin[0]));
	h1->Fit(pk1,"MQR","same",maxBin[1]-(FIT_WIDTH_PERCENT*maxBin[1]),maxBin[1]+(FIT_WIDTH_PERCENT*maxBin[1]));
	h2->Fit(pk2,"MQR","same",maxBin[2]-(FIT_WIDTH_PERCENT*maxBin[2]),maxBin[2]+(FIT_WIDTH_PERCENT*maxBin[2]));
	h3->Fit(pk3,"MQR","same",maxBin[3]-(FIT_WIDTH_PERCENT*maxBin[3]),maxBin[3]+(FIT_WIDTH_PERCENT*maxBin[3]));
	h4->Fit(pk4,"MQR","same",maxBin[4]-(FIT_WIDTH_PERCENT*maxBin[4]),maxBin[4]+(FIT_WIDTH_PERCENT*maxBin[4]));
	h5->Fit(pk5,"MQR","same",maxBin[5]-(FIT_WIDTH_PERCENT*maxBin[5]),maxBin[5]+(FIT_WIDTH_PERCENT*maxBin[5]));
	h6->Fit(pk6,"MQR","same",maxBin[6]-(FIT_WIDTH_PERCENT*maxBin[6]),maxBin[6]+(FIT_WIDTH_PERCENT*maxBin[6]));

	//Store Fit Parameters in Arrays

	amplitude[0][j] = pk0->GetParameter(0);
	centroid[0][j] = pk0->GetParameter(1);
	sigma[0][j] = pk0->GetParameter(2);
	resolution[0][j] = (sigma[0][j]*2.3548)/(centroid[0][j])*100;

	amplitude[1][j] = pk1->GetParameter(0);
	centroid[1][j] = pk1->GetParameter(1);
	sigma[1][j] = pk1->GetParameter(2);
	resolution[1][j] = (sigma[1][j]*2.3548)/(centroid[1][j])*100;

	amplitude[2][j] = pk2->GetParameter(0);
	centroid[2][j] = pk2->GetParameter(1);
	sigma[2][j] = pk2->GetParameter(2);
	resolution[2][j] = (sigma[2][j]*2.3548)/(centroid[2][j])*100;

	amplitude[3][j] = pk3->GetParameter(0);
	centroid[3][j] = pk3->GetParameter(1);
	sigma[3][j] = pk3->GetParameter(2);
	resolution[3][j] = (sigma[3][j]*2.3548)/(centroid[3][j])*100;

	amplitude[4][j] = pk4->GetParameter(0);
	centroid[4][j] = pk4->GetParameter(1);
	sigma[4][j] = pk4->GetParameter(2);
	resolution[4][j] = (sigma[4][j]*2.3548)/(centroid[4][j])*100;

	amplitude[5][j] = pk5->GetParameter(0);
	centroid[5][j] = pk5->GetParameter(1);
	sigma[5][j] = pk5->GetParameter(2);
	resolution[5][j] = (sigma[5][j]*2.3548)/(centroid[5][j])*100;

	amplitude[6][j] = pk6->GetParameter(0);
	centroid[6][j] = pk6->GetParameter(1);
	sigma[6][j] = pk6->GetParameter(2);
	resolution[6][j] = (sigma[6][j]*2.3548)/(centroid[6][j])*100;

	//Only plot fits within in resolution range
	//Not really needed anymore but still left it in

	for(i=0; i<7; i++){
		if( (resolution[i][j] < RES_HI) && (resolution[i][j] > RES_LO) ){
			if( i == 0 ){
				res_v_time0->SetPoint(j,time[j]/60/60,resolution[i][j]);
				chan_v_time0->SetPoint(j,time[j]/60/60,centroid[i][j]);
			}else if( i == 1){
				res_v_time1->SetPoint(j,time[j]/60/60,resolution[i][j]);
				chan_v_time1->SetPoint(j,time[j]/60/60,centroid[i][j]);
			}else if( i == 2){
				res_v_time2->SetPoint(j,time[j]/60/60,resolution[i][j]);
				chan_v_time2->SetPoint(j,time[j]/60/60,centroid[i][j]);
			}else if( i == 3){
				res_v_time3->SetPoint(j,time[j]/60/60,resolution[i][j]);
				chan_v_time3->SetPoint(j,time[j]/60/60,centroid[i][j]);
			}else if( i == 4){
				res_v_time4->SetPoint(j,time[j]/60/60,resolution[i][j]);
				chan_v_time4->SetPoint(j,time[j]/60/60,centroid[i][j]);
			}else if( i == 5){
				res_v_time5->SetPoint(j,time[j]/60/60,resolution[i][j]);
				chan_v_time5->SetPoint(j,time[j]/60/60,centroid[i][j]);
			}else if( i == 6){
				res_v_time6->SetPoint(j,time[j]/60/60,resolution[i][j]);
				chan_v_time6->SetPoint(j,time[j]/60/60,centroid[i][j]);
			}
		}
	}

	c1->cd(1);
	h0->Draw();
	c1->cd(2);
	h1->Draw();
	c1->cd(3);
	h2->Draw();
	c1->cd(4);
	h3->Draw();
	c1->cd(5);
	h4->Draw();
	c1->cd(6);
	h5->Draw();
	c1->cd(7);
	h6->Draw();
	

	c1->Update();
//	sleep(1);

} // End of File loop

c1->Clear();
c1->Divide(1,2,0,0);
//////////////////// Panel 1 ///////////////////////////////////////
c1->cd(1);
gPad->SetGrid();

chan_v_time->SetTitle("Channel vs. Time  (661.657 keV)");
chan_v_time0->SetMarkerColor(kRed);
chan_v_time1->SetMarkerColor(kMagenta);
chan_v_time2->SetMarkerColor(kBlue);
chan_v_time3->SetMarkerColor(kCyan);
chan_v_time4->SetMarkerColor(kGreen);
chan_v_time5->SetMarkerColor(kYellow);
chan_v_time6->SetMarkerColor(kBlack);

chan_v_time->Add(chan_v_time0,"p*");
chan_v_time->Add(chan_v_time1,"p*");
chan_v_time->Add(chan_v_time2,"p*");
chan_v_time->Add(chan_v_time3,"p*");
chan_v_time->Add(chan_v_time4,"p*");
//chan_v_time->Add(chan_v_time5,"p*");
chan_v_time->Add(chan_v_time6,"p*");
chan_v_time->Draw("a*");



//////////////////// Panel 2 ///////////////////////////////////////
c1->cd(2);
gPad->SetGrid();

res_v_time->SetTitle("Resolution vs. Time  (661.657 keV)");
res_v_time0->SetMarkerColor(kRed);
res_v_time1->SetMarkerColor(kMagenta);
res_v_time2->SetMarkerColor(kBlue);
res_v_time3->SetMarkerColor(kCyan);
res_v_time4->SetMarkerColor(kGreen);
res_v_time5->SetMarkerColor(kYellow);
res_v_time6->SetMarkerColor(kBlack);

res_v_time->Add(res_v_time0,"p*");
res_v_time->Add(res_v_time1,"p*");
res_v_time->Add(res_v_time2,"p*");
res_v_time->Add(res_v_time3,"p*");
res_v_time->Add(res_v_time4,"p*");
//res_v_time->Add(res_v_time5,"p*");
res_v_time->Add(res_v_time6,"p*");
res_v_time->Draw("a*");

c1->Update();
printf("\n\n");


} // End of Main























