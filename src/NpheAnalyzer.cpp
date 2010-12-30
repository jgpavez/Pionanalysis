/*
 * NpheAnalyzer.c
 *
 *  Created on: Dec 28, 2010
 *      Author: jgpavez
 */

#ifndef __CINT__
#include "TIdentificator.h"
#include "TClasTool.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TROOT.h"
#endif
//#include "Analyser/analysis_lib/include/massConst.h"
void calculateNPhe()
{
#ifdef __CINT__
	gROOT->Reset();
	gSystem->Load("libClasTool.so");
	gSystem->Load("libTIdentificator.so");
#endif

	TClasTool *input = new TClasTool();
	input->InitDSTReader("ROOTDSTR");
	input->AddFile("/home/utfsm/jgpavez/Pionanalysis/DATA/clas_42066_22.pass2.root");
	Double_t nPhe;

	TIdentificator *ident = new TIdentificator(input);

	TH1F *npheHist = new TH1F("nphe","nphe if firts particle is electron",1000,0.0,500);

	Int_t nEntries = input->GetEntries();

	for ( Int_t k = 0; k < nEntries; k++ ){
		input->Next();
		if ( ident->GetCategorization(0) == "electron" ){
			Int_t nRows = input->GetNRows("EVNT");
			for ( Int_t i = 0; i < nRows; i++ ){
				nPhe = ident->Nphe(i);
				if ( nPhe > 0.0 )
					npheHist->Fill(nPhe);
			}
		}
	}
	delete ident; ident = NULL;
	delete input; input = NULL;
}

void somePlots(){
#ifdef __CINT__
	gROOT->Reset();
	gSystem->Load("libClasTool.so");
	gSystem->Load("libTIdentificator.so");
#endif

	TClasTool *input = new TClasTool();
	input->InitDSTReader("ROOTDSTR");
	input->AddFile("/home/utfsm/jgpavez/Pionanalysis/DATA/clas_42066_22.pass2.root");
	Double_t nPhe;
	TH2D *hist = new TH2D("TimeCorr4","Histogram for time 4 correction vs Momentum",80,0,2.5,200,-20,20);
	TIdentificator *ident = new TIdentificator(input);

	Int_t nEntries = input->GetEntries();
	Double_t correction;
	for ( Int_t k = 0; k < nEntries; k++ ){
		input->Next();
				Int_t nRows = input->GetNRows("EVNT");
				if ( (ident->GetCategorization(0) == "electron") && ((ident->Q2() > 1.) && (ident->W() > 2.) && (ident->Yb() < 0.85)) && (ident->FidCheckCut() == 1)){
					for ( Int_t i = 0; i < nRows; i++ ){
						if ( (ident->GetCategorization(i) == "low energy pion +") && (ident->FidCheckCutPiPlus(i) == 1)){
							correction = ident->TimeCorr4(0.13957018,i);
							if ( correction < 1)
								hist->Fill(ident->Momentum(i),correction);
						}
					}
				}
	}
	delete ident; ident = NULL;
	delete input; input = NULL;

}
