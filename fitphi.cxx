// -------------------------------------------------------------------------
// SIMPLE SCRIPT FOR PHI FIT

#ifndef __CINT__
#include <map>
#include "TROOT.h"
#include "TTree.h"
#include "TNtuple_generator.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TCut.h"
#include "TH1F.h"
#include "TCanvas.h"
#include <iostream>
#include <fstream>
#include <cstdio>
#include "TChain.h"
#include "TProfile.h"
#include <cmath>
#include "RooFit.h"
#endif


const Int_t fN_Q2 = 12;  			// Number of Q2 bins : 12
const Int_t fN_XB = 3; 				// Number of Nu bins : 9		// low limit of the Zh bin : 0.4
const Int_t fN_ZH = 3; 				// Number of Zh bins: 3
const Int_t fN_PT2 = 5; 			// Number of theta bins (in respect to the virtual photon) : 20
const Int_t fN_PHI_PQ = 12;

void fitphi()
{
	TFile *phihist = new TFile("~/Pionanalysis/DATA/phi_hists.root");
	TFile *phihistfit = new TFile("~/Pioanalysis/DATA/phi_hist_fit.root","RECREATE");
	for (int i = 0; i < fN_Q2; i++){
		for ( int j = 0 ; j < fN_XB ; j++){
			for ( int k = 0; k < fN_ZH; k++){
				for ( int p = 0; p < fN_PT2; p++){
					phihist->cd();
					TH1D *phiplot = (TH1D *)phihist->Get(Form("pt2hist_%d_%d_%d_%d",i,j,k,p));
					//A+B*cos(phi)+C*cos(2*phi)
					RooRealVar phi("x","x",0,-4,4);
					RooRealVar A("A","A",0,-4,4);
					RooRealVar B("B","B",0,-4,4);
					RooRealVar C("C","C",0,-4.4);

					RooPlot *phiframe = phi.Frame();
					phiframe->SetName(Form("phihist_%d_%d",i,j));
					//ptframe->SetTitle(frame_title);

		            RooDataHist dh("dh", "dh", phi, Import(*phiplot));
					dh.plotOn(phiframe,DataError(RooAbsData::SumW2));

					const char *phi_formula = "A + B*cos(phi) + C*cos(2*phi)";
					RooGenericPdf formula("formula","phi formula",phi_formula,
				            RooArgSet(phi,A,B,C));
					formula.plotOn(phiframe);
					formula.fitTo(dm,SumW2Error(kTRUE),Save());

					formula.plotOn(phiframe);
					std::cout<<phiframe->chiSquare()<<std::endl;
					sum_chi += phiframe->chiSquare();
					phihistfit->cd();
					phiframe->Write(Form("phihist_%d_%d",i,j));
					delete phiplot;

				}
			}
		}
	}

}
