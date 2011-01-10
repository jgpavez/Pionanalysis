//______________________________________________________________________________________
//
// TNtuple_generator class implementation
//

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
#endif

TNtuple_generator::TNtuple_generator(): f_location("/home/jgpave/Pionanalysis/DATA/"), fd_ext("_data.root"), fs_ext("_simuls.root"),
										f_save("/home/utfsm/jgpavez/Pionanalysis/DATA/data_corr.root"),
										fQ2_min(1.), fQ2_max(4.), fN_Q2(12), fNU_min(2.2), fNU_max(4.2),
										fN_NU(9), fZH_min(0.4), fZH_max(0.7), fN_ZH(3),fPT2_max(4.),
										fPT2_min(0.),fN_PT2(20), fPHI_PQ_max(180.), fPHI_PQ_min(-180.), fN_PHI_PQ(36)
{

}

void TNtuple_generator::Read()
{
	// Read Function, write data from file to RC_CONSTS map.
	// the index in the map is computed by
	// Nu+10*Q2+100*Zh+10000*phi+1000000*pt2

	ifstream in;
  	in.open("test_C.dat");
  	int v1, v2, v3, v4, v5;
  	double v6;
  	while(in >> v1 >> v2 >> v3 >> v4 >> v5 >> v6) {
    	RC_CONSTS[v1 + 10*v2 + 100*v3 + 10000*v5 + 1000000*v4] = v6;
  	}
  	std::cout<<"map is created"<<std::endl;
}


void TNtuple_generator::Save_ntuple()
{
	// Function that generates a file with the ntuple of
	// Q2,Nu,Zh,Phi,Pt2,n_data,n_acc,n_thr and save it in
	// correct_data.root
	// TODO: Improve performance

	Read();
	TFile *data_corr = new TFile(f_save,"RECREATE");
	TNtuple *data_ntuple = new TNtuple("data_ntuple","data correction","target_data:Q2:NU:Zh:Phi:Pt2:n_acc:n_thr:n_data");
	Int_t i,j,k;
	TCut target_cut;
	Float_t target_data;
	const Double_t delta_q2 = (fQ2_max-fQ2_min)/fN_Q2;
  	const Double_t delta_nu = (fNU_max-fNU_min)/fN_NU;
  	const Double_t delta_Zh = (fZH_max-fZH_min)/fN_ZH;

	for ( i=0; i < fN_Q2; i++){
    	for (j=0;j < fN_NU;j++){
      		for ( k=0; k < fN_ZH; k++){
      			const Double_t Q2_min = fQ2_min+i*delta_q2;
      			const Double_t Q2_max = fQ2_min+(i+1)*delta_q2;
      			const Double_t Nu_min = fNU_min+j*delta_nu;
      			const Double_t Nu_max = fNU_min+(j+1)*delta_nu;
      			const Double_t Zh_min = fZH_min+k*delta_Zh;
      			const Double_t Zh_max = fZH_min+(k+1)*delta_Zh;

      			TCut Q2_cut = Form("Q2>%f && Q2<%f",Q2_min,Q2_max);
      			TCut Nu_cut = Form("Nu>%f && Nu<%f",Nu_min,Nu_max);
      			TCut Zh_cut = Form("Zh_pi_plus>%f && Zh_pi_plus<%f",Zh_min,Zh_max);
      			for ( Int_t metal = 0; metal < 2 ; metal++ ){
      				if (metal == 0){
      					target_cut = Form("targ_type==%d",2);
      				}
      				else {
      					target_cut = Form("targ_type==%d",1);
      				}
      			  	TCut cuts=Q2_cut && Nu_cut && Zh_cut && target_cut;
      			  	TCut cuts_simul=Q2_cut&&Nu_cut&&Zh_cut;
      			  	TCut Phi_pq_cut;

      			  	TChain *ntuple = new TChain("ntuple_pion");
      			  	ntuple->Add(f_location+"C"+fd_ext);
      			  	ntuple->Draw(">>list",cuts,"goffentrylist");
      			  	ntuple->SetEntryList((TEntryList*)gDirectory->Get("list"));
      			  	TChain *accept = new TChain("pi_accepted");
          		  	if (metal == 0) accept->Add(f_location+"C"+fs_ext);
      				else 	accept->Add(f_location+"D2"+fs_ext);
      			  	accept->Draw(">>list_acc",cuts_simul,"goffentrylist");
      			  	accept->SetEntryList((TEntryList*)gDirectory->Get("list_acc"));
      			  	TChain *thrown = new TChain("pi_thrown");
      			  	if (metal == 0) thrown->Add(f_location+"C"+fs_ext);
      				else  thrown->Add(f_location+"D2"+fs_ext);
      			 	thrown->Draw(">>list_thr",cuts_simul,"goffentrylist");
      			  	thrown->SetEntryList((TEntryList*)gDirectory->Get("list_thr"));

      			  	const Double_t delta_phi_pq = (fPHI_PQ_max-fPHI_PQ_min)/fN_PHI_PQ;

      				for(Int_t ii=0;ii<fN_PHI_PQ;ii++){
      			    	Phi_pq_cut = Form("phi_pq>%f && phi_pq<%f",fPHI_PQ_min+ii*delta_phi_pq,fPHI_PQ_min+(ii+1)*delta_phi_pq);
      			    	ntuple->Draw((const char*)Form("Pt2_pi_plus>>htmp_data(%d,%f,%f)",fN_PT2,fPT2_min,fPT2_max),Phi_pq_cut,"goff");
      			    	TH1F *htmp_data = (TH1F*)gDirectory->GetList()->FindObject("htmp_data");
      			    //	htmp_data->Sumw2();
      			    	accept->Draw((const char*)Form("Pt2_pi_plus>>htmp_acc(%d,%f,%f)",fN_PT2,fPT2_min,fPT2_max),Phi_pq_cut,"goff");
      			    	TH1F *htmp_acc = (TH1F*)gDirectory->GetList()->FindObject("htmp_acc");
      			    //	htmp_acc->Sumw2();
      			    	thrown->Draw((const char*)Form("Pt2_pi_plus>>htmp_thr(%d,%f,%f)",fN_PT2,fPT2_min,fPT2_max),Phi_pq_cut,"goff");
      			    	TH1F *htmp_thr = (TH1F*)gDirectory->GetList()->FindObject("htmp_thr");
      			    //	htmp_thr->Sumw2();

      					for (Int_t kk=1; kk<=htmp_data->GetXaxis()->GetNbins(); kk++) {
      						if (metal == 0) target_data = 2.0;
      						else target_data = 1.0;
      						data_ntuple->Fill(target_data,Q2_min+(delta_q2/2),Nu_min+(delta_nu/2),Zh_min+(delta_Zh/2),fPHI_PQ_min+(ii*delta_phi_pq)+(delta_phi_pq/2),fPT2_min+(kk*((fPT2_max-fPT2_min)/fN_PT2))+(((fPT2_max-fPT2_min)/fN_PT2)/2),htmp_acc->GetBinContent(kk),htmp_thr->GetBinContent(kk),htmp_data->GetBinContent(kk));
      			    	}

      			  	}

      			  	delete ntuple;
      			  	delete accept;
      			  	delete thrown;
      				}

      		}
    	}
  	}
	data_corr->cd();
	data_ntuple->Write();
	data_corr->Close();

}

void TNtuple_generator::MakeCorrection()
{
	// Make accceptance correction over ntuple generated with save_ntuple
	// acceptance correction = (n_data*n_thr)/n_acc
	// save a new file, corrected_data.root with the data corrected

	TFile *data = new TFile(f_save);
	TFile *correctData = new TFile(f_location + "acceptance_correction.root","RECREATE");

	TNtuple *tuple = (TNtuple *)data->Get("data_ntuple");
	TNtuple *correctionTuple = new TNtuple("corrected_data","acceptance correction applied to data","target_data:Q2:NU:Zh:Pth:Pt2:n_data");

	Float_t q2,nu,zh,phi,pt2,n_acc,n_thr,n_data,correction,target_data;

	tuple->SetBranchAddress("target_data",&target_data);
	tuple->SetBranchAddress("Q2",&q2);
	tuple->SetBranchAddress("NU",&nu);
	tuple->SetBranchAddress("Zh",&zh);
	tuple->SetBranchAddress("Phi",&phi);
	tuple->SetBranchAddress("Pt2",&pt2);
	tuple->SetBranchAddress("n_acc",&n_acc);
	tuple->SetBranchAddress("n_thr",&n_thr);
	tuple->SetBranchAddress("n_data",&n_data);

	for ( Int_t i = 0; i <= tuple->GetEntries(); i++ ){
		tuple->GetEntry(i);
		if ( n_acc == 0 ) continue;
		correction = (n_data*n_thr)/n_acc;
		correctionTuple->Fill(target_data,q2,nu,zh,phi,pt2,correction);
		}

	correctData->cd();
	correctionTuple->Write();
	correctData->Close();
	data->Close();
	delete correctData,correctData = NULL;
	delete data, data = NULL;
}

void TNtuple_generator::Transverse_momentum_broadening()
{
	// Tranverse momentum broadening computation
	// use data from acceptance_correction.root generated
	// by MakeCorrection

	TFile *data = new TFile(f_location + "acceptance_correction.root");
	TFile *broad = new TFile(f_location + "broadening2.root","RECREATE");
	TNtuple *data_ntuple = (TNtuple *)data->Get("corrected_data");
	TNtuple *broadening_t = new TNtuple("broadening","broadening calculated","broadening:Q2:NU:ZH");

	Double_t meanD2,meanC;
	Double_t broadening;
	TCut pt2Cut;
	TH1D *proycx,*proydx;
	Double_t delta_q2,delta_nu,delta_zh;
	delta_q2 = (fQ2_max-fQ2_min)/3;
	delta_nu = (fNU_max-fNU_min)/3;
	delta_zh = (fZH_max-fZH_min)/3;

	for ( Int_t i = 0; i < 3; i++){
		for ( Int_t k = 0; k < 3; k++){
			for ( Int_t j = 0; j < 3; j++){
				TCut cutc = Form("target_data == 2.0 && (Q2 > %f) && (Q2 < %f) && (NU > %f) && (NU < %f) && (Zh > %f) && (Zh < %f)",
									fQ2_min+i*delta_q2,fQ2_min+(i+1)*delta_q2,fNU_min+k*delta_nu,fNU_min+(k+1)*delta_nu,fZH_min+j*delta_zh ,fZH_min+(j+1)*delta_zh);
				TCut cutd2 = Form("target_data == 1.0 && (Q2 > %f) && (Q2 < %f) && (NU > %f) && (NU < %f) && (Zh > %f) && (Zh < %f)",
									fQ2_min+i*delta_q2,fQ2_min+(i+1)*delta_q2,fNU_min+k*delta_nu,fNU_min+(k+1)*delta_nu,fZH_min+j*delta_zh ,fZH_min+(j+1)*delta_zh);

				data_ntuple->Draw(Form("n_data:Pt2>>htmpc(%d,%f,%f)",fN_PT2,fPT2_min,fPT2_max),cutc,"profilegoff");
				data_ntuple->Draw(Form("n_data:Pt2>>htmpd2(%d,%f,%f)",fN_PT2,fPT2_min,fPT2_max),cutd2,"profilegoff");

				proycx = (TH1D*)((TProfile*)gDirectory->Get("htmpc"))->ProjectionX("htmpcp");
				proydx = (TH1D*)((TProfile*)gDirectory->Get("htmpd2"))->ProjectionX("htmpd2p");

				meanC = proycx->GetMean();
				meanD2 = proydx->GetMean();
				broadening = meanC - meanD2;

				broadening_t->Fill(broadening,i+1,k+1,j+1);
				delete gDirectory->Get("htmpcp");
				delete gDirectory->Get("htmpd2p");
			}
		}
	}
	broad->cd();
	broadening_t->Write();
	broad->Close(); delete broad;
	data->Close();  delete data;
}
//_____________________________________________________________________________________________________________
