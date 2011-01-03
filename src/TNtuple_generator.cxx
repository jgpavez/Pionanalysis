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

TNtuple_generator::TNtuple_generator(): f_location("/home/utfsm/jgpavez/Pionanalysis/DATA/"), fd_ext("_data.root"), fs_ext("_simuls.root"), f_save("/home/utfsm/jgpavez/Pionanalysis/DATA/data_corr.root"), fQ2_min(1.), fQ2_max(4.), fN_Q2(12), fNU_min(2.2), fNU_max(4.2), fN_NU(9), fZH_min(0.4), fZH_max(0.7), fN_ZH(3),fPT2_max(4.),
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
	data->Close();
	correctData->cd();
	correctionTuple->Write();
	correctData->Close();
	delete correctData;
	delete data;
	data = NULL;
	correctData = NULL;
}

void TNtuple_generator::Transverse_momentum_broadening()
{
	// Tranverse momentum broadening computation
	// use data from acceptance_correction.root generated
	// by MakeCorrection

	TFile *data = new TFile(f_location + "acceptance_correction.root");
	TFile *resized = new TFile(f_location + "resized_hist.root");
	TFile *broad = new TFile(f_location + "broadening.root","RECREATE");
	TNtuple *data_ntuple = (TNtuple *)data->Get("corrected_data");
	TNtuple *broadening_t = new TNtuple("broadening","broadening calculated","broadening");
	Double_t meanD2,meanC;
	Float_t broadening;
	TProfile *resizedQ2 = (TProfile *)resized->Get("resizeQ2");
	TProfile *resizedNu = (TProfile *)resized->Get("resizeNu");
	TProfile *resizedZh = (TProfile *)resized->Get("resizeZh");
	TCut pt2Cut;
	for ( Int_t i = 1; i <= 3; i++){
		for ( Int_t k = 1; k <= 3; k++){
			for ( Int_t j = 1; j <= 3; j++){
				TCut cutc = Form("target_data == 2.0 && (Q2 > %f) && (Q2 < %f) && (NU > %f) && (NU < %f) && (Zh > %f) && (Zh < %f)",resizedQ2->GetBinLowEdge(i),resizedQ2->GetBinLowEdge(i) + resizedQ2->GetBinWidth(i),resizedNu->GetBinLowEdge(k) ,resizedNu->GetBinLowEdge(k) + resizedNu->GetBinWidth(k) ,resizedZh->GetBinLowEdge(j) ,resizedZh->GetBinLowEdge(j) + resizedZh->GetBinWidth(j));
				TCut cutd2 = Form("target_data == 1.0 && (Q2 > %f) && (Q2 < %f) && (NU > %f) && (NU < %f) && (Zh > %f) && (Zh < %f)",resizedQ2->GetBinLowEdge(i),resizedQ2->GetBinLowEdge(i) + resizedQ2->GetBinWidth(i),resizedNu->GetBinLowEdge(k) ,resizedNu->GetBinLowEdge(k) + resizedNu->GetBinWidth(k) ,resizedZh->GetBinLowEdge(j) ,resizedZh->GetBinLowEdge(j) + resizedZh->GetBinWidth(j));
				data_ntuple->Draw(Form("n_data:Pt2>>htmpc(%d,%f,%f)",fN_PT2,fPT2_min,fPT2_max),cutc,"profilegoff");
				data_ntuple->Draw(Form("n_data:Pt2>>htmpd2(%d,%f,%f)",fN_PT2,fPT2_min,fPT2_max),cutd2,"profilegoff");
				std::cout<<cutc<<std::endl;
				std::cout<<cutd2<<std::endl;
				((TProfile*)gDirectory->Get("htmpc"))->Sumw2();
				((TProfile*)gDirectory->Get("htmpd2"))->Sumw2();
				meanC = ((TProfile*)gDirectory->Get("htmpc"))->GetMean();
				meanD2 = ((TProfile*)gDirectory->Get("htmpd2"))->GetMean();
				broadening = (pow(meanC,2)) - (pow(meanD2,2));
				std::cout<<broadening<<std::endl;
				broadening_t->Fill(broadening);

			}
		}
	}
	broad->cd();
	broadening_t->Write();
	broad->Close(); delete broad;
	data->Close();  delete data;
}

void TNtuple_generator::ResizeHist()
{
	// Bin resizing for values in acceptance_correction.root
	// This recalculate bins, and generate histograms for Q2,NU,Zh
	// also save its in resized_hist.root

	TFile *data = new TFile(f_location + "acceptance_correction.root");
	TFile *resize = new TFile(f_location + "resized_hist.root","RECREATE");

	TNtuple *tuple = (TNtuple *)data->Get("corrected_data");
	tuple->Draw(Form("n_data:Q2>>profQ2(%d,%f,%f)",fN_Q2,fQ2_min,fQ2_max),"","profilegoff");
	TProfile* beforeresProfileQ2 = (TProfile*)gDirectory->Get("profQ2");
	beforeresProfileQ2->Rebin(4,"resizeQ2");
	TProfile *resizeProfileQ2 = (TProfile *)gDirectory->Get("resizeQ2");
	tuple->Draw(Form("n_data:NU>>profNu(%d,%f,%f)",fN_NU,fNU_min,fNU_max),"","profilegoff");
	TProfile* beforeresProfileNu = (TProfile*)gDirectory->Get("profNu");
	beforeresProfileNu->Rebin(3,"resizeNu");
	TProfile *resizeProfileNu = (TProfile *)gDirectory->Get("resizeNu");
	tuple->Draw(Form("n_data:Zh>>resizeZh(%d,%f,%f)",fN_ZH,fZH_min,fZH_max),"","profilegoff");
	TProfile *resizeProfileZh = (TProfile *)gDirectory->Get("resizeZh");

	resize->cd();
	resizeProfileQ2->Write();
	resizeProfileNu->Write();
	resizeProfileZh->Write();
	resize->Close();
	delete resize; resize = NULL;
	data->Close();
	delete data; data = NULL;
}

//_____________________________________________________________________________________________________________
