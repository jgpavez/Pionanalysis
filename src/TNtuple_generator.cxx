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
#include "TMath.h"
#include <cmath>
#endif

TNtuple_generator::TNtuple_generator(): f_location("/home/jgpave/Pionanalysis/DATA/"), fd_ext("_data.root"), fs_ext("_simuls.root"),
										f_save("/home/jgpave/Pionanalysis/DATA/data_corr_new.root"),
										fQ2_min(1.), fQ2_max(4.), fN_Q2(4), fNU_min(2.2), fNU_max(4.2),
										fN_NU(3), fZH_min(0.4), fZH_max(0.7), fN_ZH(3),fPT2_max(4.),
										fPT2_min(0.),fN_PT2(5), fPHI_PQ_max(180.), fPHI_PQ_min(-180.), fN_PHI_PQ(12)
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
	TNtuple *data_ntuple = new TNtuple("data_ntuple","data correction","target_data:Q2:NU:X:Zh:Phi:Pt2:Pt:n_acc:n_thr:n_data");
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
					Float_t xb = (Q2_min+(delta_q2/2))/(2*(Nu_min+(delta_nu/2))*0.938272013);

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
      							data_ntuple->Fill(target_data,Q2_min+(delta_q2/2),Nu_min+(delta_nu/2), xb, Zh_min+(delta_Zh/2),fPHI_PQ_min+(ii*delta_phi_pq)+(delta_phi_pq/2),htmp_acc->GetBinCenter(kk), TMath::Sqrt(htmp_acc->GetBinCenter(kk)), htmp_acc->GetBinContent(kk),htmp_thr->GetBinContent(kk),htmp_data->GetBinContent(kk));
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
	TFile *correctData = new TFile(f_location + "acceptance_correction_new.root","RECREATE");

	TNtuple *tuple = (TNtuple *)data->Get("data_ntuple");
	TNtuple *correctionTuple = new TNtuple("corrected_data","acceptance correction applied to data","target_data:Q2:NU:X:Zh:Phi:Pt2:Pt:n_data:data_error");

	Float_t q2,nu,xb,zh,phi,pt2,pt,n_acc,n_thr,n_data,correction,target_data;

	tuple->SetBranchAddress("target_data",&target_data);
	tuple->SetBranchAddress("Q2",&q2);
	tuple->SetBranchAddress("NU",&nu);
	tuple->SetBranchAddress("X",&xb);
	tuple->SetBranchAddress("Zh",&zh);
	tuple->SetBranchAddress("Phi",&phi);
	tuple->SetBranchAddress("Pt2",&pt2);
	tuple->SetBranchAddress("Pt",&pt);
	tuple->SetBranchAddress("n_acc",&n_acc);
	tuple->SetBranchAddress("n_thr",&n_thr);
	tuple->SetBranchAddress("n_data",&n_data);

	for ( Int_t i = 0; i <= tuple->GetEntries(); i++ ){
		tuple->GetEntry(i);
		if ( n_acc == 0 || n_thr == 0 ) continue;
		Float_t e2_data = TMath::Abs(n_data);
		Float_t e_data = TMath::Sqrt(e2_data);
		Float_t e2_thr = TMath::Abs(n_thr);
		Float_t e_thr = TMath::Sqrt(e2_thr);
		Float_t e2_acc = TMath::Abs(n_acc);
		Float_t e_acc = TMath::Sqrt(e2_acc);

		//Divide acc/thr

		Float_t w1 = n_acc/n_thr;
		//Binomial Error
		Float_t b221 = n_thr*n_thr;
		Float_t err1;
		if ( n_acc != n_thr )
			err1 = TMath::Abs( ( (1.-2.*w1)*e_acc*e_acc + w1*w1*e_thr*e_thr )/(b221) );
		else
			err1 = 0;

	//	std::cout<<w1<<std::endl;
		//Divide data/w1
		correction = n_data/w1;
		//Normal error
		Float_t b222 = w1*w1;
		Float_t err2 = (e_data*e_data*w1*w1 + err1*err1*n_data*n_data)/(b222*b222);
		std::cout<<err2<<std::endl;
		correctionTuple->Fill(target_data,q2,nu,xb,zh,phi,pt2,pt,correction,err2);
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

void TNtuple_generator::calculateBins(){
	TFile *acceptance = new TFile(f_location + "acceptance_correction.root");
	TNtuple *old_ntuple = (TNtuple *)acceptance->Get("corrected_data");
	TFile *dataxb = new TFile(f_location + "dataxb.root","RECREATE");
	TNtuple *new_ntuple = new TNtuple("ntuple","ntuple","Q2:NU:Zh:Pt2:X:Phi:n_data");

	const Double_t delta_q2 = (fQ2_max-fQ2_min)/fN_Q2;
  	const Double_t delta_nu = (fNU_max-fNU_min)/fN_NU;
  	const Double_t delta_Zh = (fZH_max-fZH_min)/fN_ZH;
  	const Double_t delta_phi = (fPHI_PQ_max-fPHI_PQ_min)/fN_PHI_PQ;
  	Float_t xb;
	for ( int i = 0; i < fN_Q2; i++){
		for ( int j = 0; j < fN_NU; j++){
			for ( int k = 0; k < fN_ZH; k++){
				for ( int p = 0; p < fN_PHI_PQ; p++){
					const Double_t Q2_min = fQ2_min+i*delta_q2;
					const Double_t Nu_min = fNU_min+j*delta_nu;
					const Double_t Zh_min = fZH_min+k*delta_Zh;
					const Double_t Phi_min = fPHI_PQ_min + p*delta_phi;
					TCut cut1 = Form("Q2 > %f && Q2 < %f && NU > %f && NU < %f && Zh > %f && Zh < %f && Pth > %f && Pth < %f",
							Q2_min, Q2_min+(delta_q2),Nu_min, Nu_min+(delta_nu),Zh_min, Zh_min+(delta_Zh),Phi_min, Phi_min+(delta_phi));
					xb = (Q2_min+(delta_q2/2))/(2*(Nu_min+(delta_nu/2))*0.938272013);
					old_ntuple->Draw(Form("n_data:Pt2>>pt2hist(%d,%f,%f)",fN_PT2,fPT2_min,fPT2_max),cut1,"profilegoff");
					TH1D* proycx = (TH1D*)((TProfile*)gDirectory->Get("pt2hist"))->ProjectionX("htmpcp");

					for ( int pt = 1; pt <= 5; pt++){
						new_ntuple->Fill(Q2_min+(delta_q2/2),Nu_min+(delta_nu/2),Zh_min+(delta_Zh/2),proycx->GetBinCenter(pt),xb,Phi_min+(delta_phi/2),proycx->GetBinContent(pt));
					}
					delete proycx;
					delete gDirectory->Get("pt2hist");
				}
			}
		}
	}
	dataxb->cd();
	new_ntuple->Write();
	dataxb->Close();
	acceptance->Close();
	delete acceptance;
	delete dataxb;
}

void TNtuple_generator::createPhiPlots(){
	TFile *dataxb = new TFile(f_location + "acceptance_correction_new.root");
	TNtuple *tuple = (TNtuple *)dataxb->Get("corrected_data");
	TFile *phi_hists = new TFile(f_location + "phi_hists.root","RECREATE");
	phi_hists->cd();
	const Double_t XB_MIN = (fQ2_min/(2*fNU_min*0.938272013));
	const Double_t XB_MAX = (fQ2_max/(2*fNU_max*0.938272013));
	const Double_t delta_q2 = (fQ2_max-fQ2_min)/4;
  	const Double_t delta_xb = (XB_MAX - XB_MIN)/3;
  	const Double_t delta_Zh = (fZH_max-fZH_min)/fN_ZH;
  	const Double_t delta_pt2 = (fPT2_max - fPT2_min)/5;
	for ( int i = 0; i < 4; i++){
		for ( int j = 0; j < 3; j++){
			for ( int k = 0; k < fN_ZH; k++){
				for ( int p = 0; p < 5; p++){
					const Double_t Q2_min = fQ2_min+i*delta_q2;
					const Double_t Xb_min = XB_MIN+j*delta_xb;
					const Double_t Zh_min = fZH_min+k*delta_Zh;
					const Double_t Pt2_min = fPT2_min + p*delta_pt2;
					TCut cut1 = Form("Q2 > %f && Q2 < %f && X > %f && X < %f && Zh > %f && Zh < %f && Pt2 > %f && Pt2 < %f",
							Q2_min, Q2_min+(delta_q2),Xb_min, Xb_min+(delta_xb),Zh_min, Zh_min+(delta_Zh),Pt2_min, Pt2_min+(delta_pt2));
					tuple->Draw(Form("data_error:Phi>>errhist_%d_%d_%d_%d(%d,%f,%f)",i,j,k,p,12,fPHI_PQ_min,fPHI_PQ_max),cut1,"goffprofile");
					tuple->Draw(Form("n_data:Phi>>pt2hist_%d_%d_%d_%d(%d,%f,%f)",i,j,k,p,12,fPHI_PQ_min,fPHI_PQ_max),cut1,"goffprofile");
					TH1D* proycx = (TH1D*)((TProfile*)gDirectory->Get(Form("pt2hist_%d_%d_%d_%d",i,j,k,p)))->ProjectionX(Form("htmpcp_%d_%d_%d_%d",i,j,k,p));
					TH1D* proyerr = (TH1D*)((TProfile*)gDirectory->Get(Form("errhist_%d_%d_%d_%d",i,j,k,p)))->ProjectionX(Form("htmperr_%d_%d_%d_%d",i,j,k,p));
					for ( int e = 1; e < 12; e++)
						proycx->SetBinError(e,proyerr->GetBinContent(e));
					proycx->Write();
					delete proycx;
					delete proyerr;
					delete gDirectory->Get(Form("pt2hist_%d_%d_%d_%d",i,j,k,p));
					delete gDirectory->Get(Form("errhist_%d_%d_%d_%d",i,j,k,p));
				}
			}
		}
	}

	dataxb->Close();
	phi_hists->Close();
	delete phi_hists;
	delete dataxb;
}
//_____________________________________________________________________________________________________________
