#ifndef __NTUPLE_GENERATOR__
#define __NTUPLE_GENERATOR__

#include "TROOT.h"
#include <iostream>
#include <map>

//________________________________________________________________________________________
//
// TNtuple_generator
// Class that generates a ntuple with Q2,NU,PT2,Zh,Phi,n_data,n_acc,n_thr
// and save it to a file,
//

class TNtuple_generator{
private:

	const TString f_location;   	 // The location of data files with ntuples inside
	const TString fd_ext;
	const TString fs_ext;
	const TString f_save;

	std::map<int, double> RC_CONSTS;

	const Double_t fQ2_min;			// low limit of the Q2 bin : 1
	const Double_t fQ2_max; 		// high limit of the Q2 bin : 4
	const Int_t fN_Q2;  			// Number of Q2 bins : 12
	const Double_t fNU_min; 		// low limit of the Nu bin : 2.2
	const Double_t fNU_max; 		// high limit of the Nu bin : 4.2
	const Int_t fN_NU; 				// Number of Nu bins : 9
	const Double_t fZH_max; 		// low limit of the Zh bin : 0.4
	const Double_t fZH_min; 		// high limit of the Zh bin : 0.7
	const Int_t fN_ZH; 				// Number of Zh bins: 3
	const Double_t fPT2_max; 		// low limit of the theta bin (in respect to the virtual photon) : -0
	const Double_t fPT2_min; 		// high limit of the theta bin (in respect to the virtual photon) : 4
	const Int_t fN_PT2; 			// Number of theta bins (in respect to the virtual photon) : 20
	const Double_t fPHI_PQ_max; 	// low limit of the phi bin (in respect to the virtual photon) : -180
	const Double_t fPHI_PQ_min;	  	// high limit of the phi bin (in respect to the virtual photon) : 180
	const Int_t fN_PHI_PQ; 			// Number of phi bins (in respect to the virtual photon) : 36

public:

	TNtuple_generator();
	void Read();
	void Save_ntuple();
	Int_t Return_index(Int_t, Int_t, Int_t, Int_t,Int_t) const;
	Double_t Transverse_momentum_broadening();
	void MakeCorrection();
	void ResizeHist();

};
inline Int_t TNtuple_generator::Return_index(Int_t i1, Int_t i2, Int_t i3, Int_t i4, Int_t i5) const
	{
		// Inline function to calculate index of RC_CONSTS map
		return (i1+1) + (i2+1)*10 + (i3+1)*100 + (i4+1)*10000 + i5*1000000;
	}
#endif  //__NTUPLE_GENERATOR__
