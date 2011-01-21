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

	const TString f_location;   			 // The location of data files with ntuples inside
	const TString fd_ext;
	const TString fs_ext;
	const TString f_save;

	std::map<int, double> RC_CONSTS;

	static const double fQ2_min 	 = 	1;			// low limit of the Q2 bin : 1
	static const double fQ2_max 	 = 	4; 			// high limit of the Q2 bin : 4
	static const unsigned int 	 fN_Q2		 =  12;  		// Number of Q2 bins : 12
	static const double fNU_min 	 = 	2.2; 		// low limit of the Nu bin : 2.2
	static const double fNU_max 	 =  4.2; 		// high limit of the Nu bin : 4.2
	static const unsigned int 	 fN_NU		 =  9; 			// Number of Nu bins : 9
	static const double fZH_min 	 =  0.4;		// low limit of the Zh bin : 0.4
	static const double fZH_max 	 =  0.7; 		// high limit of the Zh bin : 0.7
	static const unsigned int 	 fN_ZH		 =  3; 		// Number of Zh bins: 3
	static const double fPT2_min	 =  0.; 		// low limit of the theta bin (in respect to the virtual photon) : -0
	static const double fPT2_max	 =  4.; 		// high limit of the theta bin (in respect to the virtual photon) : 4
	static const unsigned int	 fN_PT2  	 =  20; 		// Number of theta bins (in respect to the virtual photon) : 20
    static const double fPHI_PQ_min = -180.; 		// low limit of the phi bin (in respect to the virtual photon) : -180
    static const double fPHI_PQ_max =  180.;	  	// high limit of the phi bin (in respect to the virtual photon) : 180
    static const unsigned int	 fN_PHI_PQ	 =  36; 		// Number of phi bins (in respect to the virtual photon) : 36
    static const double fXB_min 	 = 	0.242225;
    static const double fXB_max 	 = 	0.507519;
    static const unsigned int	 fN_XB		 =  3;
    static const double fPT_min  	 =	0.;
    static const double fPT_max	 =  2.;
    static const int 	 fN_PT 		 =  5;

	void Read();

public:

	TNtuple_generator();
	void Save_ntuple();
	void Transverse_momentum_broadening() const;
	void MakeCorrection() const;
	void calculateBins() const;
	void RecalculateBins() const;
	Int_t Return_index(Int_t, Int_t, Int_t, Int_t,Int_t) const;

};
inline Int_t TNtuple_generator::Return_index(Int_t i1, Int_t i2, Int_t i3, Int_t i4, Int_t i5) const
	{
		// Inline function to calculate index of RC_CONSTS map
		return (i1+1) + (i2+1)*10 + (i3+1)*100 + (i4+1)*10000 + i5*1000000;
	}

#endif  //__NTUPLE_GENERATOR__
