#ifndef __CINT__

#include "TROOT.h"
#include "Riostream.h"
#include "TApplication.h"
#include "TNtuple_generator.h"
#include "TCanvas.h"

extern void InitGui();
VoidFuncPtr_t initfuncs[] = {InitGui,0};

TROOT my_app("myapp","myapp",initfuncs);

int main(int argc, char **argv)
{
	if (argc == 1){
		std::cout<<"Seleccione funcion:"<<std::endl;
		std::cout<<"1) save_ntuple"<<std::endl;
		std::cout<<"2) make_correction"<<std::endl;
		std::cout<<"3) calculate tranverse_momentum_broadening"<<std::endl;
		return 0;
	}
	gROOT->SetStyle("Plain");
	TApplication *rootapp = new TApplication("TNtuple generation",&argc,argv);
	TNtuple_generator *gen = new TNtuple_generator();
	switch (atoi(argv[1])){
		case(1):gen->Transverse_momentum_broadening();
		break;
		case(2):gen->MakeCorrection();
		break;
		case(3):gen->Transverse_momentum_broadening();
		break;
		}
//	rootapp->Run();
	delete rootapp;
	delete gen;
	return 0;
}

#endif
