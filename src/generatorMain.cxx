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
	gROOT->SetStyle("Plain");
	TApplication *rootapp = new TApplication("TNtuple generation",&argc,argv);
	TNtuple_generator *gen = new TNtuple_generator();
	gen->save_ntuple();
	rootapp->Run();
	return 0;
}

#endif
