/*
 * NpheAnalyzer.c
 *
 *  Created on: Dec 28, 2010
 *      Author: jgpavez
 */

#ifndef __CINT__
#include "TIdentificator.h"
#include "TClastool.h"
#include "TH2F.h"
#include "TProfile.h"
#endif
void calculateNPhe()
{
#ifdef __CINT__
	gRoot->Reset();
	gSystem->Load("libClasTool.so");
	gSystem->Load("libTIdentificator.so");
#endif

	TClasTool
}
