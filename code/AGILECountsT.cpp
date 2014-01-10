/***************************************************************************
                          AGILECountsT.cpp  -  description
                             -------------------
    copyright            : (C) 2014 Andrea Bulgarelli, Tomaso Contessi, Andrew Chen
    email                : bulgarelli@iasfbo.inaf.it
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software for non commercial purpose              *
 *   and for public research institutes; you can redistribute it and/or    *
 *   modify it under the terms of the GNU General Public License.          *
 *   For commercial purpose see appropriate license terms                  *
 *                                                                         *
 ***************************************************************************/
 
 
#include "AGILECountsT.h"
#include <fstream>

AGILECountsT::AGILECountsT(string archiveevt) {
	evtfilter = new EVTFilter(archiveevt); 
	this->archiveevt = archiveevt;
}

AGILECountsT::~AGILECountsT() {
	delete evtfilter;
}

bool AGILECountsT::prequery(double tstart, double tstop,  GammaExtractParams& params) {
	double mdim = params["mdim"];
	double la = params["la"];
	double ba = params["ba"];
	
	evtfilter->setPostfilter1(mdim, la, ba);
	
	int phasecode = params["phasecode"];
	int filtercode =  params["filtercode"];
	double emin = params["emin"];
	double emax = params["emax"];
	double albrad =  params["albrad"];
	double fovradmin =  params["fovradmin"];
	double fovradmax =  params["fovradmax"];
	
	if(evtfilter->prequery(tstart, tstop, phasecode, filtercode, emin, emax, albrad, fovradmin, fovradmax )) {
		;		
	} else 
		return false;
	return true;
}

bool AGILECountsT::EvalCounts(double tstart, double tstop,  GammaExtractParams& params, uint32_t *resultingCts) {
	
		
	
	double mdim = params["mdim"];
	double la = params["la"];
	double ba = params["ba"];
	
	evtfilter->setPostfilter1(mdim, la, ba);
	
	int phasecode = params["phasecode"];
	int filtercode =  params["filtercode"];
	double emin = params["emin"];
	double emax = params["emax"];
	double albrad =  params["albrad"];
	double fovradmin =  params["fovradmin"];
	double fovradmax =  params["fovradmax"];
	
	if(evtfilter->query(tstart, tstop, phasecode, filtercode, emin, emax, albrad, fovradmin, fovradmax )) {
		*resultingCts = evtfilter->time.size(); 
		//cout << "nrows cts: " << *resultingCts << endl;		
	} else 
		return false;
	return true;
}

bool AGILECountsT::WritePhotonList(string outfile) {
	ofstream mapText(outfile.c_str(), std::ofstream::app);
	streamsize prec = mapText.precision();
	
	for(long i=0; i< evtfilter->time.size(); i++) {
		mapText.setf(ios::fixed);
		mapText << setprecision(5) << evtfilter->time[i];
		mapText.unsetf(ios::floatfield);
		mapText.unsetf(ios::fixed);
		mapText <<" "<< setprecision(prec) <<" "<< evtfilter->ra[i] <<" "<< evtfilter->dec[i] <<" "<< (int) evtfilter->theta[i] <<" ";
		mapText << (int) evtfilter->energy[i] <<" "<< (int) evtfilter->ph_earth[i] <<" "<< (int) evtfilter->evstatus[i] << endl;						
	}	
	mapText.close();
}
