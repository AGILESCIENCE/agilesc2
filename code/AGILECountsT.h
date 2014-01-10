/***************************************************************************
                          AGILECountsT.h  -  description
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

#ifndef _AGILECOUNTST_H
#define _AGILECOUNTST_H

#include "CalibUtils.h"
#include "GenmapParams.h"
#include "EVTFilter.h"

class AGILECountsT {
public:
	AGILECountsT(string archiveevt);
	~AGILECountsT();
	bool EvalCounts(double tstart, double tstop, GammaExtractParams& params, uint32_t* resultingCts);
	bool prequery(double tmin, double tmax, GammaExtractParams& params);
	
	bool WritePhotonList(string out);
	
protected:
	string archiveevt;
	EVTFilter* evtfilter;
};

#endif
