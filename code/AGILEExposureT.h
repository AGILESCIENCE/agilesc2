/***************************************************************************
                          AGILEExposureT.h  -  description
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

#ifndef _AGILEEXPOSURET_H
#define _AGILEEXPOSURET_H

#include <cstdlib>
#include <cstring>
#include <cmath>

#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>

#include "pil.h"

#include "wcshdr.h"
#include "wcsmath.h"
#include "wcstrig.h"
#include "sph.h"

#include "CalibUtils.h"
#include "GenmapParams.h"

#include "LOGFilter.h"



class AGILEExposureT {
public:
	AGILEExposureT(string archivelog, string sarFile, uint32_t timestep, double emin, double emax, double index);
	~AGILEExposureT();
	bool EvalExposure(double tstart, double tstop, GammaExtractParams& params, double* resultingExp);

protected:
	double Area(double xbin, double ybin, double theta, int projection);
	double Alikesinaa(double input);
	inline double AG_expmapgen_area(double xbin, double ybin, double theta, int projection);	
	
	///return the exposure	
	double Exposure(LOGFilter* filter,  GammaExtractParams& params);
	
	LOGFilter* logfilter;
	AeffGridAverage* raeff;
	string archivelog;
	
	const double c_binFactor = 0.0003046174197867085688996857673060958405;
	const double c_angleFactor = 0.0174532925199432954743716805978692718782;
};

#endif