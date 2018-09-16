/***************************************************************************
                          main.cpp
                             -------------------
    copyright            : (C) 2013 AGILE Team
    email                : bulgarelli@iasfbo.inaf.it
    contributors		 : Andrew Chen, Alberto Pellizzoni, Alessio Trois (IASF-Milano),
    					   Andrea Bulgarelli (IASF-Bologna), Tomaso Contessi (Nuove Idee sas)

 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include <string.h>

//#define DEBUG 1

#include <FitsUtils.h>
#include <Selection.h>
#include <Eval.h>
#include <PilParams.h>

using std::cout;
using std::endl;

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

#include "AGILEExposureT.h"
#include "AGILECountsT.h"

const char* startString = {
	"################################################################\n"
	"###                   Task AG_apT v1.0.0 - A.B.              ###"
};

const char* endString = {
	"### Task AG_apT exiting .................................... ###\n"
	"################################################################"
};

const PilDescription paramsDescr[] = {
	{ PilString, "outfile", "Output file name" },
	{ PilString, "logfile", "Grid log index file name" },
	{ PilString, "evtfile", "Event file index file name" },
	{ PilString, "sarFileName", "Effective area file name" },
	{ PilString, "edpFileName", "Energy dispersion file name" },
	{ PilString, "timelist", "Time intervals file name" },
	{ PilReal, "mres", "Bin size (degrees)" },
	{ PilReal, "la", "Longitude of map center (galactic)" },
	{ PilReal, "ba", "Latitude of map center (galactic)" },
	{ PilReal, "lonpole", "Rotation of map (degrees)" },
	{ PilReal, "albrad", "Radius of earth albedo (degrees)" },
	{ PilReal, "y_tol", "Boresight movement tolerance (degrees)" },
	{ PilReal, "roll_tol", "Roll tolerance (degrees)" },
	{ PilReal, "earth_tol", "Earth tolerance (degrees)" },
	{ PilInt, "phasecode", "Orbital phase code" },
	{ PilInt, "timestep", "LOG file step size" },
	{ PilReal, "index", "Spectral index" },
	{ PilReal, "tmin", "Initial time (sec)" },
	{ PilReal, "tmax", "Final time (sec)" },
	{ PilReal, "emin", "Minimum energy" },
	{ PilReal, "emax", "Maximum energy" },
	{ PilReal, "fovradmin", "Min radius of field of view (degrees)" },
	{ PilReal, "fovradmax", "Max radius of field of view (degrees)" },
	{ PilInt, "filtercode", "Event filter code" },
	{ PilReal, "timeslot", "Time slot" },
	{ PilNone, "", "" }
};


using namespace std;

static string String(const Interval& intv)
{
    stringstream str(ios_base::out);
    str.precision(6);
    str << fixed << intv.Start() << ".." << intv.Stop();
    return str.str();
}

static string String(const Intervals& intvs)
{
    stringstream str(ios_base::out);
    const char* sep = "";
    for (int i=0; i<intvs.Count(); ++i) {
        str << sep << String(intvs[i]);
        sep = ", ";
    }
    return str.str();
}



struct timespec start, stop;
struct timespec startg, stopg;




static void PrintVector(const VecF& arr)
{
    cout << "[" << arr.Size() << "] ";
    for (int i=0; i<arr.Size(); ++i)
        cout << " " << arr[i];
    cout << endl;
}


int main(int argc,char **argv)
{
	cout << startString << endl;

	PilParams params(paramsDescr);
	if (!params.Load(argc, argv))
		return EXIT_FAILURE;

	cout << endl << "INPUT PARAMETERS:" << endl;
	params.Print();

	double emin = params["emin"];
	double emax = params["emax"];
	double index = params["index"];
	const char* sarFile =  params["sarFileName"];
	const char* outfile = params["outfile"];
	double deltaT = params["timeslot"];
	const char* logfile = params["logfile"];
	const char* evtfile = params["evtfile"];
	int filtercode = params["filtercode"];
	double y_tol = params["y_tol"];
	double earth_tol = params["earth_tol"];
	double albrad = params["albrad"];

	double la = params["la"];
	double ba = params["ba"];
	double lonpole = params["lonpole"];
	double fovradmin = params["fovradmin"];
	double fovradmax = params["fovradmax"];
	int phasecode = params["phasecode"];


	double tmin = params["tmin"];
	double tmax = params["tmax"];

	int timestep = params["timestep"];

	double mdim = params["mres"];
	mdim = mdim * sqrt(2);

	double radius = params["mres"];

	double binstep = 1.0;
	const char *projection = "ARC";
	cout << "radius for evt: " << radius << " - mdim for exp: " << mdim << endl;
	cout << "Binstep: " << binstep << endl;
	cout << "Projection: " << projection << endl;

	//area calculation
	double pixel1 = DEG2RAD * DEG2RAD * fabs(mdim * mdim);
	double areapixel =  pixel1 * Sinaa(DEG2RAD*45.);
	cout << "### Area pixel " << areapixel << endl;




	Intervals intervals;
	if (!eval::LoadTimeList(params["timelist"], intervals, tmin, tmax)) {
		cerr << "Error loading timelist file '" << params["timelist"].GetStr() << "'" << endl;
		return EXIT_FAILURE;
	}


	cout << "INTERVALS N=" << intervals.Count() << ":" << endl;
	for (int i=0; i<intervals.Count(); i++)
		cout << "   " << intervals[i].String() << endl;



	cout << "GammaExtract......................evaluating the exposure"<< endl;



	ofstream expText(outfile);
	expText.setf(ios::fixed);
	double beginTime = tmin;
	double endTime = beginTime+deltaT;
	if (endTime>tmax)
		endTime = tmax;
	cout << "***** " << setprecision(25) << beginTime << " " << endTime << " " << deltaT << endl << endl;
	double totalExposure = 0;
	long  totalCounts = 0;
	Interval timeSlot;
	int status = 0;


	cout << evtfile << endl;
	AGILECountsT* ctsagile = new AGILECountsT(evtfile);

	AGILEExposureT* expagile = new AGILEExposureT(logfile, sarFile, timestep, emin, emax, index);

	string outfileevents = outfile;
	outfileevents += ".photons";

	//pre-query here

	/*
	cout << "**** prequery start " << endl;
	if(!ctsagile->prequery(tmin, tmax, params))
		cout << "evt prequery problems " << endl;
	if(!expagile->prequery(tmin, tmax, params))
		cout << "log prequery problems " << endl;
	cout << "**** prequery ok " << endl;
	*/
	do {
#ifdef DEBUG
		cout << "Time slot beginTime: " << beginTime << " endTime: " << endTime << endl;
#endif
		timeSlot.Set(beginTime, endTime);
		Intervals intervalSlots = Intersection(intervals, timeSlot);
		if (intervalSlots.Count()) {
			cout << "Selected slots:" << endl;
			for (int i=0; i<intervalSlots.Count(); ++i) {
				cout << "slot:   " << setprecision(15) << intervalSlots[i].Start() << " " << intervalSlots[i].Stop() << " (" << intervalSlots[i].Stop() - intervalSlots[i].Start() << ")  count=" << intervalSlots.Count() <<endl;
				double exp = -1;

				if (expagile->EvalExposure(intervalSlots[i].Start(), intervalSlots[i].Stop(), &exp, mdim, la, ba, lonpole, fovradmin, fovradmax, emin, emax, index, y_tol, earth_tol, phasecode, albrad, timestep)) {
					totalExposure += exp;
					cout << "exp: " << exp << " total exp: "<< totalExposure << endl;

				}

				uint32_t cts = -1;
				if(ctsagile->EvalCounts(intervalSlots[i].Start(), intervalSlots[i].Stop(), params, &cts)) {
					totalCounts += cts;
					ctsagile->WritePhotonList(outfileevents);

				}
				cout << "totalCounts: " << totalCounts << endl;



				if(exp != -1 && cts != -1 ) {
					expText << setprecision(1);
					expText << beginTime << " " << endTime << " ";
					expText << setprecision(2);
					expText << exp/areapixel << " ";
					expText << cts << endl;

				} else {
					cerr << "problems in the query " << cts << " " << exp << endl;
				}
				if(exp == 0 && cts != 0)
					cout << "WARNING: " << beginTime << " " << endTime << " " << endl;
			}
		}
		//else
		//	cerr << "No intervals selected" << endl;
		beginTime = endTime;
		endTime += deltaT;
		if (tmax<endTime)
			endTime = tmax;
	} while (beginTime<tmax);

	expText.close();

	cout << "Total Exposure: " << totalExposure << endl;
	cout << "Total Counts: " << totalCounts << endl;

    if (status) {
        if (status != 105)
            remove(outfile);
        cout << "GammaExtract.................... exiting ERROR: " << status;
        if (status>0)
            fits_report_error(stdout, status);
        cout << endl;
    }
    else
        cout << "GammaExtract.................... exiting SUCCESS" << endl;

    return status;
}
