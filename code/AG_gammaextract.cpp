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



class AppScreen
{
public:
    AppScreen()
    {
        cout << " "<< endl;
        cout << " "<< endl;
        cout << "#################################################################"<< endl;
        cout << "### GammaExtract V1.0 - 30/12/2013 - A.B., A.T., A.C., T.C. ###"<< endl;
        cout << "#################################################################"<< endl;
        cout << "#################################################################"<< endl;
        cout << endl;
    }

    ~AppScreen()
    {
        cout << endl << endl << endl;
        cout << "#################################################################"<< endl;
        cout << "##########  Task GammaExtract.......... exiting ###############"<< endl;
        cout << "#################################################################"<< endl << endl;
    }
};



static void PrintVector(const VecF& arr)
{
    cout << "[" << arr.Size() << "] ";
    for (int i=0; i<arr.Size(); ++i)
        cout << " " << arr[i];
    cout << endl;
}


int main(int argc,char **argv)
{


    GammaExtractParams params;
    if (!params.Load(argc, argv))
        return -1;
	
    Intervals intvs;
    double tmin = params["tmin"];
    double tmax = params["tmax"];
    const char* intFileName = params["timelist"];

    if (strcmp(intFileName, "None")) {
        intvs = ReadIntervals(intFileName);
        tmin = intvs.Min();
        tmax = intvs.Max();
        params["tmin"] = tmin;
        params["tmax"] = tmax;
    }
    else {
        Interval intv(tmin, tmax);
        intvs.Add(intv);
    }

    params.Print();
    if (intvs.Count()>1) {
        cout << intvs.Count() << " intervals:" << endl;
        for (int i=0; i<intvs.Count(); ++i)
            cout << "   " << String(intvs[i]) << endl;
    }


    

	cout << "GammaExtract......................evaluating the exposure"<< endl;

	double emin = params["emin"];
	double emax = params["emax"];
	double index = params["index"];
	const char* sarFile =  params["sarFileName"];
	const char* outfile = params["outfile"];
	double deltaT = params["timeslot"];
	const char* logfile = params["logfile"];
	const char* evtfile = params["evtfile"];
	int timestep = params["timestep"];
	int filtercode = params["filtercode"];

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
	
	string outfileevents = params.outfile;
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
		timeSlot.Set(beginTime, endTime);
		Intervals intervalSlots = Intersection(intvs, timeSlot);
		if (intervalSlots.Count()) {

			//cout << "Selected intervals" << endl;
			for (int i=0; i<intervalSlots.Count(); ++i) {
				cout << "slot:   " << setprecision(15) << intervalSlots[i].Start() << " " << intervalSlots[i].Stop() << " (" << intervalSlots[i].Stop() - intervalSlots[i].Start() << ") " << endl;
					double exp = -1;
					
					if (expagile->EvalExposure(intervalSlots[i].Start(), intervalSlots[i].Stop(), params, &exp)) {
						totalExposure += exp;
					}
					
					uint32_t cts = -1;
					if(ctsagile->EvalCounts(intervalSlots[i].Start(), intervalSlots[i].Stop(), params, &cts)) {
						totalCounts += cts;
						ctsagile->WritePhotonList(outfileevents);
						
					}
					
					if(exp != -1 && cts != -1 ) {
						expText << setprecision(1);
						expText << beginTime << " " << endTime << " ";
						expText << setprecision(2);
						expText << cts << " ";
						expText << exp << endl;
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
