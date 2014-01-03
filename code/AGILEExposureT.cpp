/***************************************************************************
                          AGILEExposureT.cpp  -  description
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
 
 
#include "AGILEExposureT.h"

AGILEExposureT::AGILEExposureT(string archivelog, string sarFile, uint32_t timestep, double emin, double emax, double index) {
	logfilter = new LOGFilter(archivelog, timestep); 
	this->archivelog = archivelog;
	cout << "AeffGridAverage " << emin << " " << emax << " " << index << endl;
	raeff = new  AeffGridAverage(sarFile.c_str(), emin, emax, index);
}

AGILEExposureT::~AGILEExposureT() {
	delete logfilter;
}

double AGILEExposureT::Area(double xbin, double ybin, double theta, int projection)
{
    if (projection==ARC)
        return c_binFactor * xbin * ybin * Sinaa(c_angleFactor * theta);
    else /// AIT
        return c_binFactor * xbin * ybin;
}

double AGILEExposureT::Alikesinaa(double input){
	if (input == 0) return(1.0); else return(sin(input)/input);
}

inline double AGILEExposureT::AG_expmapgen_area(double xbin, double ybin, double theta, int projection) {
    if (projection==ARC)
        return 0.0003046174197867085688996857673060958405 * xbin * ybin * Alikesinaa(0.0174532925199432954743716805978692718782 * theta);
    else
        return 0.0003046174197867085688996857673060958405 * xbin * ybin; 
}

bool AGILEExposureT::EvalExposure(double tstart, double tstop, GammaExtractParams& params, double *resultingExp) {
	logfilter->reset();
	int phasecode = params["phasecode"];
	
	if(logfilter->query( tstart, tstop, phasecode )) {
		//cout << "nrows: " << logfilter->time.size() << endl;
		*resultingExp = Exposure(logfilter,  params);
	} else 
		return false;
	return true;
}



double AGILEExposureT::Exposure(LOGFilter* filter, GammaExtractParams& params)
{
    double lp = 0, bp = 0;
    double learth, bearth;
    double lp0 = 0, bp0 = 0;
    double x = 0, y = 0;
    double theta = 0, phi = 0, phi2 = 0;
    double eul[5];


    double theta2;
    double lng, lat;
    double area = 0;
    int ait;

    struct prjprm prj;
    prjini(&prj);

    double mdim = params["mdim"];
    //1.0 because mres=1 (one bin for each map)
    double mres = 1.0;
        
    switch (params.projection) {
    case ARC:
    	//OLD x = -(params.mdim/2)+(params.mres*(i+0.5));
    	//where i=0 and mres=1
    	//=> x = -(mdim/2) + (1 * 0.5)
    	//x = -mdim;
    	x = -(mdim/2)+(1*(0+0.5));
    	//y = -mdim;
    	y = -(mdim/2)+(1*(0+0.5));
        ait = 0;
        theta2 = 90.0-sqrt(x*x+y*y);
        phi2 = atan2d(-y, -x);
        eul[0] = params["la"];
        eul[1] = 90.0-double(params["ba"]);
        eul[2] = params["lonpole"];
        eul[3] = cosd(eul[1]);
        eul[4] = sind(eul[1]);

        sphx2s(eul, 1, 1, 0, 0, &phi2, &theta2, &lng, &lat);

        area = AG_expmapgen_area(mres, mres, 90-theta2, params.projection);

        /*cout << " x: " << x;
		cout << " y: " << y;
		cout << " theta2: " << theta2;
		cout << " phi2: " << theta2;
		cout << " lng " << lng;
		cout << " lat " << lat;
		for(int kii=0; kii<5; kii++) cout << " eul[" << kii << "] " << eul[kii];
		cout << setprecision(15) << " area " << area;
		cout << endl;
		*/
        break;
    case AIT:
        prj.r0 = RAD2DEG;
        prj.phi0 = params["la"];
        prj.theta0 = params["ba"];
        aitset(&prj);
        x = mdim;
        y = -mdim;
        aitx2s(&prj, 1, 1, 0, 0, &x, &y, &lng, &lat, &ait);
        area = AG_expmapgen_area(mres, mres, 0, params.projection);
    }

    long n = 0;
    double time = 0;

    long allnrows = filter->time.size();

    long* change = new long[allnrows];


    double A = 0.0; /// The resulting value
    
    double timestep = params["timestep"];

    //for (long nrows=0; nrows < allnrows; nrows++) {
        long count = 0;
        double earth_ra0 = filter->earth_ra[0], earth_dec0 = filter->earth_dec[0];
        double ra_y0 = filter->ra_y[0], dec_y0 = filter->dec_y[0];
        for (long k = 1; k<allnrows; ++k) {
        	/*
        	cout << "0\t" << k-1 << "\t" << filter->phase[k-1] << "\t" << filter->ra_y[k-1] << "\t" << filter->dec_y[k-1] << "\t" << filter->earth_ra[k-1] << "\t" << filter->earth_dec[k-1] << endl;
            cout << "YTolTest " << params.YTolTest(filter->ra_y[k-1], filter->dec_y[k-1], ra_y0, dec_y0) << endl;
            //cout << params.RollTolTest(&psi[k-1]) << endl;
            cout << "EarthTolTest " << params.EarthTolTest(filter->earth_ra[k-1], filter->earth_dec[k-1], earth_ra0, earth_dec0) << endl;
            */
        	if ((filter->phase[k-1] != filter->phase[k])
                    || params.YTolTest(filter->ra_y[k-1], filter->dec_y[k-1], ra_y0, dec_y0)
                    || params.EarthTolTest(filter->earth_ra[k-1], filter->earth_dec[k-1], earth_ra0, earth_dec0)) {
                //|| params.RollTolTest(&psi[k-1])
                change[count++] = k;
                earth_ra0 = filter->earth_ra[k-1];
                earth_dec0 = filter->earth_dec[k-1];
                ra_y0 = filter->ra_y[k-1];
                dec_y0 = filter->dec_y[k-1];
                /*
                cout << "----" << endl;
				cout << " count " << count;
				cout << " change[count] " << change[count];
				cout << " earth_ra0 " << earth_ra0;
				cout << " earth_dec0 " << earth_dec0;
				cout << " ra_y0 " << ra_y0;
				cout << " dec_y0 " << dec_y0;
				cout << endl;
				*/
            }
        }
        
        //cout << "----" << endl;
    	//cout << " count " << count << endl;
    	
        if (count == 0)
            change[0] = allnrows;

        long lowrow = 0;
        long highrow = 0;

        for (long k = 0; k<=count; ++k) {
            time = 0;
            if ( k == 0 ) {
                lowrow = 0;
                highrow = change[0];
            }
            else if ( k == count) {
                lowrow = change[k-1];
                highrow = allnrows;
            }
            else {
                lowrow = change[k-1];
                highrow = change[k];
            }

            for (n = lowrow; n<highrow; n++)
                time += filter->livetime[n] * timestep;
                
            //cout << "lowrow " << lowrow << " highrow " << highrow << " time " << time << endl;

            Euler(filter->ra_y[lowrow], filter->dec_y[lowrow], &lp, &bp, 1);
            Euler(filter->earth_ra[lowrow], filter->earth_dec[lowrow], &learth, &bearth, 1);
            
            if (ait == 0) {
                theta = SphDistDeg(lng, lat, lp, bp);
                phi = 0;
                /*cout << "theta " << theta << " " << lng << " " << lat << " " << lp << " " << bp << endl;
                cout << " learth " << learth;
                cout << " bearth " << bearth; 
                cout << " area " << area;
                cout << " time " << time;
                cout << endl;
                cout << params.FovTest(theta) << " ";
                cout << params.AlbTest(lng, lat, learth, bearth);
                cout << endl;
                */
                if (params.FovTest(theta) && params.AlbTest(lng, lat, learth, bearth))
                    A += 1e-3*time*(raeff->AvgVal(theta, phi))*area;
                //cout << " raeff.AvgVal(theta, phi) " << raeff.AvgVal(theta, phi) << endl;
                //cout << " A " << A << endl;
            }
        }
    //}

    delete[] change;
    return A;
}
