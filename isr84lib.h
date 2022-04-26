
#ifndef _ISR84LIB_H
#define _ISR84LIB_H

//*****************************************************************************************************
//*                                                                                                   *
//*  This code is free software; you can redistribute it and/or modify it at your will.               *
//*  It is my hope however that if you improve it in any way you will find a way to share it too.     *
//*                                                                                                   *
//*  If you have any comments, questions, suggestions etc. mail me (jgray77@gmail.com)                *
//*                                                                                                   *
//*  This program is distributed AS-IS in the hope that it will be useful, but WITHOUT ANY WARRANTY;  *
//*  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.        *
//*                                                                                                   *
//*****************************************************************************************************
//
//
//===================================================================================================
//	Israel Local Grids <==> WGS84 conversions
//===================================================================================================
//
// The Israel New Grid (ITM) is a Transverse Mercator projection of the GRS80 ellipsoid.
// The Israel Old Grid (ICS) is a Cassini-Soldner projection of the modified Clark 1880 ellipsoid.
//
// To convert from a local grid to WGS84 you first do a "UTM to Lat/Lon" conversion using the
// known formulas but with the local grid data (Central Meridian, Scale Factor and False
// Easting and Northing). This results in Lat/Long in the local ellipsoid coordinate system.
// Afterwards you do a Molodensky transformation from this ellipsoid to WGS84.
//
// To convert from WGS84 to a local grid you first do a Molodensky transformation from WGS84
// to the local ellipsoid, after which you do a Lat/Lon to UTM conversion, again with the data of
// the local grid instead of the UTM data.
//
// The UTM to Lat/Lon and Lat/Lon to UTM conversion formulas were taken as-is from the
// excellent article by Prof.Steven Dutch of the University of Wisconsin at Green Bay:
//		http://www.uwgb.edu/dutchs/UsefulData/UTMFormulas.htm
//
// The [abridged] Molodensky transformations were taken from
//		http://home.hiwaay.net/~taylorc/bookshelf/math-science/geodesy/datum/transform/molodensky/
// and can be found in many sources on the net.
//
// Additional sources:
// ===================
// 1. dX,dY,dZ values:  http://www.geo.hunter.cuny.edu/gis/docs/geographic_transformations.pdf
//
// 2. ITM data:  http://www.mapi.gov.il/geodesy/itm_ftp.txt
//    for the meridional arc false northing, the value is given at
//    http://www.mapi.gov.il/reg_inst/dir2b.doc
//    (this doc also gives a different formula for Lat/lon -> ITM, but not the reverse)
//
// 3. ICS data:  http://www.mapi.gov.il/geodesy/ics_ftp.txt
//    for the meridional arc false northing, the value is given at several places as the
//    correction value for Garmin GPS sets, the origin is unknown.
//    e.g. http://www.idobartana.com/etrexkb/etrexisr.htm
//
// Notes:
// ======
// 1. The conversions between ICS and ITM are
//			ITM Lat = ICS Lat - 500000
//			ITM Lon = ICS Lon + 50000
//	  e.g. ITM 678000,230000 <--> ICS 1178000 180000
//
//	  Since the formulas for ITM->WGS84 and ICS->WGS84 are different, the results will differ.
//    For the above coordinates we get the following results (WGS84)
//		ITM->WGS84 32.11'43.945" 35.18'58.782"
//		ICS->WGS84 32.11'43.873" 35.18'58.200"
//      Difference    ~3m            ~15m
//
// 2. If you have, or have seen, formulas that contain the term Sin(1"), I recommend you read
//    Prof.Dutch's enlightening explanation about it in his link above.
//
//===================================================================================================
#ifndef __LIBRARY_H__
#define __LIBRARY_H__

#include <stdio.h>
#include <math.h>

double pi() { return 3.141592653589793; }
double sin2(double x) { return sin(x)*sin(x); }
double cos2(double x) { return cos(x)*cos(x); }
double tan2(double x) { return tan(x)*tan(x); }
double tan4(double x) { return tan2(x)*tan2(x); }

enum eDatum {
    eWGS84=0,
    eGRS80=1,
    eCLARK80M=2
};

typedef struct _DTM {
    double a;	// a  Equatorial earth radius
    double b;	// b  Polar earth radius
    double f;	// f= (a-b)/a  Flatenning
    double esq;	// esq = 1-(b*b)/(a*a)  Eccentricity Squared
    double e;	// sqrt(esq)  Eccentricity
    // deltas to WGS84
    double dX;
    double dY;
    double dZ;
} DATUM,*PDATUM;

static DATUM Datum[3] = {

        // WGS84 data
        {
                6378137.0,				// a
                6356752.3142,			// b
                0.00335281066474748,  	// f = 1/298.257223563
                0.006694380004260807,	// esq
                0.0818191909289062, 	// e
                // deltas to WGS84
                0,
                0,
                0
        },

        // GRS80 data
        {
                6378137.0,				// a
                6356752.3141,			// b
                0.0033528106811823,		// f = 1/298.257222101
                0.00669438002290272,	// esq
                0.0818191910428276,		// e
                // deltas to WGS84
                -48,
                55,
                52
        },

        // Clark 1880 Modified data
        {
                6378300.789,			// a
                6356566.4116309,		// b
                0.003407549767264,		// f = 1/293.466
                0.006803488139112318,	// esq
                0.08248325975076590,	// e
                // deltas to WGS84
                -235,
                -85,
                264
        }
};

enum gGrid {
    gICS=0,
    gITM=1
};

typedef struct _GRD {
    double lon0;
    double lat0;
    double k0;
    double false_e;
    double false_n;
} GRID,*PGRID;

static GRID Grid[2] = {

        // ICS data
        {
                0.6145667421719,			// lon0 = central meridian in radians of 35.12'43.490"
                0.55386447682762762,		// lat0 = central latitude in radians of 31.44'02.749"
                1.00000,					// k0 = scale factor
                170251.555,					// false_easting
                2385259.0					// false_northing
        },

        // ITM data
        {
                0.61443473225468920,		// lon0 = central meridian in radians 35.12'16.261"
                0.55386965463774187,		// lat0 = central latitude in radians 31.44'03.817"
                1.0000067,					// k0 = scale factor
                219529.584,					// false_easting
                2885516.9488				// false_northing = 3512424.3388-626907.390
                // MAPI says the false northing is 626907.390, and in another place
                // that the meridional arc at the central latitude is 3512424.3388
        }
};

// prototypes of local Grid <=> Lat/Lon conversions
//void Grid2LatLon(int N,int E,double& lat,double& lon,gGrid from,eDatum to);
//void LatLon2Grid(double lat,double lon,int& N,int& E,eDatum from,gGrid to);
//// prototype of Moldensky transformation function
//void Molodensky(double ilat,double ilon,double& olat,double& olon,eDatum from,eDatum to);

//=================================================
// Israel New Grid (ITM) to WGS84 conversion
//=================================================
void itm2wgs84(int N,int E,double* lat,double *lon);

//=================================================
// WGS84 to Israel New Grid (ITM) conversion
//=================================================
void wgs842itm(double lat,double lon,int* N,int* E);

//=================================================
// Israel Old Grid (ICS) to WGS84 conversion
//=================================================
void ics2wgs84(int N,int E,double* lat,double *lon);

//=================================================
// WGS84 to Israel Old Grid (ICS) conversion
//=================================================
void wgs842ics(double lat,double lon,int* N,int* E);

//====================================
// Local Grid to Lat/Lon conversion
//====================================
void Grid2LatLon(int N,int E,double* lat,double* lon,enum gGrid from,enum eDatum to);

//====================================
// Lat/Lon to Local Grid conversion
//====================================
void LatLon2Grid(double lat,double lon,int* N,int* E,enum eDatum from,enum gGrid to);

//======================================================
// Abridged Molodensky transformation between 2 datums
//======================================================
static void Molodensky(double ilat,double ilon,double* olat,double* olon,enum eDatum from,enum eDatum to);


#endif

#endif //_ISR84LIB_H
