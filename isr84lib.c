#include "isr84lib.h"


//=================================================
// Israel New Grid (ITM) to WGS84 conversion
//=================================================
void itm2wgs84(int N,int E,double* lat,double *lon)
{
// 1. Local Grid (ITM) -> GRS80
    double lat80,lon80;
    Grid2LatLon(N,E,&lat80,&lon80,gITM,eGRS80);

// 2. Molodensky GRS80->WGS84
    double lat84,lon84;
    Molodensky(lat80,lon80,&lat84,&lon84,eGRS80,eWGS84);

// final results
    *lat = lat84*180/pi();
    *lon = lon84*180/pi();
}

//=================================================
// WGS84 to Israel New Grid (ITM) conversion
//=================================================
void wgs842itm(double lat,double lon,int* N,int* E)
{
    double latr = lat*pi()/180;
    double lonr = lon*pi()/180;

// 1. Molodensky WGS84 -> GRS80
    double lat80,lon80;
    Molodensky(latr,lonr,&lat80,&lon80,eWGS84,eGRS80);

// 2. Lat/Lon (GRS80) -> Local Grid (ITM)
    LatLon2Grid(lat80,lon80,N,E,eGRS80,gITM);
}

//=================================================
// Israel Old Grid (ICS) to WGS84 conversion
//=================================================
void ics2wgs84(int N,int E,double* lat,double *lon)
{
// 1. Local Grid (ICS) -> Clark_1880_modified
    double lat80,lon80;
    Grid2LatLon(N,E,&lat80,&lon80,gICS,eCLARK80M);

// 2. Molodensky Clark_1880_modified -> WGS84
    double lat84,lon84;
    Molodensky(lat80,lon80,&lat84,&lon84,eCLARK80M,eWGS84);

// final results
    *lat = lat84*180/pi();
    *lon = lon84*180/pi();
}

//=================================================
// WGS84 to Israel Old Grid (ICS) conversion
//=================================================
void wgs842ics(double lat,double lon,int* N,int* E)
{
    double latr = lat*pi()/180;
    double lonr = lon*pi()/180;

// 1. Molodensky WGS84 -> Clark_1880_modified
    double lat80,lon80;
    Molodensky(latr,lonr,&lat80,&lon80,eWGS84,eCLARK80M);

// 2. Lat/Lon (Clark_1880_modified) -> Local Grid (ICS)
    LatLon2Grid(lat80,lon80,N,E,eCLARK80M,gICS);
}

//====================================
// Local Grid to Lat/Lon conversion
//====================================
void Grid2LatLon(int N,int E,double* lat,double* lon,enum gGrid from,enum eDatum to)
{
//================
// GRID -> Lat/Lon
//================

    double y = N + Grid[from].false_n;
    double x = E - Grid[from].false_e;
    double M = y / Grid[from].k0;

    double a = Datum[to].a;
    double b = Datum[to].b;
    double e = Datum[to].e;
    double esq = Datum[to].esq;

    double mu = M / (a*(1 - e*e/4 - 3*pow(e,4)/64 - 5*pow(e,6)/256));

    double ee = sqrt(1-esq);
    double e1 = (1-ee)/(1+ee);
    double j1 = 3*e1/2 - 27*e1*e1*e1/32;
    double j2 = 21*e1*e1/16 - 55*e1*e1*e1*e1/32;
    double j3 = 151*e1*e1*e1/96;
    double j4 = 1097*e1*e1*e1*e1/512;

// Footprint Latitude
    double fp =  mu + j1*sin(2*mu) + j2*sin(4*mu) + j3*sin(6*mu) + j4*sin(8*mu);

    double sinfp = sin(fp);
    double cosfp = cos(fp);
    double tanfp = sinfp/cosfp;
    double eg = (e*a/b);
    double eg2 = eg*eg;
    double C1 = eg2*cosfp*cosfp;
    double T1 = tanfp*tanfp;
    double R1 = a*(1-e*e) / pow(1-(e*sinfp)*(e*sinfp),1.5);
    double N1 = a / sqrt(1-(e*sinfp)*(e*sinfp));
    double D = x / (N1*Grid[from].k0);

    double Q1 = N1*tanfp/R1;
    double Q2 = D*D/2;
    double Q3 = (5 + 3*T1 + 10*C1 - 4*C1*C1 - 9*eg2*eg2)*(D*D*D*D)/24;
    double Q4 = (61 + 90*T1 + 298*C1 + 45*T1*T1 - 3*C1*C1 - 252*eg2*eg2)*(D*D*D*D*D*D)/720;
// result lat
    *lat = fp - Q1*(Q2-Q3+Q4);

    double Q5 = D;
    double Q6 = (1 + 2*T1 + C1)*(D*D*D)/6;
    double Q7 = (5 - 2*C1 + 28*T1 - 3*C1*C1 + 8*eg2*eg2 + 24*T1*T1)*(D*D*D*D*D)/120;
// result lon
    *lon = Grid[from].lon0 + (Q5 - Q6 + Q7)/cosfp;
}

//====================================
// Lat/Lon to Local Grid conversion
//====================================
void LatLon2Grid(double lat,double lon,int* N,int* E,enum eDatum from,enum gGrid to)
{
// Datum data for Lat/Lon to TM conversion
    double a = Datum[from].a;
    double e = Datum[from].e; 	// sqrt(esq);
    double b = Datum[from].b;

//===============
// Lat/Lon -> TM
//===============
    double slat1 = sin(lat);
    double clat1 = cos(lat);
    double clat1sq = clat1*clat1;
    double tanlat1sq = slat1*slat1 / clat1sq;
    double e2 = e*e;
    double e4 = e2*e2;
    double e6 = e4*e2;
    double eg = (e*a/b);
    double eg2 = eg*eg;

    double l1 = 1 - e2/4 - 3*e4/64 - 5*e6/256;
    double l2 = 3*e2/8 + 3*e4/32 + 45*e6/1024;
    double l3 = 15*e4/256 + 45*e6/1024;
    double l4 = 35*e6/3072;
    double M = a*(l1*lat - l2*sin(2*lat) + l3*sin(4*lat) - l4*sin(6*lat));
//double rho = a*(1-e2) / pow((1-(e*slat1)*(e*slat1)),1.5);
    double nu = a / sqrt(1-(e*slat1)*(e*slat1));
    double p = lon - Grid[to].lon0;
    double k0 = Grid[to].k0;
// y = northing = K1 + K2p2 + K3p4, where
    double K1 = M*k0;
    double K2 = k0*nu*slat1*clat1/2;
    double K3 = (k0*nu*slat1*clat1*clat1sq/24)*(5 - tanlat1sq + 9*eg2*clat1sq + 4*eg2*eg2*clat1sq*clat1sq);
// ING north
    double Y = K1 + K2*p*p + K3*p*p*p*p - Grid[to].false_n;

// x = easting = K4p + K5p3, where
    double K4 = k0*nu*clat1;
    double K5 = (k0*nu*clat1*clat1sq/6)*(1 - tanlat1sq + eg2*clat1*clat1);
// ING east
    double X = K4*p + K5*p*p*p + Grid[to].false_e;

// final rounded results
    *E = (int)(X+0.5);
    *N = (int)(Y+0.5);
}

//======================================================
// Abridged Molodensky transformation between 2 datums
//======================================================
static void Molodensky(double ilat,double ilon,double* olat,double* olon,enum eDatum from,enum eDatum to)
{
// from->WGS84 - to->WGS84 = from->WGS84 + WGS84->to = from->to
    double dX = Datum[from].dX - Datum[to].dX;
    double dY = Datum[from].dY - Datum[to].dY;
    double dZ = Datum[from].dZ - Datum[to].dZ;

    double slat = sin(ilat);
    double clat = cos(ilat);
    double slon = sin(ilon);
    double clon = cos(ilon);
    double ssqlat = slat*slat;

//dlat = ((-dx * slat * clon - dy * slat * slon + dz * clat)
//        + (da * rn * from_esq * slat * clat / from_a)
//        + (df * (rm * adb + rn / adb )* slat * clat))
//       / (rm + from.h);

    double from_f = Datum[from].f;
    double df = Datum[to].f - from_f;
    double from_a = Datum[from].a;
    double da = Datum[to].a - from_a;
    double from_esq = Datum[from].esq;
    double adb = 1.0 / (1.0 - from_f);
    double rn = from_a / sqrt(1 - from_esq * ssqlat);
    double rm = from_a * (1 - from_esq) / pow((1 - from_esq * ssqlat),1.5);
    double from_h = 0.0; // we're flat!

    double dlat = (-dX*slat*clon - dY*slat*slon + dZ*clat
                   + da*rn*from_esq*slat*clat/from_a +
                   + df*(rm*adb + rn/adb)*slat*clat) / (rm+from_h);

// result lat (radians)
    *olat = ilat+dlat;

// dlon = (-dx * slon + dy * clon) / ((rn + from.h) * clat);
    double dlon = (-dX*slon + dY*clon) / ((rn+from_h)*clat);
// result lon (radians)
    *olon = ilon+dlon;
}
