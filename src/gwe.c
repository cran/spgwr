/*  Copyright by Roger Bivand (C) 2003-5 (with thanks to Chris Brunsdon
 *  for access to his Fortran code)
 */
#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <R_ext/Applic.h>

void gw_gcdist(double *lon1, double *lon2, double *lat1, double *lat2, 
		double *dist);

void gw_dists(double *u, double *v, double *uout, double *vout, 
		int *n, double *dists, int *lonlat); 

void gw_adapt(double *u, double *v, double *uout, double *vout, int *n1, 
		int *n2, double *bw, double *qin, double *d, int *lonlat) 
{
	int N1 = *n1, N2 = *n2, i, index;
	double q = *qin;
	double uo[1], vo[1];
	
	index = floor((N1-1)*q + 0.5); /* + 1 */

	for (i=0; i<N2; i++) {
	    	uo[0] = uout[i];
	    	vo[0] = vout[i];
		gw_dists(u, v, uo, vo, n1, d, lonlat);

		R_rsort(d, N1);
		bw[i] = d[index];
	}
}

void gw_dists(double *u, double *v, double *uout, double *vout, 
		int *n, double *dists, int *lonlat)
{
	int N = *n, j;
	double gc[1], res;
		
	if (lonlat[0] == 0) {
		for (j=0; j<N; j++) {
			res = pythag((u[j]-uout[0]), (v[j]-vout[0]));
			if (R_FINITE(res)) dists[j] = res;
			else dists[j] = (double) 0;
		}
	} else {
		for (j=0; j<N; j++) {
			gw_gcdist(u+j, uout, v+j, vout, gc);
			if (R_FINITE(gc[0])) dists[j] = gc[0];
			else dists[j] = (double) 0;
		}
	}
} 


/* http://home.att.net/~srschmitt/greatcircle.html */

void gw_gcdist(double *lon1, double *lon2, double *lat1, double *lat2, 
		double *dist) {
	
    double F, G, L, sinG2, cosG2, sinF2, cosF2, sinL2, cosL2, S, C;
    double w, R, a, f, D, H1, H2;
    double lat1R, lat2R, lon1R, lon2R, DE2RA;
    
    DE2RA = M_PI/180;
    a = 6378.137;              /* WGS-84 equatorial radius in km */
    f = 1.0/298.257223563;     /* WGS-84 ellipsoid flattening factor */
    
    lat1R = lat1[0]*DE2RA;
    lat2R = lat2[0]*DE2RA;
    lon1R = lon1[0]*DE2RA;
    lon2R = lon2[0]*DE2RA;
    
    F = ( lat1R + lat2R )/2.0;
    G = ( lat1R - lat2R )/2.0;
    L = ( lon1R - lon2R )/2.0;

    sinG2 = R_pow_di( sin( G ), 2 );
    cosG2 = R_pow_di( cos( G ), 2 );
    sinF2 = R_pow_di( sin( F ), 2 );
    cosF2 = R_pow_di( cos( F ), 2 );
    sinL2 = R_pow_di( sin( L ), 2 );
    cosL2 = R_pow_di( cos( L ), 2 );

    S = sinG2*cosL2 + cosF2*sinL2;
    C = cosG2*cosL2 + sinF2*sinL2;

    w = atan( sqrt( S/C ) );
    R = sqrt( S*C )/w;

    D = 2*w*a;
    H1 = ( 3*R - 1 )/( 2*C );
    H2 = ( 3*R + 2 )/( 2*S );

    dist[0] = D*( 1 + f*H1*sinF2*cosG2 - f*H2*cosF2*sinG2 ); 

}

