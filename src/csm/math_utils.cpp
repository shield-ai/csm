#include "csm_all.h"

void possible_interval(
	const double *p_i_w, LDP ld, 
	double max_angular_correction_deg, double max_linear_correction, int*from, int*to, int*start_cell) 
{
	double angle_res = (ld->max_theta-ld->min_theta)/ld->nrays;

	/* Delta for the angle */
	double delta = fabs(deg2rad(max_angular_correction_deg)) +
	        fabs(atan(max_linear_correction/norm_d(p_i_w)));

	/* Dimension of the cell range */
	int range = (int) ceil(delta/angle_res);

	/* To be turned into an interval of cells */
	double start_theta = atan2(p_i_w[1], p_i_w[0]);
	
	/* Make sure that start_theta is in the interval [min_theta,max_theta]. 
	   For example, -1 is not in [0, 2pi] */
	if(start_theta<ld->min_theta) start_theta += 2*M_PI;
	if(start_theta>ld->max_theta) start_theta -= 2*M_PI;
	
	*start_cell  = (int)
		((start_theta - ld->min_theta) / (ld->max_theta-ld->min_theta) * ld->nrays);

	*from = minmax(0,ld->nrays-1, *start_cell-range);
	*to =   minmax(0,ld->nrays-1, *start_cell+range);

	if(0)
	printf("from: %d to: %d delta: %f start_theta: %f min/max theta: [%f,%f] range: %d start_cell: %d\n",
		*from, *to,
		delta,start_theta,ld->min_theta,ld->max_theta, range, *start_cell);
}

#ifndef ENABLE_OPTIMIZATION
double normalize_0_2PI(double t) {
	if(std::isnan(t)) {
		sm_error("Passed NAN to normalize_0_2PI().\n");
		return std::numeric_limits<double>::quiet_NaN();
	}
	while(t<0) t+=2*M_PI;
	while(t>=2*M_PI) t-=2*M_PI;
	return t;
}

/* Executes ray tracing for a segment. p0 and p1 are the segments extrema, eye is the position
of the eye, and direction is the direction of the ray coming out of the eye. Returns true
if the ray intersects the segment, and in that case *range contains the length of the ray. */
int segment_ray_tracing(const double p0[2], const double p1[2], const double eye[2], double direction, double*range) {
	
	*range = std::numeric_limits<double>::quiet_NaN();
	
	// p0 - p1
	double arrow[2] = {p1[0]-p0[0],p1[1]-p0[1]};
	// Normal to segment line
	double S[2] = { -arrow[1], arrow[0]};
	// Viewing direction
	double N[2] = { cos(direction), sin(direction)};
	// If S*N = 0 then they cannot cross
	double S_dot_N = dot_d(S,N);
	if( S_dot_N == 0) return 0;
	// Rho of the line in polar coordinates (multiplied by |S|)
	double line_rho = dot_d(p0,S);
	// Rho of the eye  (multiplied by |S|)
	double eye_rho = dot_d(eye,S);
	// Black magic
	double dist = (line_rho - eye_rho) / S_dot_N;
	if(dist<=0) return 0;
	
	// Now we check whether the crossing point
	// with the line lies within the segment
	
	// Crossing point
	double crossing[2] = {eye[0] + N[0]*dist, eye[1]+N[1]*dist};
	// Half of the segment
	double midpoint[2] = { 0.5*(p1[0]+p0[0]),0.5*(p1[1]+p0[1])};
	
	double seg_size = distance_d(p0, p1);
	double dist_to_midpoint = distance_d(crossing, midpoint);
	
	if(dist_to_midpoint > seg_size/2 )
		return 0;
	
	*range = dist;
	return 1;
}

double segment_alpha(const double p0[2], const double p1[2]) {
	double arrow[2] = {p1[0]-p0[0],p1[1]-p0[1]};
	// Normal to segment line
	double S[2] = { -arrow[1], arrow[0]};
	return atan2(S[1], S[0]);
}


double max_in_array(const double*v, int n) {
	assert(n>0);
	double m = v[0];
	int i; 
	for(i=0;i<n;i++)
		if(v[i]>m) m = v[i];
	return m;
}
#endif
