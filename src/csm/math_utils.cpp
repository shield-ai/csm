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
