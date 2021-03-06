#include <math.h>
#include <iostream>

#include "icp.h"
#include <csm/flags.h>
#include <vector>

double hoare_selection(double *data, int start, int end, int k);

/** expects cartesian valid */
void visibilityTest(LDP laser_ref, const gsl_vector*u) {

	//double theta_from_u[laser_ref->nrays];
	std::vector<double> theta_from_u(laser_ref->nrays, 0.0);
	
	int j;
	for(j=0;j<laser_ref->nrays;j++) {
		if(!ld_valid_ray(laser_ref,j)) continue;
		theta_from_u[j] = 
			atan2(gvg(u,1)-laser_ref->points[j].p[1],
			      gvg(u,0)-laser_ref->points[j].p[0]);
	}
	
    sm_debug("\tvisibility: Found outliers: ");
	for(j=1;j<laser_ref->nrays;j++) {
		if(!ld_valid_ray(laser_ref,j)||!ld_valid_ray(laser_ref,j-1)) continue;
		if(theta_from_u[j]<theta_from_u[j-1]) {
			laser_ref->valid[j] = 0;
            sm_debug("%d ",j);
		}
	}
    sm_debug("\n");
}


/** 
	If multiple points in laser_sens match to the same point in laser_ref, 
	only the nearest one wins.

	Uses: laser_sens->corr, laser_sens->p
	Modifies: laser_sens->corr
 */
void kill_outliers_double(struct sm_params*params) {
	const double threshold = 3; /* TODO: add as configurable */

	LDP laser_ref  = params->laser_ref;
	LDP laser_sens = params->laser_sens;

	//double dist2_i[laser_sens->nrays];
	std::vector<double> dist2_i(laser_sens->nrays, 0.0);
	//double dist2_j[laser_ref->nrays];
	std::vector<double> dist2_j(laser_ref->nrays, 1000000.);
	
	int i;
	for(i=0;i<laser_sens->nrays;i++) {
		if(!ld_valid_corr(laser_sens, i)) continue;
		int j1 = laser_sens->corr[i].j1;
		dist2_i[i] = laser_sens->corr[i].dist2_j1;
		dist2_j[j1] = (std::min)(dist2_j[j1], dist2_i[i]);
	}
	
#ifndef ENABLE_OPTIMIZATION
  int nkilled=0;
#endif
	for(i=0;i<laser_sens->nrays;i++) {
		if(!ld_valid_corr(laser_sens, i)) continue;
		int j1 = laser_sens->corr[i].j1;
		if(dist2_i[i] > (threshold*threshold)*dist2_j[j1]) {
			laser_sens->corr[i].valid=0;
#ifndef ENABLE_OPTIMIZATION
      nkilled++;
#endif
		}
	}
    sm_debug("\tkill_outliers_double: killed %d correspondences\n",nkilled);
}
	
/** 
	Trims the corrispondences using an adaptive algorithm 

	Assumes cartesian coordinates computed. (points and points_w)
	 
	So, to disable this:
		outliers_maxPerc = 1
		outliers_adaptive_order = 1 

*/
void kill_outliers_trim(struct sm_params*params,  double*total_error) {

	LDP laser_ref  = params->laser_ref;
	LDP laser_sens = params->laser_sens;
	
	/* dist2, indexed by k, contains the error for the k-th correspondence */
	int k = 0; 
	//double dist2[laser_sens->nrays];
	std::vector<double> dist2(laser_sens->nrays, 0.0);
		
	int i;
	//double dist[laser_sens->nrays];
	std::vector<double> dist(laser_sens->nrays, 0.0);
	/* for each point in laser_sens */
	for(i=0;i<laser_sens->nrays;i++) {
		/* which has a valid correspondence */
		if(!ld_valid_corr(laser_sens, i)) { dist[i]=std::numeric_limits<double>::quiet_NaN(); continue; }
		double *p_i_w = laser_sens->points_w[i].p;
		
		int j1 = laser_sens->corr[i].j1;
		int j2 = laser_sens->corr[i].j2;
		/* Compute the distance to the corresponding segment */
		dist[i]=  dist_to_segment_d(
			laser_ref->points[j1].p, laser_ref->points[j2].p, p_i_w);
		dist2[k] = dist[i];
		k++;	
	}
	
	

#if 0	
	double dist2_copy[k]; for(i=0;i<k;i++) dist2_copy[i] = dist2[i];
#endif 


	double maxPerc = params->outliers_maxPerc;
	double adaptive_order = params->outliers_adaptive_order;
	double adaptive_mult = params->outliers_adaptive_mult;
	if(*total_error < 0.0) // PHL - hack to tell kill_outliers_trim to go easy on outliers (gets re-initialized to 0 below)
	{
		maxPerc = params->outliers_maxPerc_reduced;
		adaptive_order = params->outliers_adaptive_order_reduced;
		adaptive_mult = params->outliers_adaptive_mult_reduced;
	}


	/* two errors limits are defined: */
		/* In any case, we don't want more than outliers_maxPerc% */
	int order = (int)floor(k*maxPerc);
			order = (std::max)(0, (std::min)(order, k-1));

	/* The dists for the correspondence are sorted
	   in ascending order */
    std::sort(dist2.begin(), dist2.begin()+k);
	
		/* Then we take a order statics (o*K) */
		/* And we say that the error must be less than alpha*dist(o*K) */
		int order2 = (int)floor(k*adaptive_order);
			order2 = (std::max)(0, (std::min)(order2, k-1));
		double error_limit2 = adaptive_mult*dist2[order2];
	
	double error_limit = (std::min)(dist2[order], error_limit2);

#if 0
  double error_limit1_ho = hoare_selection(dist2_copy, 0, k-1, order);
  double error_limit2_ho = error_limit2;
  if((error_limit1_ho != error_limit1) || (error_limit2_ho != error_limit2)) {
    printf("%f == %f    %f  == %f\n",
        error_limit1_ho, error_limit1, error_limit2_ho, error_limit2);
  }
#endif

  sm_debug("\ticp_outliers: maxPerc %f error_limit: fix %f adaptive %f \n",
      params->outliers_maxPerc,dist2[order],error_limit2);


	*total_error = 0;
#ifndef ENABLE_OPTIMIZATION
  int nvalid = 0;
#endif
	for(i=0;i<laser_sens->nrays;i++) {
		if(!ld_valid_corr(laser_sens, i)) continue;
		if(dist[i] > error_limit) {
			laser_sens->corr[i].valid = 0;
			laser_sens->corr[i].j1 = -1;
			laser_sens->corr[i].j2 = -1;
		} else {
			*total_error += dist[i];
#ifndef ENABLE_OPTIMIZATION
      nvalid++;
#endif
		}
	}
    sm_debug("\ticp_outliers: valid %d/%d (limit: %f) mean error = %f \n",nvalid,k,error_limit, *total_error/nvalid);
}
#if 0
double hoare_selection(double *data, int start, int end, int k)
{
  int pivotIndex, i, j;
  float pivotValue, tmp;

  while(start < end) {
	  //select a random pivot
	pivotIndex = start + (end-start)/2;
	 pivotValue = data[pivotIndex];
				 //sort the array into two sub arrays around the pivot value
	 i = start;
	 j = end;
	 while(i <= j) {
		while(data[i] < pivotValue)
		  i++;
		while(data[j] > pivotValue)
		  j--;
		if(i<=j) {
		  tmp = data[i];
		  data[i] = data[j];
		  data[j] = tmp;
		  i++;
		  j--;
		}
	 }
			 //continue search in left sub array
	 if(j < k)
		start = i;
	 if(k < i) //continue search in right sub array
		end = j;
  }
  return(data[k]);
}

#endif
