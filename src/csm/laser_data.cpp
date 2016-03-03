#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "csm_all.h"

/* -------------------------------------------------- */

LDP ld_alloc_new(int nrays) {
	LDP ld = (LDP) malloc(sizeof(struct laser_data));
	ld_alloc(ld, nrays);
	return ld;
}

template <typename T>
T* alloc_array(int n, T def) {
  T *v = (T*) malloc(sizeof(T)*n);
  for (int i=0; i<n; ++i) {
    v[i] = def;
  }
  return v;
}

template <typename T>
T* alloc_zero_array(int n) {
  T *v = (T*) malloc(sizeof(T)*n);
  memset(v, 0, (sizeof(T)*n));
  return v;
}

void ld_alloc(LDP ld, int nrays) {
	ld->nrays = nrays;
	
	ld->valid        = alloc_zero_array<int>(nrays);
	ld->readings     = alloc_array(nrays, std::numeric_limits<float>::quiet_NaN());
	ld->readings_sigma = alloc_array(nrays, std::numeric_limits<double>::quiet_NaN());
	ld->theta        = alloc_array(nrays, std::numeric_limits<double>::quiet_NaN());
	
	ld->min_theta = std::numeric_limits<double>::quiet_NaN();
	ld->max_theta = std::numeric_limits<double>::quiet_NaN();
	
	ld->cluster      = alloc_array(nrays, -1);
	ld->alpha        = alloc_array(nrays, std::numeric_limits<double>::quiet_NaN());
	ld->alpha_valid  = alloc_zero_array<int>(nrays);

	ld->true_alpha   = alloc_array(nrays, std::numeric_limits<double>::quiet_NaN());
	
	ld->up_bigger    = alloc_zero_array<int>(nrays);
	ld->up_smaller   = alloc_zero_array<int>(nrays);
	ld->down_bigger  = alloc_zero_array<int>(nrays);
	ld->down_smaller = alloc_zero_array<int>(nrays);

	ld->corr = (struct correspondence*) 
		malloc(sizeof(struct correspondence)*nrays);

	int i;
	for(i=0;i<ld->nrays;i++) {
		ld->corr[i].valid = 0;
		ld->corr[i].j1 = -1;
		ld->corr[i].j2 = -1;
	}
	
	for(i=0;i<3;i++) {
		ld->odometry[i] = 
		ld->estimate[i] = 
		ld->true_pose[i] = std::numeric_limits<double>::quiet_NaN();
	}
	
	ld->points = (point2d*) malloc(nrays * sizeof(point2d));
	ld->points_w = (point2d*) malloc(nrays * sizeof(point2d));
	
	for(i=0;i<nrays;i++) {
		ld->points[i].p[0] = 
		ld->points[i].p[1] = 
		ld->points[i].rho = 
		ld->points[i].phi = std::numeric_limits<double>::quiet_NaN();
		ld->points_w[i] = ld->points[i];
	}
	
	strcpy(ld->hostname, "CSM");
}

void ld_free(LDP ld) {
	ld_dealloc(ld);
	free(ld);
}

void ld_dealloc(struct laser_data*ld){	
	free(ld->valid);
	free(ld->readings);
	free(ld->readings_sigma);
	free(ld->theta);
	free(ld->cluster);
	free(ld->alpha);
	free(ld->alpha_valid);
	free(ld->true_alpha);
	free(ld->up_bigger);
	free(ld->up_smaller);
	free(ld->down_bigger);
	free(ld->down_smaller);
	free(ld->corr);
	
/*	int i;
	for(i=0;i<ld->nrays;i++)
		gsl_vector_free(ld->p[i]);
	free(ld->p);*/

	free(ld->points);
	free(ld->points_w);
}


void ld_compute_cartesian(LDP ld) {
	int i;
	for(i=0;i<ld->nrays;i++) {
/*		if(!ld_valid_ray(ld,i)) continue;*/
		double x = std::cos(ld->theta[i])*ld->readings[i];
		double y = std::sin(ld->theta[i])*ld->readings[i];
		
		ld->points[i].p[0] = x, 
		ld->points[i].p[1] = y;
		ld->points[i].rho = std::numeric_limits<double>::quiet_NaN();
		ld->points[i].phi = std::numeric_limits<double>::quiet_NaN();
	}
}


void ld_compute_world_coords(LDP ld, const double *pose) {
	double cos_theta = std::cos(pose[2]); 
	double sin_theta = std::sin(pose[2]);

	point2d * points = ld->points;
	point2d * points_w = ld->points_w;
	for(int i=0;i<ld->nrays;i++) {
		if(!ld_valid_ray(ld,i)) continue;
		double x = points[i].p[0], 
		       y = points[i].p[1]; 
		
		if(std::isnan(x) || std::isnan(y)) {
			sm_error("ld_compute_world_coords(): I expected that cartesian coords were already computed: ray #%d: %f %f.\n", i, x, y);
		}
		
		points_w[i].p[0] = cos_theta * x -sin_theta*y + pose[0];
		points_w[i].p[1] = sin_theta * x +cos_theta*y + pose[1];
		/* polar coordinates */
		x = points_w[i].p[0];
		y = points_w[i].p[1];
		points_w[i].rho = sqrt( x*x+y*y);
		points_w[i].phi = atan2(y, x);
	}
}



int ld_num_valid_correspondences(LDP ld) {
	int i; 
	int num = 0;
	for(i=0;i<ld->nrays;i++) {
		if(ld->corr[i].valid)
			num++;
	}
	return num;
}


int ld_valid_fields(LDP ld)  {
	if(!ld) {
		sm_error("NULL pointer given as laser_data*.\n");	
		return 0;
	}
	
	const int min_nrays = 10;
	const int max_nrays = 10000;
	if(ld->nrays < min_nrays || ld->nrays > max_nrays) {
		sm_error("Invalid number of rays: %d\n", ld->nrays);
		return 0;
	}
	if(std::isnan(ld->min_theta) || std::isnan(ld->max_theta)) {
		sm_error("Invalid min / max theta: min_theta = %f max_theta = %f\n",
			ld->min_theta, ld->max_theta);
		return 0;
	}
	const double min_fov = deg2rad(20.0); 
	const double max_fov = 2.01 * M_PI;
	double fov = ld->max_theta - ld->min_theta;
	if( fov < min_fov || fov > max_fov) {
		sm_error("Strange FOV: %f rad = %f deg \n", fov, rad2deg(fov));
		return 0;
	}
	if(fabs(ld->min_theta - ld->theta[0]) > 1e-8) {
		sm_error("Min_theta (%f) should be theta[0] (%f)\n",
			ld->min_theta, ld->theta[0]);
		return 0;
	}
	if(fabs(ld->max_theta - ld->theta[ld->nrays-1]) > 1e-8) {
		sm_error("Min_theta (%f) should be theta[0] (%f)\n",
			ld->max_theta, ld->theta[ld->nrays-1]);
		return 0;
	}
	/* Check that there are valid rays */
	const double min_reading = 0;
	const double max_reading = 100;
	int i; for(i=0;i<ld->nrays;i++) {
		if(ld->valid[i]) {
			double r = ld->readings[i];
			if(std::isnan(r) || std::isnan(ld->theta[i])) {
				sm_error("Ray #%d: r = %f  theta = %f but valid is %d\n",
					i, r, ld->theta[i], ld->valid[i]);
				return 0;
			}
			if( !( min_reading < r && r < max_reading ) ) {
				sm_error("Ray #%d: %f is not in interval (%f, %f) \n",
					i, r, min_reading, max_reading);
				return 0;
			}		
		} else {
			/* ray not valid, but checking theta anyway */
			if(std::isnan(ld->theta[i])) {
				sm_error("Ray #%d: valid = %d  but theta = %f\n",
					i,  ld->valid[i], ld->theta[i]);
				return 0;
			}

			if(ld->cluster[i] != -1 ) {
				sm_error("Invalid ray #%d has cluster %d\n.", i, ld->cluster[i]);
				return 0;
			}
		}
		if(ld->cluster[i] < -1 ) {
			sm_error("Ray #%d: Invalid cluster value %d\n.", i, ld->cluster[i]);
			return 0;
		}
		
		if(!std::isnan(ld->readings_sigma[i]) && ld->readings_sigma[i] < 0) {
			sm_error("Ray #%d: has invalid readings_sigma %f \n", i, ld->readings_sigma[i]);
			return 0;
		}
		
	}
	/* Checks that there is at least 10% valid rays */
	int num_valid   = count_equal(ld->valid, ld->nrays, 1);
	if (num_valid < ld->nrays * 0.10) {
    int num_invalid = count_equal(ld->valid, ld->nrays, 0);
		sm_error("Valid: %d/%d invalid: %d.\n", num_valid, ld->nrays, num_invalid);
		return 0;
	}

	return 1;
}


/** Computes an hash of the correspondences */
unsigned int ld_corr_hash(LDP ld){
	unsigned int hash = 0;
	unsigned int i    = 0;

	for(i = 0; i < (unsigned)ld->nrays; i++) {
		int str = ld_valid_corr(ld, (int)i) ? (ld->corr[i].j1 + 1000*ld->corr[i].j2) : -1;
		hash ^= ((i & 1) == 0) ? (  (hash <<  7) ^ (str) ^ (hash >> 3)) :
		                         (~((hash << 11) ^ (str) ^ (hash >> 5)));
	}

	return (hash & 0x7FFFFFFF);
}
