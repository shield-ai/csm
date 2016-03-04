#ifndef H_LASER_DATA
#define H_LASER_DATA

// KUKA: WINDOWS COMPILATION
#ifndef WINDOWS
#include <sys/time.h>
#else
#include <time.h>
#include <windows.h>
#endif

#include <stdio.h>

#include "restrict.h"

struct correspondence;

typedef struct {
	double p[2];
	double rho, phi;
} point2d;

struct laser_data {
	int nrays;
	double  min_theta;
	double  max_theta;
	
	double * restrict theta;
	
	int    * restrict valid;
	double * restrict readings;
	
	int    * restrict cluster;
	
	double * restrict alpha;
	int    * restrict alpha_valid;

	double * restrict readings_sigma;

	double * restrict true_alpha;
	
	struct correspondence*  restrict corr;

	double true_pose[3];		
	double odometry[3];	
	double estimate[3];	
	

	/** Cartesian representation */
	point2d *  restrict points;
	/** Cartesian representation, in "world" (laser_ref) coordinates. 
	    Computed using ld_compute_world_coords() */
	point2d *  restrict points_w;

	/** Timestamp */
	struct timeval tv;
	char hostname[32];


	/* Jump tables needed by find_correspondences_tricks(). */
	int * restrict up_bigger, 
	    * restrict up_smaller, 
	    * restrict down_bigger, 
	    * restrict down_smaller;	
};

struct correspondence {
	/** 1 if this correspondence is valid  */
	int valid; 
	/** Closest point in the other scan.  */
	int j1;
	/** Second closest point in the other scan.  */
	int j2;
	/** Type of correspondence (point to point, or point to line) */
	enum { corr_pp = 0, corr_pl = 1} type;
	/** Squared distance from p(i) to point j1 */
	double dist2_j1; 
};

typedef struct laser_data* LDP;

/** This returns a new structure, with all fields initialized */
LDP ld_alloc_new(int nrays);

/** This DOES free() the pointer  */
void ld_free(LDP);

/** This allocs the fields in the given structure. Use ld_alloc_new(), not this. */
void ld_alloc(LDP, int nrays);

/** This does NOT free the pointer. Don't use -- use ld_alloc_new()/ld_free() instead. */
void ld_dealloc(LDP);

/** Fills the x,y fields in "points" by transforming (theta, reading) to cartesian */
void ld_compute_cartesian(LDP);

/** Computes the "points_w" coordinates by roto-translating "points" */
void ld_compute_world_coords(LDP, const double *pose);

/** Fills the fields: *up_bigger, *up_smaller, *down_bigger, *down_smaller.*/
void ld_create_jump_tables(LDP);

/** Computes an hash of the correspondences */
unsigned int ld_corr_hash(LDP);

/** Returns the number of valid correspondences. */
int ld_num_valid_correspondences(LDP);

/** Do an extensive sanity check about the data contained in the structure. */
int ld_valid_fields(LDP);

/** A simple clustering algorithm. Sets the `cluster' field in the structure. */
void ld_simple_clustering(LDP ld, double threshold);

/** A cool orientation estimation algorithm. Needs cluster. */
void ld_compute_orientation(LDP ld, int size_neighbourhood, double sigma);

void possible_interval(
	const double *p_i_w, LDP laser_sens, 
	double max_angular_correction_deg, double max_linear_correction, int*from, int*to, int*start_cell);


#include "laser_data_inline.h"

#endif

