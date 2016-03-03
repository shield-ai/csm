#include "egsl.h"

val egsl_vFda(size_t rows, size_t cols, const double *a) {
	val v = egsl_alloc(rows, cols);

	size_t i; size_t j;
	for(i=0;i<rows;i++)
	for(j=0;j<cols;j++) {
		*egsl_atmp(v,i,j) = a[j+i*cols];
	}
	return v;
}

val egsl_vFa(size_t rows, const double*a) {
	val v = egsl_alloc(rows,1);
	size_t i;
	for(i=0;i<rows;i++)
		*egsl_atmp(v,i,0) =  a[i];
	return v;
}

gsl_matrix* egsl_v2gslm(val v){
	gsl_matrix * m = egsl_gslm(v); 
	gsl_matrix * m2 = gsl_matrix_alloc(m->rows(),m->cols());
	gsl_matrix_memcpy(m2,m);
	return m2;
}
