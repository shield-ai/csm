#ifndef H_MATH_UTILS
#define H_MATH_UTILS


/* Sometimes I really don't understand compilers.. */ 
#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

/** Inlinable local functions **/
inline int minmax(int from, int to, int x) {
	return (std::max)((std::min)(x,to),from);
}

/** Templated math functions */
/** Square value **/
template <typename T>
inline T square(T x) {
  return x*x;
}

/** Calculate distance squared **/
template <typename T>
inline T distance_squared_d(const T a[2], const T b[2]) {
	return square(a[0]-b[0]) + square(a[1]-b[1]);
}

/** Calculate normal **/
template <typename T>
inline T norm_d(const T p[2]) {
	return std::sqrt(p[0]*p[0]+p[1]*p[1]);
}

/** Degrees to radians */
template <typename T>
constexpr T deg2rad(T deg) {
	return deg * (M_PI / 180);
}

/** Radians to degrees */
template <typename T>
constexpr T rad2deg(T rad) {
	return rad * (180 / M_PI);	
}

/** Copies n doubles from from to to */
template <typename T>
void copy_d(const T *from, int n, T *to) {
	int i; for(i=0;i<n;i++) to[i] = from[i];
}

/** Returns true if any value in d is NAN */
template <typename T>
int any_nan(const T *d, int n) {
	int i; for(i=0;i<n;i++) 
		if(std::isnan(d[i]))
			return 1;
	return 0;
}

/** Count numbers of items in array v equal to value */
template <typename T>
int count_equal(const T *v, int n, T value) {
	int num = 0, i;
	for(i=0;i<n;i++) if(value == v[i]) num++;
	return num;
}

/** Returns an angle difference in the [-pi, pi] range */
template <typename T>
T angleDiff(T a, T b) {
	T t = a - b;
	while(t<-M_PI) t+= 2*M_PI;
	while(t>M_PI)  t-= 2*M_PI;
	return t;
}

/** These are the operators defined in Smith & Cheeseman  */
/** Note: ominus_d and oplus_d coperators combined for performance */
template <typename T>
void pose_diff_d(const T pose2[3], const T pose1[3], T res[3]) {
  T temp[2];
	T c = std::cos(pose1[2]);
	T s = std::sin(pose1[2]);
	temp[0] = -c*pose1[0]-s*pose1[1];
	temp[1] =  s*pose1[0]-c*pose1[1];
	s = -s;
	res[0]=temp[0]+c*pose2[0]-s*pose2[1];
	res[1]=temp[1]+s*pose2[0]+c*pose2[1];
	res[2]=-pose1[2]+pose2[2];
	
	while(res[2] > +M_PI) res[2] -= 2*M_PI;
	while(res[2] < -M_PI) res[2] += 2*M_PI;
}

/** Projects (p[0],p[1]) on the LINE passing through (ax,ay)-(bx,by). If distance!=0, distance is set
to the distance from the point to the segment */
template <typename T>
void projection_on_line_d(const T a[2], const T b[2], const T p[2], T res[2], T *distance) {
	T t0 = a[0]-b[0];
	T t1 = a[1]-b[1];
	T one_on_r = 1 / sqrt(t0*t0+t1*t1);
	/* normal */
	T nx = t1  * one_on_r ;
	T ny = -t0 * one_on_r ;
	T c= nx, s = ny; 
	T rho = c*a[0]+s*a[1];

	res[0] =   c*rho + s*s*p[0] - c*s*p[1] ;
	res[1] =   s*rho - c*s*p[0] + c*c*p[1] ;	
	
	if(distance)
		*distance = fabs(rho-(c*p[0]+s*p[1]));
}

template <typename T>
void projection_on_segment_d(const T a[2], const T b[2], const T x[2], T proj[2]) 
{
	projection_on_line_d(a,b,x,proj,static_cast<double*>(nullptr));
	if ((proj[0]-a[0])*(proj[0]-b[0]) +
	    (proj[1]-a[1])*(proj[1]-b[1]) < 0 ) {
		/* the projection is inside the segment */
	} else 
		if(distance_squared_d(a,x) < distance_squared_d(b,x)) 
			copy_d(a,2,proj);
		else
			copy_d(b,2,proj);
}

/** Distance of x from its projection on segment a-b */
template <typename T>
T dist_to_segment_d(const T a[2], const T b[2], const T x[2]) {
	T proj[2]; T distance;
	projection_on_line_d(a,b,x,proj, &distance);
	if ((proj[0]-a[0])*(proj[0]-b[0]) +
	    (proj[1]-a[1])*(proj[1]-b[1]) < 0 ) {
		/* the projection is inside the segment */
		return distance;
	} else 
		return sqrt((std::min)(distance_squared_d(a,x), distance_squared_d(b,x)));
}

/** Templated utility function **/
/** A function to print poses and covariances in a friendly way */
static char tmp_buf[1024];
template <typename T>
const char* friendly_pose(const T *pose) {
	sprintf(tmp_buf, "(%4.2f mm, %4.2f mm, %4.4f deg)",
		1000*pose[0],1000*pose[1],rad2deg(pose[2]));
	return tmp_buf;
}

#endif
