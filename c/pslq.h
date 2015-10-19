#include <math.h>

#define ABS(x) (((x)<0)?-(x):(x))
#define MIN(a, b) ( ((a)<(b)) ? (a) : (b) )
#define MAX(a, b) ( ((a)>(b)) ? (a) : (b) )
#define SQR(a) ( (a)*(a) )

/*
   This is the PSLQ algorithm as usually implemented except that
   the matrix B is replaced by its transpose in order to optimize
   permutations of columns.
*/
#define pslq(x,n,r) { \
 \
  double PSLQ_hygienic_s[n]; \
  double PSLQ_hygienic_y[n]; \
  double PSLQ_hygienic_a_[(n)*(n)] = { 0.0 } ; \
  double PSLQ_hygienic_b_[(n)*(n)] = { 0.0 }; \
  double PSLQ_hygienic_h_[(n)*((n)-1)] = { 0.0 }; \
  double *PSLQ_hygienic_a[n]; \
  double *PSLQ_hygienic_b[n]; \
  double *PSLQ_hygienic_h[n]; \
  double PSLQ_hygienic_t; \
  double *PSLQ_hygienic_swap; \
  int PSLQ_hygienic_i,PSLQ_hygienic_j,PSLQ_hygienic_k; \
  const double PSLQ_hygienic_teps = 16.0 * 1.0e-10; \
  const double PSLQ_hygienic_gam = 1.2; \
 \
  for(PSLQ_hygienic_i=0;PSLQ_hygienic_i<(n);PSLQ_hygienic_i++) { \
      PSLQ_hygienic_a[PSLQ_hygienic_i] = PSLQ_hygienic_a_ + PSLQ_hygienic_i*(n); \
      PSLQ_hygienic_b[PSLQ_hygienic_i] = PSLQ_hygienic_b_ + PSLQ_hygienic_i*(n); \
      PSLQ_hygienic_h[PSLQ_hygienic_i] = PSLQ_hygienic_h_ + PSLQ_hygienic_i*((n)-1); \
      PSLQ_hygienic_a[PSLQ_hygienic_i][PSLQ_hygienic_i] = 1.0; \
      PSLQ_hygienic_b[PSLQ_hygienic_i][PSLQ_hygienic_i] = 1.0; \
  } \
 \
  PSLQ_hygienic_t = SQR((x)[(n)-1]); PSLQ_hygienic_s[(n)-1] = ABS( (x)[(n)-1] ); \
  for(PSLQ_hygienic_i=(n)-2;PSLQ_hygienic_i>=0;PSLQ_hygienic_i--) { \
    PSLQ_hygienic_t += SQR((x)[PSLQ_hygienic_i]); PSLQ_hygienic_s[PSLQ_hygienic_i] = __builtin_sqrt( PSLQ_hygienic_t ); \
  } \
 \
  PSLQ_hygienic_t = PSLQ_hygienic_s[0]; \
  for(PSLQ_hygienic_i=0;PSLQ_hygienic_i<(n);PSLQ_hygienic_i++) { \
    PSLQ_hygienic_y[PSLQ_hygienic_i] = (x)[PSLQ_hygienic_i] / PSLQ_hygienic_t; \
  } \
  PSLQ_hygienic_s[0] = 1.0; \
  for(PSLQ_hygienic_i=1;PSLQ_hygienic_i<(n);PSLQ_hygienic_i++) { \
     PSLQ_hygienic_s[PSLQ_hygienic_i] /= PSLQ_hygienic_t; \
  } \
 \
  for(PSLQ_hygienic_i=0; PSLQ_hygienic_i<(n); PSLQ_hygienic_i++) { \
    for(PSLQ_hygienic_j=0; PSLQ_hygienic_j<=MIN(PSLQ_hygienic_i,(n)-2); PSLQ_hygienic_j++) { \
      PSLQ_hygienic_h[PSLQ_hygienic_i][PSLQ_hygienic_j] = (PSLQ_hygienic_i==PSLQ_hygienic_j)? PSLQ_hygienic_s[PSLQ_hygienic_j+1]/PSLQ_hygienic_s[PSLQ_hygienic_j] : \
            -PSLQ_hygienic_y[PSLQ_hygienic_i]*PSLQ_hygienic_y[PSLQ_hygienic_j] \
            /(PSLQ_hygienic_s[PSLQ_hygienic_j]*PSLQ_hygienic_s[PSLQ_hygienic_j+1]); \
    } \
  } \
 \
  for(PSLQ_hygienic_i=1;PSLQ_hygienic_i<(n);PSLQ_hygienic_i++) { \
    for(PSLQ_hygienic_j=PSLQ_hygienic_i-1; PSLQ_hygienic_j>=0; PSLQ_hygienic_j--) { \
      PSLQ_hygienic_t = __builtin_nearbyint(PSLQ_hygienic_h[PSLQ_hygienic_i][PSLQ_hygienic_j] / PSLQ_hygienic_h[PSLQ_hygienic_j][PSLQ_hygienic_j] ); \
      PSLQ_hygienic_y[PSLQ_hygienic_j] += PSLQ_hygienic_t * PSLQ_hygienic_y[PSLQ_hygienic_i]; \
      for(PSLQ_hygienic_k=0;PSLQ_hygienic_k<=PSLQ_hygienic_j;PSLQ_hygienic_k++) { PSLQ_hygienic_h[PSLQ_hygienic_i][PSLQ_hygienic_k] -= PSLQ_hygienic_t * PSLQ_hygienic_h[PSLQ_hygienic_j][PSLQ_hygienic_k]; } \
      for(PSLQ_hygienic_k=0;PSLQ_hygienic_k<(n);PSLQ_hygienic_k++) { \
        PSLQ_hygienic_a[PSLQ_hygienic_i][PSLQ_hygienic_k] -= PSLQ_hygienic_t * PSLQ_hygienic_a[PSLQ_hygienic_j][PSLQ_hygienic_k]; \
        PSLQ_hygienic_b[PSLQ_hygienic_j][PSLQ_hygienic_k] += PSLQ_hygienic_t * PSLQ_hygienic_b[PSLQ_hygienic_i][PSLQ_hygienic_k]; \
      } \
    } \
  } \
 \
  int PSLQ_hygienic_m; \
  int PSLQ_hygienic_process = 1; \
 \
  while(PSLQ_hygienic_process) { \
    double PSLQ_hygienic_m_val = -1.0; \
    double g = PSLQ_hygienic_gam; \
    PSLQ_hygienic_m = -1; \
 \
    for(PSLQ_hygienic_i=0;PSLQ_hygienic_i<(n)-1;PSLQ_hygienic_i++, g*=PSLQ_hygienic_gam) { \
      PSLQ_hygienic_t = g*ABS(PSLQ_hygienic_h[PSLQ_hygienic_i][PSLQ_hygienic_i]); \
      if (PSLQ_hygienic_t>PSLQ_hygienic_m_val) { PSLQ_hygienic_m_val = PSLQ_hygienic_t; PSLQ_hygienic_m = PSLQ_hygienic_i; } \
    } \
 \
    PSLQ_hygienic_t = PSLQ_hygienic_y[PSLQ_hygienic_m]; PSLQ_hygienic_y[PSLQ_hygienic_m] = PSLQ_hygienic_y[PSLQ_hygienic_m+1]; \
          PSLQ_hygienic_y[PSLQ_hygienic_m+1] = PSLQ_hygienic_t; \
    PSLQ_hygienic_swap = PSLQ_hygienic_a[PSLQ_hygienic_m]; PSLQ_hygienic_a[PSLQ_hygienic_m] = PSLQ_hygienic_a[PSLQ_hygienic_m+1]; PSLQ_hygienic_a[PSLQ_hygienic_m+1] = PSLQ_hygienic_swap; \
    PSLQ_hygienic_swap = PSLQ_hygienic_b[PSLQ_hygienic_m]; PSLQ_hygienic_b[PSLQ_hygienic_m] = PSLQ_hygienic_b[PSLQ_hygienic_m+1]; PSLQ_hygienic_b[PSLQ_hygienic_m+1] = PSLQ_hygienic_swap; \
    PSLQ_hygienic_swap = PSLQ_hygienic_h[PSLQ_hygienic_m]; PSLQ_hygienic_h[PSLQ_hygienic_m] = PSLQ_hygienic_h[PSLQ_hygienic_m+1]; PSLQ_hygienic_h[PSLQ_hygienic_m+1] = PSLQ_hygienic_swap; \
 \
    if (PSLQ_hygienic_m<(n)-2) { \
      double PSLQ_hygienic_t0, PSLQ_hygienic_t1, PSLQ_hygienic_t2, PSLQ_hygienic_t3, PSLQ_hygienic_t4; \
      PSLQ_hygienic_t0 = sqrt( SQR(PSLQ_hygienic_h[PSLQ_hygienic_m][PSLQ_hygienic_m]) + SQR(PSLQ_hygienic_h[PSLQ_hygienic_m][PSLQ_hygienic_m+1]) ); \
      PSLQ_hygienic_t1 = PSLQ_hygienic_h[PSLQ_hygienic_m][PSLQ_hygienic_m] / PSLQ_hygienic_t0; \
      PSLQ_hygienic_t2 = PSLQ_hygienic_h[PSLQ_hygienic_m][PSLQ_hygienic_m+1] / PSLQ_hygienic_t0; \
      for(PSLQ_hygienic_i=PSLQ_hygienic_m;PSLQ_hygienic_i<(n);PSLQ_hygienic_i++) { \
        PSLQ_hygienic_t3 = PSLQ_hygienic_h[PSLQ_hygienic_i][PSLQ_hygienic_m]; \
        PSLQ_hygienic_t4 = PSLQ_hygienic_h[PSLQ_hygienic_i][PSLQ_hygienic_m+1]; \
        PSLQ_hygienic_h[PSLQ_hygienic_i][PSLQ_hygienic_m] = PSLQ_hygienic_t1*PSLQ_hygienic_t3 + PSLQ_hygienic_t2*PSLQ_hygienic_t4; \
        PSLQ_hygienic_h[PSLQ_hygienic_i][PSLQ_hygienic_m+1] = PSLQ_hygienic_t1*PSLQ_hygienic_t4 - PSLQ_hygienic_t2*PSLQ_hygienic_t3; \
      } \
    } \
 \
    for(PSLQ_hygienic_i=PSLQ_hygienic_m+1;PSLQ_hygienic_i<(n);PSLQ_hygienic_i++) { \
      for(PSLQ_hygienic_j=MIN(PSLQ_hygienic_i-1,PSLQ_hygienic_m+1); PSLQ_hygienic_j >=0; PSLQ_hygienic_j--) { \
        PSLQ_hygienic_t = __builtin_nearbyint(PSLQ_hygienic_h[PSLQ_hygienic_i][PSLQ_hygienic_j]/PSLQ_hygienic_h[PSLQ_hygienic_j][PSLQ_hygienic_j]); \
        PSLQ_hygienic_y[PSLQ_hygienic_j] += PSLQ_hygienic_t*PSLQ_hygienic_y[PSLQ_hygienic_i]; \
        for(PSLQ_hygienic_k=0;PSLQ_hygienic_k<=PSLQ_hygienic_j;PSLQ_hygienic_k++) { PSLQ_hygienic_h[PSLQ_hygienic_i][PSLQ_hygienic_k] -= PSLQ_hygienic_t*PSLQ_hygienic_h[PSLQ_hygienic_j][PSLQ_hygienic_k]; } \
        for(PSLQ_hygienic_k=0;PSLQ_hygienic_k<(n);PSLQ_hygienic_k++) { \
          PSLQ_hygienic_a[PSLQ_hygienic_i][PSLQ_hygienic_k] -= PSLQ_hygienic_t*PSLQ_hygienic_a[PSLQ_hygienic_j][PSLQ_hygienic_k]; \
          PSLQ_hygienic_b[PSLQ_hygienic_j][PSLQ_hygienic_k] += PSLQ_hygienic_t*PSLQ_hygienic_b[PSLQ_hygienic_i][PSLQ_hygienic_k]; \
        } \
      } \
    } \
 \
    PSLQ_hygienic_m_val = -1.0e308; \
    for(PSLQ_hygienic_j=0;PSLQ_hygienic_j<(n)-1;PSLQ_hygienic_j++) { \
      PSLQ_hygienic_t = ABS( PSLQ_hygienic_h[PSLQ_hygienic_j][PSLQ_hygienic_j] ); \
      if (PSLQ_hygienic_t > PSLQ_hygienic_m_val) PSLQ_hygienic_m_val = PSLQ_hygienic_t; \
    } \
 \
    for(PSLQ_hygienic_i=0;PSLQ_hygienic_i<(n);PSLQ_hygienic_i++) { \
      PSLQ_hygienic_t = ABS( PSLQ_hygienic_y[PSLQ_hygienic_i] ); \
      if (PSLQ_hygienic_t<PSLQ_hygienic_teps) { \
          PSLQ_hygienic_m = PSLQ_hygienic_i; \
          PSLQ_hygienic_process = 0; \
          break; \
      } \
    } \
  } \
 \
  for(PSLQ_hygienic_i=0;PSLQ_hygienic_i<(n);PSLQ_hygienic_i++) { (r)[PSLQ_hygienic_i] = PSLQ_hygienic_b[PSLQ_hygienic_m][PSLQ_hygienic_i]; } \
}
