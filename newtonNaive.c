#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <time.h>

#define TOL 1e-3
#define SZ 1000000lu
#define IX 400000

void newton_iter(double re_z0, double im_z0, int *degree_ptr, int *attr_indices, size_t *n_iter) {
  // Initial guess, initialize iteration counter
  complex double zVal = re_z0 + im_z0*I;
  *n_iter = 0;
    
  // Commence spaghetti
  switch (*degree_ptr) {
    case 1:
        for (;;) {
            if ((fabs(creal(zVal)) > 1e10) || (imagfabs(cimag(zVal)) > 1e10)) {
                *attr_indices = 0;
                break;
            }
            if (cabs(zVal) < TOL) {
                *attr_indices = 0;
                break;
            }
            if (cabs(zVal - 1.) < TOL) {
                *attr_indices = 1;
                break;
            }
            zVal = (complex double) 1.;
            (*n_iter)++;
        }
        break;
    case 2:
        for (;;) {
            if ((fabs(creal(zVal)) > 1e10) || (fabs(cimag(zVal)) > 1e10)) {
                *attr_indices = 0;
                break;
            }
            if (cabs(zVal) < TOL) {
                *attr_indices = 0;
                break;
            }
            if (cabs(zVal - 1.) < TOL) {
                *attr_indices = 1;
                break;
            }
            if (cabs(zVal + 1.) < TOL) {
                *attr_indices = 2;
                break;
            }
            zVal = (zVal + 1./zVal)/2.;
            (*n_iter)++;
        }
        break;
    case 3:
        for (;;) {
            if ((fabs(creal(zVal)) > 1e10) || (fabs(cimag(zVal)) > 1e10)) {
                *attr_indices = 0;
                break;
            }
            if (cabs(zVal) < TOL) {
                *attr_indices = 0;
                break;
            }
            if (cabs(zVal - 1.) < TOL) {
                *attr_indices = 1;
                break;
            }
            if (cabs(zVal + (0.5 - 0.8660254037844386*I)) < TOL) {
                *attr_indices = 2;
                break;
            }
            if (cabs(zVal + (0.5 + 0.8660254037844386*I)) < TOL) {
                *attr_indices = 3;
                break;
            }
            zVal = (2.*zVal + 1./(zVal*zVal))/3.;
            (*n_iter)++;
        }
        break;
    case 4:
        for (;;) {
            if ((fabs(creal(zVal)) > 1e10) || (fabs(cimag(zVal)) > 1e10)) {
                *attr_indices = 0;
                break;
            }
            if (cabs(zVal) < TOL) {
                *attr_indices = 0;
                break;
            }
            if (cabs(zVal - 1.) < TOL) {
                *attr_indices = 1;
                break;
            }
            if (cabs(zVal - I) < TOL) {
                *attr_indices = 2;
                break;
            }
            if (cabs(zVal + 1.) < TOL) {
                *attr_indices = 3;
                break;
            }
            if (cabs(zVal + I) < TOL) {
                *attr_indices = 4;
                break;
            }
            zVal = (3.*zVal + 1./(zVal*zVal*zVal))/4.;
            (*n_iter)++;
        }
        break;
    case 5:
        for (;;) {
            if ((fabs(creal(zVal)) > 1e10) || (fabs(cimag(zVal)) > 1e10)) {
                *attr_indices = 0;
                break;
            }
            if (cabs(zVal) < TOL) {
                *attr_indices = 0;
                break;
            }
            if (cabs(zVal - 1.) < TOL) {
                *attr_indices = 1;
                break;
            }
            if (cabs(zVal - (0.309016994374948 + 0.951056516295153*I)) < TOL) {
                *attr_indices = 2;
                break;
            }
            if (cabs(zVal + 0.809016994374947 - 0.587785252292473*I) < TOL) {
                *attr_indices = 3;
                break;
            }
            if (cabs(zVal + 0.809016994374947 + 0.587785252292473*I) < TOL) {
                *attr_indices = 4;
                break;
            }
            if (cabs(zVal - 0.309016994374948 + 0.951056516295153*I) < TOL) {
                *attr_indices = 5;
                break;
            }
            zVal = (4.*zVal + 1./(zVal*zVal*zVal*zVal))/5.;
            (*n_iter)++;
        }
        break;
    case 6:
        for (;;) {
            if ((fabs(creal(zVal)) > 1e10) || (fabs(cimag(zVal)) > 1e10)) {
                *attr_indices = 0;
                break;
            }
            if (cabs(zVal) < TOL) {
                *attr_indices = 0;
                break;
            }
            if (cabs(zVal - 1.) < TOL) {
                *attr_indices = 1;
                break;
            }
            if (cabs(zVal - (0.5 + 0.8660254037844386*I)) < TOL) {
                *attr_indices = 2;
                break;
            }
            if (cabs(zVal + 0.5 - 0.8660254037844386*I) < TOL) {
                *attr_indices = 3;
                break;
            }
            if (cabs(zVal + 1.) < TOL) {
                *attr_indices = 4;
                break;
            }
            if (cabs(zVal + 0.5 + 0.8660254037844386*I) < TOL) {
                *attr_indices = 5;
                break;
            }
            if (cabs(zVal - 0.5 + 0.8660254037844386*I) < TOL) {
                *attr_indices = 6;
                break;
            }
            zVal = (5.*zVal + 1./(zVal*zVal*zVal*zVal*zVal))/6.;
            (*n_iter)++;
        }
        break;
    case 7:
        for (;;) {
            if ((fabs(creal(zVal)) > 1e10) || (fabs(cimag(zVal)) > 1e10)) {
                *attr_indices = 0;
                break;
            }
            if (cabs(zVal) < TOL) {
                *attr_indices = 0;
                break;
            }
            if (cabs(zVal - 1.) < TOL) {
                *attr_indices = 1;
                break;
            }
            if (cabs(zVal - (0.623489801858734 + 0.781831482468030*I)) < TOL) {
                *attr_indices = 2;
                break;
            }
            if (cabs(zVal + 0.222520933956314 - 0.974927912181824*I) < TOL) {
                *attr_indices = 3;
                break;
            }
            if (cabs(zVal + 0.900968867902419 - 0.433883739117558*I) < TOL) {
                *attr_indices = 4;
                break;
            }
            if (cabs(zVal + 0.900968867902419 + 0.433883739117558*I) < TOL) {
                *attr_indices = 5;
                break;
            }
            if (cabs(zVal + 0.222520933956314 + 0.974927912181824*I) < TOL) {
                *attr_indices = 6;
                break;
            }
            if (cabs(zVal - 0.623489801858734 + 0.781831482468030*I) < TOL) {
                *attr_indices = 7;
                break;
            }
            zVal = (6.*zVal + 1./(zVal*zVal*zVal*zVal*zVal*zVal))/7.;
            (*n_iter)++;
        }
        break;
    case 8:
        for (;;) {
            if ((fabs(creal(zVal)) > 1e10) || (fabs(cimag(zVal)) > 1e10)) {
                *attr_indices = 0;
                break;
            }
            if (cabs(zVal) < TOL) {
                *attr_indices = 0;
                break;
            }
            if (cabs(zVal - 1.) < TOL) {
                *attr_indices = 1;
                break;
            }
            if (cabs(zVal - (0.707106781186548 + 0.707106781186547*I)) < TOL) {
                *attr_indices = 2;
                break;
            }
            if (cabs(zVal - I) < TOL) {
                *attr_indices = 3;
                break;
            }
            if (cabs(zVal + 0.707106781186548 - 0.707106781186547*I) < TOL) {
                *attr_indices = 4;
                break;
            }
            if (cabs(zVal + 1.) < TOL) {
                *attr_indices = 5;
                break;
            }
            if (cabs(zVal + 0.707106781186548 + 0.707106781186547*I) < TOL) {
                *attr_indices = 6;
                break;
            }
            if (cabs(zVal + I) < TOL) {
                *attr_indices = 7;
                break;
            }
            if (cabs(zVal - 0.707106781186548 + 0.707106781186547*I) < TOL) {
                *attr_indices = 8;
                break;
            }
            zVal = (7.*zVal + 1./(zVal*zVal*zVal*zVal*zVal*zVal*zVal))/8.;
            (*n_iter)++;
        }
        break;
    case 9:
        for (;;) {
            if ((fabs(creal(zVal)) > 1e10) || (fabs(cimag(zVal)) > 1e10)) {
                *attr_indices = 0;
                break;
            }
            if (cabs(zVal) < TOL) {
                *attr_indices = 0;
                break;
            }
            if (cabs(zVal - 1.) < TOL) {
                *attr_indices = 1;
                break;
            }
            if (cabs(zVal - (0.766044443118978 + 0.642787609686539*I)) < TOL) {
                *attr_indices = 2;
                break;
            }
            if (cabs(zVal - (0.173648177666930 + 0.984807753012208*I)) < TOL) {
                *attr_indices = 3;
                break;
            }
            if (cabs(zVal + 0.500000000000000 - 0.866025403784439*I) < TOL) {
                *attr_indices = 4;
                break;
            }
            if (cabs(zVal + 0.939692620785908 - 0.342020143325669*I) < TOL) {
                *attr_indices = 5;
                break;
            }
            if (cabs(zVal + 0.939692620785908 + 0.342020143325669*I) < TOL) {
                *attr_indices = 6;
                break;
            }
            if (cabs(zVal + 0.500000000000000 + 0.866025403784439*I) < TOL) {
                *attr_indices = 7;
                break;
            }
            if (cabs(zVal - 0.173648177666930 + 0.984807753012208*I) < TOL) {
                *attr_indices = 8;
                break;
            }
            if (cabs(zVal - 0.766044443118978 + 0.642787609686539*I) < TOL) {
                *attr_indices = 9;
                break;
            }
            zVal = (8.*zVal + 1./(zVal*zVal*zVal*zVal*zVal*zVal*zVal*zVal))/9.;
            (*n_iter)++;
        }
        break;
    // No further cases. Hyyyyyype.

    default:
      fprintf(stderr, "unexpected degree\n");
      exit(1);
  }
}

int main() {
    int degree = 5;
    int *attractor = (char*) malloc(SZ*sizeof(int));
    size_t *convergence = (size_t*) malloc(SZ*sizeof(size_t));

    for ( int jx = 0; jx < SZ; ++jx ) {
        newton_iter(-2. + 4./(SZ - 1)*IX,  2. - 4./(SZ - 1)*jx, &degree, &attractor[jx], &convergence[jx]);
    }

    printf("We find zero #%d of x^%d - 1 in %zu iteration(s) when our initial guess is %f + %fi\n", attractor[11], degree, convergence[11], -2. + 4./(SZ - 1)*IX, 2. - 4./(SZ - 1)*11);
    //printf("We find zero #%c of x^%d - 1 in %zu iteration(s) when our initial guess is %f + %fi\n", attractor[SZ/2], degree - ZEROCHARVAL, convergence[SZ/2], -2. + 4./(SZ - 1)*IX, 2. - 4./(SZ - 1)*SZ/2);
    //printf("(vårt värde är %f + %fi)\n", crealf(zVal), cimagf(zVal));

    free(attractor);
    free(convergence);
    return 0;
}