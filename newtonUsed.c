#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <time.h>

#define ZEROCHARVAL 48
#define TOLSQ 1e-6
#define ABSSQ(zfl) (creal(zfl)*creal(zfl) + cimag(zfl)*cimag(zfl))
#define SZ 1000000lu
#define IX 4000

void newton_iter(const double re_z0, const double im_z0, const char *degree_ptr, char *attr_indices, size_t *n_iter) {
  double realdum, imagdum;
  complex double zDum;
  // Initial guess, initialize iteration counter
  complex double zVal = re_z0 + im_z0*I;
  *n_iter = 0;
    
  // Commence spaghetti
  switch (*degree_ptr) {
    case '1':
      *attr_indices = '1';
      zVal = (complex double) 1.;
      *n_iter = 1;
      break;
    case '2':
      for (;;) {
        realdum = fabs(creal(zVal));
        imagdum = fabs(cimag(zVal));
        if ((realdum > 1e10) || (imagdum > 1e10)) {
          *attr_indices = '0';
          break;
        }
        if (realdum*realdum + imagdum*imagdum < TOLSQ) {
          *attr_indices = '0';
          break;
        }
        if (ABSSQ(zVal*zVal - 1.) < 4e-6) {
          if (ABSSQ(zVal - 1.) < TOLSQ) {
            *attr_indices = '1';
            break;
          }
          if (ABSSQ(zVal + 1.) < TOLSQ) {
            *attr_indices = '2';
            break;
          }
        }
        zVal = (zVal + 1./zVal)/2.;
        (*n_iter)++;
      }
      break;
    case '3':
      for (;;) {
        realdum = fabs(creal(zVal));
        imagdum = fabs(cimag(zVal));
        if ((realdum > 1e10) || (imagdum > 1e10)) {
          *attr_indices = '0';
          break;
        }
        if (realdum*realdum + imagdum*imagdum < TOLSQ) {
          *attr_indices = '0';
          break;
        }
        zDum = zVal*zVal;
        if (ABSSQ(zDum*zVal - 1.) < 1.6e-5) {
          if (ABSSQ(zVal - 1.) < TOLSQ) {
            *attr_indices = '1';
            break;
          }
          if (ABSSQ(zVal + (0.5 - 0.8660254037844386*I)) < TOLSQ) {
            *attr_indices = '2';
            break;
          }
          if (ABSSQ(zVal + (0.5 + 0.8660254037844386*I)) < TOLSQ) {
            *attr_indices = '3';
            break;
          }
        }
        zVal = (2.*zVal + 1./zDum)/3.;
        (*n_iter)++;
      }
      break;
    case '4':
      for (;;) {
        realdum = fabs(creal(zVal));
        imagdum = fabs(cimag(zVal));
        if ((realdum > 1e10) || (imagdum > 1e10)) {
          *attr_indices = '0';
          break;
        }
        if (realdum*realdum + imagdum*imagdum < TOLSQ) {
          *attr_indices = '0';
          break;
        }
        zDum = zVal*zVal;
        if (ABSSQ(zDum*zDum - 1.) < 6.4e-5) {
          if (ABSSQ(zVal - 1.) < TOLSQ) {
            *attr_indices = '1';
            break;
          }
          if (ABSSQ(zVal - I) < TOLSQ) {
            *attr_indices = '2';
            break;
          }
          if (ABSSQ(zVal + 1.) < TOLSQ) {
            *attr_indices = '3';
            break;
          }
          if (ABSSQ(zVal + I) < TOLSQ) {
            *attr_indices = '4';
            break;
          }
        }
        zVal = (3.*zVal + 1./(zDum*zVal))/4.;
        (*n_iter)++;
      }
      break;
    case '5':
      for (;;) {
        realdum = fabs(creal(zVal));
        imagdum = fabs(cimag(zVal));
        if ((realdum > 1e10) || (imagdum > 1e10)) {
          *attr_indices = '0';
          break;
        }
        if (realdum*realdum + imagdum*imagdum < TOLSQ) {
          *attr_indices = '0';
          break;
        }
        zDum = zVal*zVal;
        zDum *= zDum;
        if (ABSSQ(zDum*zVal - 1.) < 2.56e-4) {
          if (ABSSQ(zVal - 1.) < TOLSQ) {
            *attr_indices = '1';
            break;
          }
          if (ABSSQ(zVal - (0.309016994374948 + 0.951056516295153*I)) < TOLSQ) {
            *attr_indices = '2';
            break;
          }
          if (ABSSQ(zVal + 0.809016994374947 - 0.587785252292473*I) < TOLSQ) {
            *attr_indices = '3';
            break;
          }
          if (ABSSQ(zVal + 0.809016994374947 + 0.587785252292473*I) < TOLSQ) {
            *attr_indices = '4';
            break;
          }
          if (ABSSQ(zVal - 0.309016994374948 + 0.951056516295153*I) < TOLSQ) {
            *attr_indices = '5';
            break;
          }
        }
        zVal = (4.*zVal + 1./zDum)/5.;
        (*n_iter)++;
      }
      break;
    case '6':
      for (;;) {
        realdum = fabs(creal(zVal));
        imagdum = fabs(cimag(zVal));
        if ((realdum > 1e10) || (imagdum > 1e10)) {
          *attr_indices = '0';
          break;
        }
        if (realdum*realdum + imagdum*imagdum < TOLSQ) {
          *attr_indices = '0';
          break;
        }
        zDum = zVal*zVal*zVal;
        if (ABSSQ(zDum*zDum - 1.) < 1.024e-3) {
          if (ABSSQ(zVal - 1.) < TOLSQ) {
            *attr_indices = '1';
            break;
          }
          if (ABSSQ(zVal - (0.5 + 0.8660254037844386*I)) < TOLSQ) {
            *attr_indices = '2';
            break;
          }
          if (ABSSQ(zVal + 0.5 - 0.8660254037844386*I) < TOLSQ) {
            *attr_indices = '3';
            break;
          }
          if (ABSSQ(zVal + 1.) < TOLSQ) {
            *attr_indices = '4';
            break;
          }
          if (ABSSQ(zVal + 0.5 + 0.8660254037844386*I) < TOLSQ) {
            *attr_indices = '5';
            break;
          }
          if (ABSSQ(zVal - 0.5 + 0.8660254037844386*I) < TOLSQ) {
            *attr_indices = '6';
            break;
          }
        }
        zVal = (5.*zVal + 1./(zDum*zVal*zVal))/6.;
        (*n_iter)++;
      }
      break;
    case '7':
      for (;;) {
        realdum = fabs(creal(zVal));
        imagdum = fabs(cimag(zVal));
        if ((realdum > 1e10) || (imagdum > 1e10)) {
          *attr_indices = '0';
          break;
        }
        if (realdum*realdum + imagdum*imagdum < TOLSQ) {
          *attr_indices = '0';
          break;
        }
        zDum = zVal*zVal*zVal;
        if (ABSSQ(zDum*zDum*zVal - 1.) < 4.096e-3) {
          if (ABSSQ(zVal - 1.) < TOLSQ) {
            *attr_indices = '1';
            break;
          }
          if (ABSSQ(zVal - (0.623489801858734 + 0.781831482468030*I)) < TOLSQ) {
            *attr_indices = '2';
            break;
          }
          if (ABSSQ(zVal + 0.222520933956314 - 0.974927912181824*I) < TOLSQ) {
            *attr_indices = '3';
            break;
          }
          if (ABSSQ(zVal + 0.900968867902419 - 0.433883739117558*I) < TOLSQ) {
            *attr_indices = '4';
            break;
          }
          if (ABSSQ(zVal + 0.900968867902419 + 0.433883739117558*I) < TOLSQ) {
            *attr_indices = '5';
            break;
          }
          if (ABSSQ(zVal + 0.222520933956314 + 0.974927912181824*I) < TOLSQ) {
            *attr_indices = '6';
            break;
          }
          if (ABSSQ(zVal - 0.623489801858734 + 0.781831482468030*I) < TOLSQ) {
            *attr_indices = '7';
            break;
          }
        }
        zVal = (6.*zVal + 1./(zDum*zDum))/7.;
        (*n_iter)++;
      }
      break;
    case '8':
      for (;;) {
        realdum = fabs(creal(zVal));
        imagdum = fabs(cimag(zVal));
        if ((realdum > 1e10) || (imagdum > 1e10)) {
          *attr_indices = '0';
          break;
        }
        if (realdum*realdum + imagdum*imagdum < TOLSQ) {
          *attr_indices = '0';
          break;
        }
        zDum = zVal*zVal;
        zDum *= zDum;
        if (ABSSQ(zDum*zDum - 1.) < 1.6384e-2) {
          if (ABSSQ(zVal - 1.) < TOLSQ) {
            *attr_indices = '1';
            break;
          }
          if (ABSSQ(zVal - (0.707106781186548 + 0.707106781186547*I)) < TOLSQ) {
            *attr_indices = '2';
            break;
          }
          if (ABSSQ(zVal - I) < TOLSQ) {
            *attr_indices = '3';
            break;
          }
          if (ABSSQ(zVal + 0.707106781186548 - 0.707106781186547*I) < TOLSQ) {
            *attr_indices = '4';
            break;
          }
          if (ABSSQ(zVal + 1.) < TOLSQ) {
            *attr_indices = '5';
            break;
          }
          if (ABSSQ(zVal + 0.707106781186548 + 0.707106781186547*I) < TOLSQ) {
            *attr_indices = '6';
            break;
          }
          if (ABSSQ(zVal + I) < TOLSQ) {
            *attr_indices = '7';
            break;
          }
          if (ABSSQ(zVal - 0.707106781186548 + 0.707106781186547*I) < TOLSQ) {
            *attr_indices = '8';
            break;
          }
        }
        zVal = (7.*zVal + 1./(zDum*zVal*zVal*zVal))/8.;
        (*n_iter)++;
      }
      break;
    case '9':
      for (;;) {
        realdum = fabs(creal(zVal));
        imagdum = fabs(cimag(zVal));
        if ((realdum > 1e10) || (imagdum > 1e10)) {
          *attr_indices = '0';
          break;
        }
        if (realdum*realdum + imagdum*imagdum < TOLSQ) {
          *attr_indices = '0';
          break;
        }
        zDum = zVal*zVal;
        zDum *= zDum;
        if (ABSSQ(zDum*zDum*zVal - 1.) < 6.5536e-2) {
          if (ABSSQ(zVal - 1.) < TOLSQ) {
            *attr_indices = '1';
            break;
          }
          if (ABSSQ(zVal - (0.766044443118978 + 0.642787609686539*I)) < TOLSQ) {
            *attr_indices = '2';
            break;
          }
          if (ABSSQ(zVal - (0.173648177666930 + 0.984807753012208*I)) < TOLSQ) {
            *attr_indices = '3';
            break;
          }
          if (ABSSQ(zVal + 0.500000000000000 - 0.866025403784439*I) < TOLSQ) {
            *attr_indices = '4';
            break;
          }
          if (ABSSQ(zVal + 0.939692620785908 - 0.342020143325669*I) < TOLSQ) {
            *attr_indices = '5';
            break;
          }
          if (ABSSQ(zVal + 0.939692620785908 + 0.342020143325669*I) < TOLSQ) {
            *attr_indices = '6';
            break;
          }
          if (ABSSQ(zVal + 0.500000000000000 + 0.866025403784439*I) < TOLSQ) {
            *attr_indices = '7';
            break;
          }
          if (ABSSQ(zVal - 0.173648177666930 + 0.984807753012208*I) < TOLSQ) {
            *attr_indices = '8';
            break;
          }
          if (ABSSQ(zVal - 0.766044443118978 + 0.642787609686539*I) < TOLSQ) {
            *attr_indices = '9';
            break;
          }
        }
        zVal = (8.*zVal + 1./(zDum*zDum))/9.;
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
    const char degree = '5';
    char *attractor = (char*) malloc(SZ*sizeof(char));
    size_t *convergence = (size_t*) malloc(SZ*sizeof(size_t));

    for ( int jx = 0; jx < SZ; ++jx ) {
        newton_iter(-2. + 4./(SZ - 1)*IX,  2. - 4./(SZ - 1)*jx, &degree, &attractor[jx], &convergence[jx]);
    }

    printf("We find zero #%c of x^%d - 1 in %zu iteration(s) when our initial guess is %f + %fi\n", attractor[11], degree - ZEROCHARVAL, convergence[11], -2. + 4./(SZ - 1)*IX, 2. - 4./(SZ - 1)*11);
    //printf("We find zero #%c of x^%d - 1 in %zu iteration(s) when our initial guess is %f + %fi\n", attractor[SZ/2], degree - ZEROCHARVAL, convergence[SZ/2], -2. + 4./(SZ - 1)*IX, 2. - 4./(SZ - 1)*SZ/2);
    //printf("(vårt värde är %f + %fi)\n", crealf(zVal), cimagf(zVal));

    free(attractor);
    free(convergence);
    return 0;
}