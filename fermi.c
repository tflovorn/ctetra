#include "fermi.h"

int FindFermi(int n, int num_bands, double num_electrons, EnergyCache *Ecache, double *E_Fermi) {
    time_t start_time = time(NULL);
    double statecount_error(double E) {
        double count = NumStates(E, n, num_bands, Ecache);
        time_t this_time = time(NULL);
        printf("At E = %f got count = %f; time = %f\n", E, count, difftime(this_time, start_time));

        return count - num_electrons;
    }
    // TODO replace this with cached calls if required.
    double emin, emax;
    MinMaxVals(n, num_bands, Ecache, &emin, &emax);

    double toln = 1e-10;
    double maxiter = 300;
    int err = Bisect(E_Fermi, statecount_error, emin, emax, toln, maxiter);
    return err;
}

int Bisect(double *result, BisectFn fn, double a, double b, double toly, int maxiter) {
	double val_a = fn(a);
	double val_b = fn(b);
	if (val_a == 0.0) {
	    *result = a;
        return CTETRA_BISECT_OK;
	} else if (val_b == 0.0) {
        *result = b;
        return CTETRA_BISECT_OK;
	}
	if ((val_a < 0.0 && val_b < 0.0) || (val_a > 0.0 && val_b > 0.0)) {
		return CTETRA_BISECT_BRACKET;
	}

    double low = a;
    double high = b;
	if (val_b < 0.0) {
		low = b;
        high = a;
	}

    int iter;
    double mid, val_mid;
	for (iter = 0; iter < maxiter; iter++) {
		mid = (low + high) / 2.0;
		val_mid = fn(mid);

		if (fabs(val_mid) < toly) {
            *result = mid;
            return CTETRA_BISECT_OK;
		} else if (val_mid > 0.0) {
			// go toward low
			high = mid;
		} else {
			// go toward high
			low = mid;
		}
	}
    return CTETRA_BISECT_MAXITER;
}
