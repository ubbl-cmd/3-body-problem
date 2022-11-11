#include <iostream>
#include <math.h>

int size = 4;
double* a;
double* b;

int n;
double* cb_r;
double* cb_v;
double* cb_m;

double G = 6.67430e-20;

double veclen(double * v) {
	double res = 0;
	for (int i = 0; i < 3; i++) {
		res += v[i] * v[i];
	}
	return sqrt(res);
}

void f(double* x, double *res) {
	double* r = new double[3];
    for (int i = 0; i < n * 3; i++) {
        res[i] = x[i + n * 3];
    }
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (j == i) {
				continue;
			}
			for (int idx = 0; idx < 3; idx++) {
				r[idx] = x[j * 3 + idx] - x[3 * i + idx];
			}
			for (int idx = 0; idx < 3; idx++) {
				res[3 * i + idx + 3 * n] += r[idx] * (G * cb_m[j] / pow(veclen(r), 3));
			}
		}
	}
	delete[] r;
}


void rk(double h) {
	double *k = new double[size * 2 * n * 3];
	double * x = new double[2 * n * 3];
	for (int ki = 0; ki < size; ki++) {
		for (int i = 0; i < n * 3; i++) {
			x[i] = cb_r[i];
			for (int kis = 0; kis < ki; kis++) {
                x[i] += h *     a[ki * size + kis] * k[kis * 2 * n * 3 + i];
				x[i] += h * h * a[ki * size + kis] * k[kis * 2 * n * 3 + i + n * 3];
			}
		}
        for (int i = 0; i < n * 3; i++) {
            x[i + n * 3] = cb_v[i];
            for (int kis = 0; kis < ki; kis++) {
                x[i + n * 3] += h * h * a[ki * size + kis] * k[kis * 2 * n * 3 + i + n * 3];
            }
        }
		double *dv = new double[2 * n * 3];
		f(x, dv);
		for (int i = 0; i < 2 * n * 3; i++) {
			k[ki * 2 * n * 3 + i] = dv[i];
		}
		delete[] dv;
	}
	for (int ki = 0; ki < size; ki++) {
		for (int i = 0; i < n * 3; i++) {
			cb_v[i] += h * b[ki] * k[ki * 2 * n * 3 + i + n * 3];
	    }
    }
	for (int i = 0; i < n * 3; i++) {
		cb_r[i] += h * cb_v[i];
	}
	delete[] k;
}

int main() {
	a = new double[size * size]{
		0,   0,   0, 0,
		0.5, 0,   0, 0,
		0,   0.5, 0, 0,
		0,   0,   1, 0
	};
	b = new double[size] {
		1./6, 1./3, 1./3, 1./6
	};
/*
----------------------------------------------
	START INIT CELESTIAL OBJECTS
----------------------------------------------*/
	n = 3;
//    n = 2;
	cb_r = new double [n * 3] {
		 1.453004563709436E+08, 2.978607390059390E+07, 2.866143171530403E+04,
		1.455259266863196E+08, 2.949542371875033E+07, -4.639949931440875E+03,
		-1.359757021307868E+06, 1.333076360298289E+05, 3.057372007160987E+04

	};
	cb_v = new double [n * 3] {
		-6.397102954570147E+00,  2.906192414721306E+01, -1.203250030542335E-03,
		-5.557389967864156E+00, 2.971153423303235E+01, -1.096802103212369E-02,
		-2.335983427600637E-04, -1.571626050150792E-02,  1.314580230199784E-04
	};
	cb_m = new double [n]{
		5.9722e24,
		7.349e22,
		1.98847e30
	};
/*
----------------------------------------------
	ENDOF INIT CELESTIAL OBJECTS
----------------------------------------------*/

	double h = 60;
	double t = 0;
	double T;
	T = 31536000 / 2;
	// T=100;
	int i = 0;
	while (t < T) {
//		if (i % 500 == 1) {
			std::cout << t << " ";
			for (int i = 0; i < n * 3; i++) {
				std::cout << cb_r[i] << " ";
			}
			std::cout << "\n";
//		}
		rk(h);

		i++;
		t += h;
	}

	delete[] cb_r;
	delete[] cb_v;
	delete[] cb_m;

	delete[] a;
	delete[] b;
	return 0;
}
