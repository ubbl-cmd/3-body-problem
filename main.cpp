/*
TODO:
 Методы:
 [*] РК4
 [*] ДП8
 [*] явный метод Адамса 8-го порядка с разгоном ДП8
 [*] (КАК ОЧИСТИТЬ ПАМЯТЬ)
 Инварианты:
 [ ] положение и скорость барицентра
 [ ] момент импульса
 [ ] энергия (кинетическая + потенциальная)
 Проверка:
 [ ] Сравнение с аналитическим решением для случая двух тел (Кеплерова орбита)
 */



#include <iostream>
#include <math.h>

int size;
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
	double *k  = new double[size * 2 * n * 3];
	double *dx = new double[       2 * n * 3];
	for (int i = 0; i < 2 * n * 3; i++) {
		k[i] = 0;
		dx[i] = 0;
		// dv[i] = 0;
	}
	double *dv = new double[       2 * n * 3];
	for (int ki = 0; ki < size; ki++) {
		for (int i = 0; i < n * 3; i++) {
			dx[i]         = cb_r[i];
			dx[i + n * 3] = cb_v[i];
			for (int kis = 0; kis < ki; kis++) {
				dx[i]         += h * a[ki * size + kis] * k[kis * 2 * n * 3 + i];
				dx[i + n * 3] += h * a[ki * size + kis] * k[kis * 2 * n * 3 + i + n * 3];
			}

		}
		for (int i = 0; i < 2*n * 3; i++) {
			dv[i] = 0;
		}
		f(dx, dv);
		for (int i = 0; i < 2 * n * 3; i++) {
			k[ki * 2 * n * 3 + i] = dv[i];
		}
	}
	for (int ki = 0; ki < size; ki++) {
		for (int i = 0; i < n * 3; i++) {
			cb_r[i] += h * b[ki] * k[ki * 2 * n * 3 + i];
		}
	}
	for (int ki = 0; ki < size; ki++) {
		for (int i = 0; i < n * 3; i++) {
			cb_v[i] += h * b[ki] * k[ki * 2 * n * 3 + i + n * 3];
		}
	}
	delete[] dv;
	delete[] dx;
	delete[] k;
}

double *adams_dx;
double *adams_a;
void adams(double h) {
	for (int ki = 0; ki < size; ki++) {
		for (int idx = 0; idx < n * 3; idx++) {
			cb_r[idx] += h * adams_a[ki] * adams_dx[ki * 2 * n * 3 + idx];
			cb_v[idx] += h * adams_a[ki] * adams_dx[ki * 2 * n * 3 + idx + n*3];
		}
	}
	for (int i = 0; i < (size-1) * 2 * n * 3; i++) {
		adams_dx[i] = adams_dx[i+2*n*3];
	}
	double *dx = new double[2*n*3];
	for (int i = 0; i < n*3; i++) {
		dx[i    ] = cb_r[i];
		dx[i+n*3] = cb_v[i];
	}
	double *res = new double[2*n*3];
	for (int idx = 0; idx < 2*n*3; idx++) {
		res[idx] = 0;
	}
	f(dx, res);
	for (int idx = 0; idx < 2*n*3; idx++) {
		adams_dx[(size-1) * 2*n*3 + idx] = res[idx];
	}
	delete[] res;
	delete[] dx;
}

int main() {
/*
----------------------------------------------
	START INIT CELESTIAL OBJECTS
----------------------------------------------*/
	n = 3;
   // n = 2;
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
	double h = 50;
	double t = 0;
	double T;
	T = 31536000 / 4;
	// T=10000;
	//eu
	// size = 1;
	// a = new double[size * size] {
	// 	0
	// };
	// b = new double[size] {
	// 	1.
	// };
	//rk4
	// size = 4;
	// a = new double[size * size] {
	// 	0,   0,   0, 0,
	// 	0.5, 0,   0, 0,
	// 	0,   0.5, 0, 0,
	// 	0,   0,   1, 0
	// };
	// b = new double[size] {
	// 	1./6, 1./3, 1./3, 1./6
	// };
	// //dp8

	size = 13;
	a = new double[size * size] {

		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,

		1./18, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,

		1./48, 1.16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,

		1./32, 0, 3./32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,

		5.16, 0, -75./64, 75./64, 0, 0, 0, 0, 0, 0, 0, 0, 0,

		3./80, 0, 0, 3./16, 3./20, 0, 0, 0, 0, 0, 0, 0, 0,

		29443841./614563906, 0, 0, 77736538./6142538347, -28693883./1112000000, 23124283393./1800000000, 0, 0, 0, 0, 0, 0, 0,

		16016141./946692911, 0, 0, 61564180./158732637, 22789713./633445777, 545815736./2771057229, -180193667./1043307555, 0, 0, 0, 0, 0, 0,

		39632708./573591083, 0, 0, -433636366./683701615, -421739975./2616292301, 100302831./723423059, 790204164./839813087, 800635310./3783071287, 0, 0, 0, 0, 0,

		246121993./1340847787, 0, 0, -37695042795./15268766246, -309121744./1061227803, -12992083./490766935, 6005943433./2108947869, 393006217./1396673457, 123782331./1001029789, 0, 0, 0, 0,

		-1028468189./846180014, 0, 0, 8478235783./506512852, 1311729495./1432422823, -10304129995./1701304382, -48777925059./3047939560, 15336726248./1032824649, -45442868181./3398467696, 3065993473./597172653, 0, 0, 0,

		185892177./718116043, 0, 0, -3185094517./667107341, -477755414./1098053517, -703635378./230739211, 5731566787./1027545527, 5232866602./850066563, -4093664535./808688257, 3962137247./1805957418, 65686385./487910083, 0, 0,
		403863854./491063109, 0, 0, -5068492393./434740067, -411421997./543043805, 652783627./914296604, 11173962825./925320556, -13158990841./6184727034, 3936647629./1978049680, -160528059./685178525, 248638103./1413531060, 0, 0
	};
	b = new double[size] {
		14005451./335480064, 0, 0, 0, 0, -59238493./1068277825, 181606676./758867731, 561292985./797845732, -1041891430./1371343529, 760417239./1151165299, 118820643./751138087, -5228747749./2220607170, 1/4
	};

/*
----------------------------------------------
	ADAMS 8 + DORPRI8
----------------------------------------------*/
	adams_dx = new double[8 * 2 * n * 3];
	for (int i = 0; i < 8; i++) {
		rk(h);
		std::cout << t << " ";
		for (int idx = 0; idx < n * 3; idx++) {
			std::cout << cb_r[idx] << " ";
		}
		std::cout << "\n";
		double *dx = new double[2*n*3];
		for (int idx = 0; idx < n*3; idx++) {
			dx[idx    ] = cb_r[idx];
			dx[idx+n*3] = cb_v[idx];
		}
		double *res = new double[2*n*3];
		for (int idx = 0; idx < 2*n*3; idx++) {
			res[idx] = 0;
		}
		f(dx, res);
		for (int idx = i * 2*3*n; idx < (i+1)*2*n*3; idx++) {
			adams_dx[idx] = res[idx - i * 2*3*n];
		}
		t+=h;
		delete[] dx;
		delete[] res;
	}
	size = 8;
	adams_a = new double[size]{
		434241./120960, -1152169./120960, 2183877./120960, -2664477./120960, 2102243./120960, -1041723./120960, 295767./120960, -36799./120960
	};
	int i = 7;

/*
----------------------------------------------
	MAIN LOOP
----------------------------------------------*/
	// int i = 0;
	while (t < T) {
		// if (i % 500 == 1) {
			std::cout << t << " ";
			for (int i = 0; i < n * 3; i++) {
				std::cout << cb_r[i] << " ";
			}
			std::cout << "\n";
		// }
		// rk(h);
		adams(h);

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

/*
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
	delete[] k;*/
