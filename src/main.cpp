#include "methods/explicit_rk.hpp"
#include "methods/explicit_adams.hpp"

#include <iostream>
#include <fstream>
#include <cmath>

long double pi = 2 * acos(0.0);
struct keplerian {
	long double a;
	long double e;
	long double w;
	long double om;
	long double incl;
	long double M;
};

struct keplerian kepl;

const int DO_EVERY_N_STEPS = 500;

std::ofstream ofs;

long double* vmv(long double* a, long double* b) {
	long double* res = new long double[3];
	res[0] = a[1] * b[2] - a[2] * b[1];
	res[1] = a[2] * b[0] - a[0] * b[2];
	res[2] = a[0] * b[1] - a[1] * b[0];
	return res;
}


long double veclen(long double * v) {
	long double res = 0;
	for (int i = 0; i < 3; i++) {
		res += v[i] * v[i];
	}
	return sqrt(res);
}


long double * f(int n, long double *x, void *data) {
	const long double G = 6.67430e-20;
	long double *mass = (long double *) data;
	long double* res = new long double[n];

	long double* r = new long double[3];
	for (int i = 0; i < n / 2; i++) {
		res[i] = x[i + n / 2];
		res[i + n / 2] = 0;
	}
	for (int i = 0; i < n / 6; i++) {
		for (int j = 0; j < n / 6; j++) {
			if (j == i) {
				continue;
			}
			for (int idx = 0; idx < 3; idx++) {
				r[idx] = x[j * 3 + idx] - x[i * 3 + idx];
			}
			for (int idx = 0; idx < 3; idx++) {
				res[3 * i + idx + n / 2] += r[idx] * (G * mass[j] / pow(veclen(r), 3));
			}
		}
	}
	delete[] r;
	return res;
}

long double * f_analit(int n, long double *x, void *data) {
	const long double G = 6.67430e-20;
	long double *mass = (long double *) data;
	long double* res = new long double[n];

	long double* r = new long double[3];
	for (int i = 0; i < n / 2; i++) {
		res[i] = x[i + n / 2];
		res[i + n / 2] = 0;
	}
	for (int i = 0; i < n / 6 / 2; i++) {
		for (int j = 0; j < n / 6; j++) {
			if (j == i) {
				continue;
			}
			for (int idx = 0; idx < 3; idx++) {
				r[idx] = x[j * 3 + idx] - x[i * 3 + idx];
			}
			for (int idx = 0; idx < 3; idx++) {
				res[3 * i + idx + n / 2] += r[idx] * (G * mass[j] / pow(veclen(r), 3));
			}
		}
	}
	for (int idx = 0; idx < 3; idx++) {
		res[3 + idx] = 0;
		res[6 + 3 + idx] = 0;
	}
	delete[] r;
	return res;
}


void open_ofs(std::string filename) {
	if (!ofs.is_open()) {
		ofs.open(filename, std::ios::binary);
		if (!ofs.is_open()) {
			std::cout << "Error opening file \"" << filename << "\"!" << std::endl;
		}
	}
}


void print_trajectory(long double t, int size, long double *x, long double *mass) {
	open_ofs("trajectory");
	ofs << t << " ";
	for (int i = 0; i < size / 2; i++) {
		ofs << x[i] << " ";
	}
	ofs << std::endl;
}


long double *calculate_barycenter_coorinates(int size, long double *x, long double *mass) {
	long double *barycenter_coorinates = new long double[3] {0,0,0};
	long double mass_sum = 0;
	for (int i = 0; i < size / 6; i++) {
		mass_sum += mass[i];
	}
	for (int i = 0; i < size / 2; i++) {
		barycenter_coorinates[i % 3] += x[i] * mass[i/3] / mass_sum;
	}
	return barycenter_coorinates;
}

void print_barycenter_coorinates_and_trajectory(long double t, int size, long double *x, long double *mass) {
	open_ofs("barycenter_coorinates");
	ofs << t << " ";
	for (int i = 0; i < size / 2; i++) {
		ofs << x[i] << " ";
	}
	long double * barycenter_coorinates = calculate_barycenter_coorinates(size, x, mass);
	for (int j = 0; j < 3; j++) {
		ofs << barycenter_coorinates[j] << " ";
	}
	delete[] barycenter_coorinates;
	ofs << std::endl;
}

long double *calculate_barycenter_speed3(int size, long double *x, long double *mass) {
	long double *barycenter_speed = new long double[3] {0,0,0};
	long double mass_sum = 0;
	for (int i = 0; i < size / 6; i++) {
		mass_sum += mass[i];
	}
	for (int i = 0; i < size / 2; i++) {
		barycenter_speed[i % 3] += x[i + size / 2] * mass[i/3] / mass_sum;
	}
	return barycenter_speed;
}

long double calculate_barycenter_speed(int size, long double *x, long double *mass) {
	return veclen(calculate_barycenter_speed3(size, x, mass));
}

void print_barycenter_speed(long double t, int size, long double *x, long double *mass) {
	open_ofs("barycenter_speed");
	ofs << t << " ";
	long double barycenter_speed = calculate_barycenter_speed(size, x, mass);
	ofs << barycenter_speed << " ";
	ofs << std::endl;
}


long double calculate_kinetic_energy(int size, long double *x, long double *mass) {
	long double kinetic_energy = 0;
	long double * r = new long double[3];
	for (int i = 0; i < size / 6; i++) {
		r[0] = x[size / 2 + i * 3 + 0] * x[size / 2 + i * 3 + 0];
		r[1] = x[size / 2 + i * 3 + 1] * x[size / 2 + i * 3 + 1];
		r[2] = x[size / 2 + i * 3 + 2] * x[size / 2 + i * 3 + 2];
		long double rlen = (r[0] + r[1] + r[2]) / 2;
		kinetic_energy += mass[i] * rlen;
	}
	delete[] r;
	return kinetic_energy;
}

void print_kinetic_energy(long double t, int size, long double *x, long double *mass) {
	open_ofs("kinetic_energy");
	ofs << t << " ";
	long double kinetic_energy = calculate_kinetic_energy(size, x, mass);
	ofs << kinetic_energy << " ";
	ofs << std::endl;
}

long double calculate_potential_energy(int size, long double *x, long double *mass) {
	const long double G = 6.67430e-20;
	long double potential_energy = 0;
	long double * r = new long double[3];
	for (int i = 0; i < size / 6; i++) {
		for (int j = i+1; j < size / 6; j++) {
			if (i == j) continue;
			r[0] = x[i*3+0] - x[j*3+0];
			r[1] = x[i*3+1] - x[j*3+1];
			r[2] = x[i*3+2] - x[j*3+2];
			long double Rlen = powl(r[0] * r[0] + r[1] * r[1] + r[2] * r[2], 0.5);
			potential_energy -= G * mass[i] * mass[j] / Rlen;
		}
	}
	delete[] r;
	return potential_energy;
}

void print_potential_energy(long double t, int size, long double *x, long double *mass) {
	open_ofs("potential_energy");
	ofs << t << " ";
	long double potential_energy = calculate_potential_energy(size, x, mass);
	ofs << potential_energy << " ";
	ofs << std::endl;
}
void print_energy(long double t, int size, long double *x, long double *mass) {
	open_ofs("energy");
	ofs << t << " ";
	long double kinetic_energy = calculate_kinetic_energy(size, x, mass);
	long double potential_energy = calculate_potential_energy(size, x, mass);
	ofs << kinetic_energy << " ";
	ofs << potential_energy << " ";
	ofs << kinetic_energy + potential_energy << " ";
	ofs << std::endl;
}


long double *calculate_angular_momentum(int size, long double *x, long double *mass) {
	long double *r = new long double[3];
	long double *dist = new long double[3];
	long double *angular_momentum = new long double[3] {0,0,0};

	for (int i = 0; i < size / 6; i++) {
		for (int j = 0; j < size / 6; j++) {
			if (i == j) continue;
			r[0] = x[size/2 + i*3+0] * mass[i];
			r[1] = x[size/2 + i*3+1] * mass[i];
			r[2] = x[size/2 + i*3+2] * mass[i];
			dist[0] = -x[j*3+0] + x[i*3+0];
			dist[1] = -x[j*3+1] + x[i*3+1];
			dist[2] = -x[j*3+2] + x[i*3+2];
			long double * tmp = vmv(r, dist);
			for (int k = 0; k < 3; k++) {
				angular_momentum[k] += tmp[k];
			}
			delete[] tmp;
		}
	}
	delete[] r;
	delete[] dist;
	return angular_momentum;
}

void print_angular_momentum(long double t, int size, long double *x, long double *mass) {
	open_ofs("angular_momentum");
	ofs << t << " ";
	long double *angular_momentum = calculate_angular_momentum(size, x, mass);
	ofs << veclen(angular_momentum) << " ";
	for (int i = 0; i < 3; i++) {
		ofs << angular_momentum[i] << " ";
	}
	ofs << std::endl;
}



keplerian cartesian_to_kepler(long double *x, long double mu) {
	long double *r = new long double[3]{x[0],x[1],x[2]};
	long double *rdt = new long double[3]{x[3],x[4],x[5]};

	long double *h = vmv(r,rdt);

	long double *e = vmv(rdt, h);
	for (int i = 0; i < 3; i++) {
		e[i] /= mu;
		e[i] -= r[i]/veclen(r);
	}

	long double *n = new long double[3]{-1*h[1],h[0],0};

	long double v;
	v = acos(
		(e[0]*r[0]+e[1]*r[1]+e[2]*r[2]) /
		veclen(e) / veclen(r)
	);
	if (r[0]*rdt[0]+r[1]*rdt[1]+r[2]*rdt[2] < 0) {
		v = 2 * pi - v;
	}

	long double inclination = acos(
		h[2] / veclen(h)
	);

	long double E = 2 * atan(
		tanl(v / 2) /
		sqrt(
			(1+veclen(e)) / (1 - veclen(e))
		)
	);

	long double sigma = acos(
		n[0] / veclen(n)
	);
	if (n[1] < 0) {
		sigma = 2 * pi - sigma;
	}

	long double omega = acos(
		(n[0]*e[0]+n[1]*e[1]+n[2]*e[2]) / veclen(e) / veclen(n)
	);

	if (e[2] < 0) {
		omega = 2 * pi - omega;
	}

	long double M = E - veclen(e) * sinl(E);

	long double a = 1. / ((2./veclen(r)) - (veclen(rdt) * veclen(rdt)) / mu);
	struct keplerian kepl;
	kepl.a = a;
	kepl.e = veclen(e);
	kepl.w = omega;
	kepl.om = sigma;
    kepl.incl = inclination;
	kepl.M = M;
    delete[] e;
    delete[] n;
    delete[] h;
	return kepl;
}

long double arctan2(long double y, long double x) {
	if (x > 0) {
		return atan(y/x);
	} else if (y >= 0 && x < 0) {
		return atan(y/x) + pi;
	} else if (y < 0 && x < 0) {
		return atan(y/x) - pi;
	} else if (y > 0 && x == 0) {
		return pi / 2;
	} else if (y < 0 && x == 0) {
		return -1 * pi / 2;
	}
	return 0;
}
long double t_counter = 0;
long double E_t = 0;
long double* kepler_to_cartesian(
	long double a,
	long double e,
	long double w,
	long double om,
	long double incl,
	long double M,
	long double t,
	long double mu) {
		long double M_t = M + t * sqrt(mu / a / a / a);
		while (M_t < 0) {
			M_t += 2 * pi;
		}
		while (M_t >= 2 * pi) {
			M_t -= 2 * pi;
		}
		if (t == 0)  {
			E_t = M_t;
		} else {
			E_t += (t - t_counter) * sqrt(mu / a / a / a);
			while (t_counter <= t) {
				E_t = E_t - (E_t - e * sinl(E_t) - M_t) / (1 - e * cosl(E_t));
				t_counter += 1;
			}
		}

		long double v_t = 2 * arctan2(
			sqrt(1+e)*sinl(E_t/2), sqrt(1-e)*cosl(E_t/2)
		);

		long double rc_t = a * (1 - e * cosl(E_t));

		long double *o_t = new long double[3]{ cosl(v_t), sinl(v_t), 0.0 };
		for (int k = 0; k < 3; k++) {
			o_t[k] *= rc_t;
		}

		long double* x = new long double[3] {
			o_t[0] * (cosl(w)*cosl(om)-sinl(w)*cosl(incl)*sinl(om)) - o_t[1] * (sinl(w)*cosl(om) + cosl(w)*cosl(incl)*sinl(om)),
			o_t[0] * (cosl(w)*sinl(om)+sinl(w)*cosl(incl)*cosl(om)) + o_t[1] * (cosl(w)*cosl(incl)*cosl(om) - sinl(w)*sinl(om)),
			o_t[0] * (sinl(w)*sin(incl))                            + o_t[1] * (cosl(w)*sinl(incl))
		};
		delete[] o_t;
		return x;
}

void print_analitical_solution(long double t, int size, long double *x, long double *mass) {
	long double *r = kepler_to_cartesian(
		kepl.a,
		kepl.e,
		kepl.w,
		kepl.om,
		kepl.incl,
		kepl.M,
		t,
		1.98847e30 * 6.67430e-20);
	open_ofs("analit");
	ofs << t << " ";
	long double *diff = new long double[3] {
		x[0]-r[0],
		x[1]-r[1],
		x[2]-r[2],
	};
	ofs << veclen(diff)/veclen(r);
	delete[] diff;
	// for (int i = 0; i < 3; i++) {
	// 	ofs << x[i] << " ";
	// }
	// for (int i = 0; i < 3; i++) {
	// 	ofs << r[i]<< " ";
	// }
	ofs << std::endl;
	delete[] r;
}


void process_rk4(int size, long double *x, long double *mass, long double h, long double T, std::function<long double*(int, long double*, void*)> func, std::function<void(long double, int, long double*, long double*)> every_step_function) {
	int a_size = 4;
	long double *a_rk4 = new long double[a_size * a_size] {
		0.,  0,   0, 0,
		0.5, 0,   0, 0,
		0,   0.5, 0, 0,
		0,   0,   1, 0
	};

	long double *b_rk4 = new long double[a_size] {
		1./6, 1./3, 1./3, 1./6
	};

	explicit_rk rk4(a_size, a_rk4, b_rk4);


	long double t = 0;
	int step = 0;
	while (t <= T) {
		if (step % DO_EVERY_N_STEPS == 1) {
			every_step_function(t, size, x, mass);
		}
		rk4.rk_step(h, size, x, func, (void *) mass);


		step++;
		t += h;
	}
	delete[] a_rk4;
	delete[] b_rk4;
}


void process_dorpri8(int size, long double *x, long double *mass, long double h, long double T, std::function<long double*(int, long double*, void*)> func, std::function<void(long double, int, long double*, long double*)> every_step_function) {
	int dorpi8_size = 13;
	long double *a_dorpri8 = new long double[dorpi8_size * dorpi8_size] {
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,

		1./18, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,

		1./48, 1./16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,

		1./32, 0, 3./32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,

		5./16, 0, -75./64, 75./64, 0, 0, 0, 0, 0, 0, 0, 0, 0,

		3./80, 0, 0, 3./16, 3./20, 0, 0, 0, 0, 0, 0, 0, 0,

		29443841./614563906, 0, 0, 77736538./6922538347, -28693883./11125000000, 23124283./1800000000, 0, 0, 0, 0, 0, 0, 0,

		16016141./946692911, 0, 0, 61564180./158732637, 22789713./633445777, 545815736./2771057229, -180193667./1043307555, 0, 0, 0, 0, 0, 0,

		39632708./573591083, 0, 0, -433636366./683701615, -421739975./2616292301, 100302831./723423059, 790204164./839813087, 800635310./3783071287, 0, 0, 0, 0, 0,

		246121993./1340847787, 0, 0, -37695042795./15268766246, -309121744./1061227803, -12992083./490766935, 6005943433./2108947869, 393006217./1396673457, 123782331./1001029789, 0, 0, 0, 0,

		-1028468189./846180014, 0, 0, 8478235783./506512852, 1311729495./1432422823, -10304129995./1701304382, -48777925059./3047939560, 15336726248./1032824649, -45442868181./3398467696, 3065993473./597172653, 0, 0, 0,

		185892177./718116043, 0, 0, -3185094517./667107341, -477755414./1098053517, -703635378./230739211, 5731566787./1027545527, 5232866602./850066563, -4093664535./808688257, 3962137247./1805957418, 65686385./487910083, 0, 0,
		14005451. / 335480064, 0, 0, 0, 0, -59238493. / 1068277825, 181606767. / 758867731, 561292985. / 797845732, -1041891430. / 1371343529, 760417239. / 1151165299, 118820643. / 751138087, -528747749. / 2220607170, 1. / 4,
	};

	long double *b_dorpri8 = new long double[dorpi8_size] {
		14005451. / 335480064, 0, 0, 0, 0, -59238493. / 1068277825, 181606767. / 758867731, 561292985. / 797845732, -1041891430. / 1371343529, 760417239. / 1151165299, 118820643. / 751138087, -528747749. / 2220607170, 1. / 4
	};

	explicit_rk *dorpri8 = new explicit_rk(dorpi8_size, a_dorpri8, b_dorpri8);


	long double t = 0;
	int step = 0;
	while (t <= T) {
		if (step % DO_EVERY_N_STEPS == 1) {
			every_step_function(t, size, x, mass);
		}
		dorpri8->rk_step(h, size, x, func, (void *) mass);


		step++;
		t += h;
	}
	delete dorpri8;
}

void process_adams8(int size, long double *x, long double *mass, long double h, long double T, std::function<long double*(int, long double*, void*)> func, std::function<void(long double, int, long double*, long double*)> every_step_function) {
	int dorpi8_size = 13;
	long double *a_dorpri8 = new long double[dorpi8_size * dorpi8_size] {
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,

		1./18, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,

		1./48, 1./16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,

		1./32, 0, 3./32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,

		5./16, 0, -75./64, 75./64, 0, 0, 0, 0, 0, 0, 0, 0, 0,

		3./80, 0, 0, 3./16, 3./20, 0, 0, 0, 0, 0, 0, 0, 0,

		29443841./614563906, 0, 0, 77736538./6922538347, -28693883./11125000000, 23124283./1800000000, 0, 0, 0, 0, 0, 0, 0,

		16016141./946692911, 0, 0, 61564180./158732637, 22789713./633445777, 545815736./2771057229, -180193667./1043307555, 0, 0, 0, 0, 0, 0,

		39632708./573591083, 0, 0, -433636366./683701615, -421739975./2616292301, 100302831./723423059, 790204164./839813087, 800635310./3783071287, 0, 0, 0, 0, 0,

		246121993./1340847787, 0, 0, -37695042795./15268766246, -309121744./1061227803, -12992083./490766935, 6005943433./2108947869, 393006217./1396673457, 123782331./1001029789, 0, 0, 0, 0,

		-1028468189./846180014, 0, 0, 8478235783./506512852, 1311729495./1432422823, -10304129995./1701304382, -48777925059./3047939560, 15336726248./1032824649, -45442868181./3398467696, 3065993473./597172653, 0, 0, 0,

		185892177./718116043, 0, 0, -3185094517./667107341, -477755414./1098053517, -703635378./230739211, 5731566787./1027545527, 5232866602./850066563, -4093664535./808688257, 3962137247./1805957418, 65686385./487910083, 0, 0,
		14005451. / 335480064, 0, 0, 0, 0, -59238493. / 1068277825, 181606767. / 758867731, 561292985. / 797845732, -1041891430. / 1371343529, 760417239. / 1151165299, 118820643. / 751138087, -528747749. / 2220607170, 1. / 4,
	};

	long double *b_dorpri8 = new long double[dorpi8_size] {
		14005451. / 335480064, 0, 0, 0, 0, -59238493. / 1068277825, 181606767. / 758867731, 561292985. / 797845732, -1041891430. / 1371343529, 760417239. / 1151165299, 118820643. / 751138087, -528747749. / 2220607170, 1. / 4
	};

	explicit_rk *dorpri8 = new explicit_rk(dorpi8_size, a_dorpri8, b_dorpri8);

	int a_size = 8;
	long double *adams_a = new long double[a_size] {
		434241./120960, -1152169./120960, 2183877./120960, -2664477./120960, 2102243./120960, -1041723./120960, 295767./120960, -36799./120960
	};
	explicit_adams *adams8 = new explicit_adams(a_size, adams_a);

	long double t = h*7;
	int step = 0;
	every_step_function(t, size, x, mass);
	adams8->razgonka(h, dorpri8, size, x, func, mass);
	while (t <= T) {
		if (step % DO_EVERY_N_STEPS == 1) {
			every_step_function(t, size, x, mass);
		}
		adams8->adams_step(h, size, x, func, (void *) mass);


		step++;
		t += h;
	}
}

int main() {
	ofs.setf(std::ios::fixed);  // вывод в фиксированном формате
	ofs.precision(50);      // вывод до 6 знака после точки, включительно
	int n;
	long double *x;
	long double *mass;

	n=3;
	x = new long double[6*n] {
//		R[km/s]
		 1.453004563709436E+08, 2.978607390059390E+07,  2.866143171530403E+04,
		1.455259266863196E+08, 2.949542371875033E+07, -4.639949931440875E+03,
		-1.359757021307868E+06, 1.333076360298289E+05,  3.057372007160987E+04,
//		dR/dt[km/s^2]
		-6.397102954570147E+00,  2.906192414721306E+01, -1.203250030542335E-03,
		-5.557389967864156E+00, 2.971153423303235E+01, -1.096802103212369E-02,
		-2.335983427600637E-04, -1.571626050150792E-02,  1.314580230199784E-04
	};
	mass = new long double[n] {
//		[kg]
		5.9722e24,
		7.349e22,
		1.98847e30,
	};
	long double h = 10.0;
	long double T = 31536000.;




	/* TRAJECTORY */

	process_rk4(n * 6, x, mass, h, T, f, &print_trajectory);
	// process_dorpri8(n * 6, x, mass, h, T, f, &print_trajectory);
	// process_adams8(n * 6, x, mass, h, T, f, &print_trajectory);


	/* BARYCENTER COORDINATES AND TRAJECTORY */

	// process_rk4(n * 6, x, mass, h, T, f, &print_barycenter_coorinates_and_trajectory);
	// process_dorpri8(n * 6, x, mass, h, T, f, &print_barycenter_coorinates_and_trajectory);
	// process_adams8(n * 6, x, mass, h, T, f, &print_barycenter_coorinates_and_trajectory);


	/* BARYCENTER SPEED */

	// process_rk4(n * 6, x, mass, h, T, f, &print_barycenter_speed);
	// process_dorpri8(n * 6, x, mass, h, T, f, &print_barycenter_speed);
	// process_adams8(n * 6, x, mass, h, T, f, &print_barycenter_speed);


	/* KINETIC ENERGY */

	// process_rk4(n * 6, x, mass, h, T, f, &print_kinetic_energy);
	// process_dorpri8(n * 6, x, mass, h, T, f, &print_kinetic_energy);
	// process_adams8(n * 6, x, mass, h, T, f, &print_kinetic_energy);


	/* POTENTIAL ENERGY */

	// process_rk4(n * 6, x, mass, h, T, f, &print_potential_energy);
	// process_dorpri8(n * 6, x, mass, h, T, f, &print_potential_energy);
	// process_adams8(n * 6, x, mass, h, T, f, &print_potential_energy);


	/* ENERGY */

	// process_rk4(n * 6, x, mass, h, T, f, &print_energy);
	// process_dorpri8(n * 6, x, mass, h, T, f, &print_energy);
	// process_adams8(n * 6, x, mass, h, T, f, &print_energy);

	/* ANGULAR MOMENTUM */

	 // process_rk4(n * 6, x, mass, h, T, f, &print_angular_momentum);
	// process_dorpri8(n * 6, x, mass, h, T, f, &print_energy);
	// process_adams8(n * 6, x, mass, h, T, f, &print_energy);

	/*ANALITICAL*/

	// n = 2;
	// x = new long double[6*n] {
	// 	2.005683434185720E+07,  1.337069254521444E+08,  5.799434294798527E+07,
	// 	0,0,0,
	// 	-2.994817357192201E+01,  3.864228014740471E+00,  1.675675038999799E+00,
	// 	0,0,0,
	// };
	// mass = new long double[n] {
	// 	5.9722e24,
	// 	1.98847e30,
	// };
	// long double *bars = calculate_barycenter_speed3(n*6,x,mass);
	// for (int i = 0; i < n*3; i++){
	// 	x[n*3 + i] -= bars[i%3];
	// }
	// long double *analit_x = new long double[6]{
	// 	x[0],x[1],x[2],
	// 	x[6],x[7],x[8]
	// };
	// kepl = cartesian_to_kepler(analit_x, 1.98847e30 * 6.67430e-20);
	// process_rk4(n * 6, x, mass, h, T, f_analit, &print_analitical_solution);


	if (ofs.is_open()) {
		ofs.close();
	}
	delete[] x;
	delete[] mass;
	return 0;
}
