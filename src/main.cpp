#include "methods/explicit_rk.hpp"
#include "methods/explicit_adams.hpp"

#include <iostream>
#include <fstream>
#include <cmath>

double veclen(double * v) {
	double res = 0;
	for (int i = 0; i < 3; i++) {
		res += v[i] * v[i];
	}
	return sqrt(res);
}

double * f(int n, double *x, void *data) {
	const double G = 6.67430e-20;
	double *mass = (double *) data;
	double* res = new double[n];

	double* r = new double[3];
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

const int DO_EVERY_N_STEPS = 500;

std::ofstream ofs;

void open_ofs(std::string filename) {
	if (!ofs.is_open()) {
		ofs.open(filename, std::ios::binary);
		if (!ofs.is_open()) {
			std::cout << "Error opening file \"" << filename << "\"!" << std::endl;
		}
	}
}

void print_trajectory(double t, int size, double *x, double *mass) {
	open_ofs("trajectory");
	ofs << t << " ";
	for (int i = 0; i < size / 2; i++) {
		ofs << x[i] << " ";
	}
	ofs << std::endl;
}

double *calculate_barycenter_coorinates(int size, double *x, double *mass) {
	double *barycenter_coorinates = new double[3] {0,0,0};
	double mass_sum = 0;
	for (int i = 0; i < size / 6; i++) {
		mass_sum += mass[i];
	}
	for (int i = 0; i < size / 2; i++) {
		barycenter_coorinates[i % 3] += x[i] * mass[i/3] / mass_sum;
	}
	return barycenter_coorinates;
}

void print_barycenter_coorinates_and_trajectory(double t, int size, double *x, double *mass) {
	open_ofs("barycenter_coorinates");
	ofs << t << " ";
	for (int i = 0; i < size / 2; i++) {
		ofs << x[i] << " ";
	}
	double * barycenter_coorinates = calculate_barycenter_coorinates(size, x, mass);
	for (int j = 0; j < 3; j++) {
		ofs << barycenter_coorinates[j] << " ";
	}
	delete[] barycenter_coorinates;
	ofs << std::endl;
}

double calculate_barycenter_speed(int size, double *x, double *mass) {
	double *barycenter_speed = new double[3] {0,0,0};
	double mass_sum = 0;
	for (int i = 0; i < size / 6; i++) {
		mass_sum += mass[i];
	}
	for (int i = 0; i < size / 2; i++) {
		barycenter_speed[i % 3] += x[i + size / 2] * mass[i/3] / mass_sum;
	}
	return veclen(barycenter_speed);
}

double *calculate_barycenter_speed3(int size, double *x, double *mass) {
	double *barycenter_speed = new double[3] {0,0,0};
	double mass_sum = 0;
	for (int i = 0; i < size / 6; i++) {
		mass_sum += mass[i];
	}
	for (int i = 0; i < size / 2; i++) {
		barycenter_speed[i % 3] += x[i + size / 2] * mass[i/3] / mass_sum;
	}
	return barycenter_speed;
}

void print_barycenter_speed(double t, int size, double *x, double *mass) {
	open_ofs("barycenter_speed");
	ofs << t << " ";
	double barycenter_speed = calculate_barycenter_speed(size, x, mass);
	ofs << barycenter_speed << " ";
	ofs << std::endl;
}

double calculate_kinetic_energy(int size, double *x, double *mass) {
	double kinetic_energy = 0;
	double * r = new double[3];
	double * bars = calculate_barycenter_speed3(size, x, mass);
	for (int i = 0; i < size / 6; i++) {
		r[0] = x[size / 2 + i * 3 + 0] + bars[0];
		r[1] = x[size / 2 + i * 3 + 1] + bars[1];
		r[2] = x[size / 2 + i * 3 + 2] + bars[2];
		double rlen = veclen(r);
		kinetic_energy += mass[i] * rlen * rlen / 2;
	}
	delete[] r;
	return kinetic_energy;
}

void print_kinetic_energy(double t, int size, double *x, double *mass) {
	open_ofs("kinetic_energy");
	ofs << t << " ";
	double kinetic_energy = calculate_kinetic_energy(size, x, mass);
	ofs << kinetic_energy << " ";
	ofs << std::endl;
}

double calculate_potential_energy(int size, double *x, double *mass) {
	const double G = 6.67430e-20;
	double potential_energy = 0;
	double * r = new double[3];
	double * bar = calculate_barycenter_coorinates(size, x, mass);
	double mass_sum = 0;
	for (int i = 0; i < size / 6; i++) {
		mass_sum += mass[i];
	}
	for (int i = 0; i < size / 6; i++) {
		r[0] = (x[i * 3 + 0] - bar[0]);
		r[1] = (x[i * 3 + 1] - bar[1]);
		r[2] = (x[i * 3 + 2] - bar[2]);
		double rlen = veclen(r);
		potential_energy -= G * (mass_sum - mass[i]) * mass[i] / rlen;
	}
	delete[] r;
	delete[] bar;
	return potential_energy;
}

void print_potential_energy(double t, int size, double *x, double *mass) {
	open_ofs("potential_energy");
	ofs << t << " ";
	double potential_energy = calculate_potential_energy(size, x, mass);
	ofs << potential_energy << " ";
	ofs << std::endl;
}

void print_energy(double t, int size, double *x, double *mass) {
	open_ofs("energy");
	ofs << t << " ";
	double kinetic_energy = calculate_kinetic_energy(size, x, mass);
	double potential_energy = calculate_potential_energy(size, x, mass);
	ofs << kinetic_energy << " ";
	ofs << potential_energy << " ";
	ofs << std::endl;
}

void process_rk4(int size, double *x, double *mass, double h, double T, std::function<void(double, int, double*, double*)> every_step_function) {
	int a_size = 4;
	double *a_rk4 = new double[a_size * a_size] {
		0.,  0,   0, 0,
		0.5, 0,   0, 0,
		0,   0.5, 0, 0,
		0,   0,   1, 0
	};

	double *b_rk4 = new double[a_size] {
		1./6, 1./3, 1./3, 1./6
	};

	explicit_rk rk4(a_size, a_rk4, b_rk4);


	double t = 0;
	int step = 0;
	while (t <= T) {
		if (step % DO_EVERY_N_STEPS == 0) {
			every_step_function(t, size, x, mass);
		}
		rk4.rk_step(h, size, x, &f, (void *) mass);


		step++;
		t += h;
	}
	delete[] a_rk4;
	delete[] b_rk4;
}


void process_dorpri8(int size, double *x, double *mass, double h, double T, std::function<void(double, int, double*, double*)> every_step_function) {
	int dorpi8_size = 13;
	double *a_dorpri8 = new double[dorpi8_size * dorpi8_size] {
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
		403863854./491063109, 0, 0, -5068492393./434740067, -411421997./543043805, 652783627./914296604, 11173962825./925320556, -13158990841./6184727034, 3936647629./1978049680, -160528059./685178525, 248638103./1413531060, 0, 0
	};

	double *b_dorpri8 = new double[dorpi8_size] {
		14005451./335480064, 0, 0, 0, 0, -59238493./1068277825, 181606676./758867731, 561292985./797845732, -1041891430./1371343529, 760417239./1151165299, 118820643./751138087, -5228747749./2220607170, 1./4
	};

	explicit_rk *dorpri8 = new explicit_rk(dorpi8_size, a_dorpri8, b_dorpri8);


	double t = 0;
	int step = 0;
	while (t <= T) {
		if (step % DO_EVERY_N_STEPS == 0) {
			every_step_function(t, size, x, mass);
		}
		dorpri8->rk_step(h, size, x, &f, (void *) mass);


		step++;
		t += h;
	}
	delete dorpri8;
}

void process_adams8(int size, double *x, double *mass, double h, double T, std::function<void(double, int, double*, double*)> every_step_function) {
	int dorpi8_size = 13;
	double *a_dorpri8 = new double[dorpi8_size * dorpi8_size] {
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
		403863854./491063109, 0, 0, -5068492393./434740067, -411421997./543043805, 652783627./914296604, 11173962825./925320556, -13158990841./6184727034, 3936647629./1978049680, -160528059./685178525, 248638103./1413531060, 0, 0
	};

	double *b_dorpri8 = new double[dorpi8_size] {
		14005451./335480064, 0, 0, 0, 0, -59238493./1068277825, 181606676./758867731, 561292985./797845732, -1041891430./1371343529, 760417239./1151165299, 118820643./751138087, -5228747749./2220607170, 1./4
	};

	explicit_rk *dorpri8 = new explicit_rk(dorpi8_size, a_dorpri8, b_dorpri8);

	int a_size = 8;
	double *adams_a = new double[a_size] {
		434241./120960, -1152169./120960, 2183877./120960, -2664477./120960, 2102243./120960, -1041723./120960, 295767./120960, -36799./120960
	};
	explicit_adams *adams8 = new explicit_adams(a_size, adams_a);

	double t = h*7;
	int step = 0;
	every_step_function(t, size, x, mass);
	adams8->razgonka(h, dorpri8, size, x, f, mass);
	while (t <= T) {
		if (step % DO_EVERY_N_STEPS == 0) {
			every_step_function(t, size, x, mass);
		}
		adams8->adams_step(h, size, x, &f, (void *) mass);


		step++;
		t += h;
	}
}

int main() {
	int n = 3;
	/*
	 R[earth], R[moon], R[sun],
	 V[earth], V[moon], V[sun]
	 */
	double *x = new double[6*n] {

//		=========RADIUS VECTORS=========
		 1.453004563709436E+08, 2.978607390059390E+07,  2.866143171530403E+04,
		 1.455259266863196E+08, 2.949542371875033E+07, -4.639949931440875E+03,
		-1.359757021307868E+06, 1.333076360298289E+05,  3.057372007160987E+04,
//		=========SPEED VECTORS=========
		-6.397102954570147E+00,  2.906192414721306E+01, -1.203250030542335E-03,
		-5.557389967864156E+00, 2.971153423303235E+01, -1.096802103212369E-02,
		-2.335983427600637E-04, -1.571626050150792E-02,  1.314580230199784E-04
	};
	/*
	 M[earth], M[moon], M[sun]
	 */
	double *mass = new double[n] {
		5.9722e24,
		7.349e22,
		1.98847e30
	};
	double h = 100.0;
	double T = 31536000.;

	// double *bars = calculate_barycenter_speed3(n*6,x,mass);
	// for (int i = 0; i < n*3; i++){
	// 	x[n*3 + i] -= bars[i%3];
	// }





	/* TRAJECTORY */

	// process_rk4(n * 6, x, mass, h, T, &print_trajectory);
	// process_dorpri8(n * 6, x, mass, h, T, &print_trajectory);
	// process_adams8(n * 6, x, mass, h, T, &print_trajectory);


	/* BARYCENTER COORDINATES AND TRAJECTORY */

	// process_rk4(n * 6, x, mass, h, T, &print_barycenter_coorinates_and_trajectory);
	// process_dorpri8(n * 6, x, mass, h, T, &print_barycenter_coorinates_and_trajectory);
	// process_adams8(n * 6, x, mass, h, T, &print_barycenter_coorinates_and_trajectory);


	/* BARYCENTER SPEED */

	// process_rk4(n * 6, x, mass, h, T, &print_barycenter_speed);
	// process_dorpri8(n * 6, x, mass, h, T, &print_barycenter_speed);
	// process_adams8(n * 6, x, mass, h, T, &print_barycenter_speed);


	/* KINETIC ENERGY */

	// process_rk4(n * 6, x, mass, h, T, &print_kinetic_energy);
	// process_dorpri8(n * 6, x, mass, h, T, &print_kinetic_energy);
	// process_adams8(n * 6, x, mass, h, T, &print_kinetic_energy);


	/* POTENTIAL ENERGY */

	// process_rk4(n * 6, x, mass, h, T, &print_potential_energy);
	// process_dorpri8(n * 6, x, mass, h, T, &print_potential_energy);
	// process_adams8(n * 6, x, mass, h, T, &print_potential_energy);


	/* ENERGY */

	process_rk4(n * 6, x, mass, h, T, &print_energy);
	// process_dorpri8(n * 6, x, mass, h, T, &print_energy);
	// process_adams8(n * 6, x, mass, h, T, &print_energy);


	if (ofs.is_open()) {
		ofs.close();
	}
	delete[] x;
	delete[] mass;
	return 0;
}
