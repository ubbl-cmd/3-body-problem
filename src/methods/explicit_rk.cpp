#include "explicit_rk.hpp"

void explicit_rk::rk_step(double h, int size, double *x, std::function<double*(int, double*, void*)> f, void * data) {
	double *k     = new double[size * a_size];
	double *local = new double[size];

	double *dx;

	for (int ki = 0; ki < a_size; ki++) {
		for (int i = 0; i < size; i++) {
			local[i] = x[i];
			for (int kis = 0; kis < ki; kis++) {
				local[i] += h * a[ki * a_size + kis] * k[kis * size + i];
			}
		}
		dx = f(size, local, data);
		for (int i = 0; i < size; i++) {
			k[ki * size + i] = dx[i];
		}
		delete[] dx;
	}
	for (int ki = 0; ki < a_size; ki++) {
		for (int i = 0; i < size; i++) {
			x[i] += h * b[ki] * k[ki * size + i];
		}
	}
	delete[] local;
	delete[] k;
}
