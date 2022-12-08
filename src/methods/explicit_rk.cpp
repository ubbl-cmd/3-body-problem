#include "explicit_rk.hpp"

void explicit_rk::rk_step(long double h, int size, long double *x, std::function<long double*(int, long double*, void*)> f, void * data) {
	long double *k     = new long double[size * a_size];
	long double *local = new long double[size];

	long double *dx;

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
