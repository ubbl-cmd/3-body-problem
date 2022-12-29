#include "explicit_adams.hpp"

void explicit_adams::razgonka(long double h, explicit_rk *rk, int size, long double *x, std::function<long double*(int, long double*, void*)> f, void *data) {
	adams_dx = new long double[a_size * size];
	long double *local = new long double[size];
	long double *dx;
	for (int i = 0; i < a_size; i++) {
		for (int idx = 0; idx < size; idx++) {
			local[idx] = x[idx];
		}
		dx = f(size, local, data);
		for (int idx = 0; idx < size; idx++) {
			adams_dx[i * size + idx] = dx[idx];
		}
		delete[] dx;
		rk->rk_step(h, size, x, f, data);
	}
	delete[] local;
}

void explicit_adams::adams_step(long double h, int size, long double *x, std::function<long double*(int, long double*, void*)> f, void * data) {
	for (int ki = 0; ki < a_size; ki++) {
		for (int idx = 0; idx < size; idx++) {
			x[idx] += h * adams_a[ki] * adams_dx[ki * size + idx];
		}
	}
	for (int i = 0; i < (a_size-1) * size; i++) {
		adams_dx[i] = adams_dx[i+size];
	}
	long double *local = new long double[size];
	for (int i = 0; i < size; i++) {
		local[i] = x[i];
	}
	long double *dx;
	dx = f(size, local, data);
	for (int idx = 0; idx < size; idx++) {
		adams_dx[(a_size-1) * size + idx] = dx[idx];
	}
	delete[] dx;
	delete[] local;
}

