#include "explicit_adams.hpp"
#include <iostream>

void explicit_adams::razgonka(double h, explicit_rk *rk, int size, double *x, std::function<double*(int, double*, void*)> f, void *data) {
	adams_dx = new double[a_size * size];
	double *local = new double[size];
	double *dx;
	for (int i = 0; i < a_size; i++) {
		rk->rk_step(h, size, x, f, data);
		for (int idx = 0; idx < size; idx++) {
			local[idx] = x[idx];
		}
		dx = f(size, local, data);
		for (int idx = 0; idx < size; idx++) {
			adams_dx[i * size + idx] = dx[idx];
		}
		delete[] dx;
	}
	delete[] local;
}

void explicit_adams::adams_step(double h, int size, double *x, std::function<double*(int, double*, void*)> f, void * data) {
	for (int ki = 0; ki < a_size; ki++) {
		for (int idx = 0; idx < size; idx++) {
			x[idx] += h * adams_a[ki] * adams_dx[ki * size + idx];
		}
	}
	for (int i = 0; i < (a_size-1) * size; i++) {
		adams_dx[i] = adams_dx[i+size];
	}
	double *local = new double[size];
	for (int i = 0; i < size; i++) {
		local[i] = x[i];
	}
	double *dx;
	dx = f(size, local, data);
	for (int idx = 0; idx < size; idx++) {
		adams_dx[(a_size-1) * size + idx] = dx[idx];
	}
	delete[] dx;
	delete[] local;
}

