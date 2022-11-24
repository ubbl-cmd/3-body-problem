#ifndef EXPLICIT_ADAMS
#define EXPLICIT_ADAMS

#include <functional>
#include "explicit_rk.hpp"

class explicit_adams {
public:
	int a_size;
	double *adams_a;
	double *adams_dx;

	void adams_step(double h, int size, double *x, std::function<double*(int, double*, void*)> f, void * data);

	void razgonka(double h, explicit_rk *rk, int size, double *x, std::function<double*(int, double*, void*)> f, void *data);

	explicit_adams(int a_size, double *adams_a)
	:a_size(a_size), adams_a(adams_a) {};
};

#endif
