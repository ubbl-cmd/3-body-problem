	#ifndef EXPLICIT_ADAMS
#define EXPLICIT_ADAMS

#include <functional>
#include "explicit_rk.hpp"

class explicit_adams {
public:
	int a_size;
	long double *adams_a;
	long double *adams_dx;

	void adams_step(long double h, int size, long double *x, std::function<long double*(int, long double*, void*)> f, void * data);

	void razgonka(long double h, explicit_rk *rk, int size, long double *x, std::function<long double*(int, long double*, void*)> f, void *data);

	explicit_adams(int a_size, long double *adams_a)
	:a_size(a_size), adams_a(adams_a) {};
};

#endif
