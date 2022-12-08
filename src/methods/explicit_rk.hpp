#ifndef EXPLICIT_RK
#define EXPLICIT_RK

#include <functional>

class explicit_rk {
public:
	int a_size;
	long double* a;
	long double* b;
	void rk_step(long double h, int size, long double *x, std::function<long double*(int, long double*, void*)> f, void * data);
	explicit explicit_rk(int a_size, long double *a, long double *b)
	:a_size(a_size), a(a), b(b) {};

};

#endif
