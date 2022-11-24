#ifndef EXPLICIT_RK
#define EXPLICIT_RK

#include <functional>

class explicit_rk {
public:
	int a_size;
	double* a;
	double* b;
	void rk_step(double h, int size, double *x, std::function<double*(int, double*, void*)> f, void * data);
	explicit explicit_rk(int a_size, double *a, double *b)
	:a_size(a_size), a(a), b(b) {};

};

#endif
