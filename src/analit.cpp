#include <cmath>
#include <iostream>

long double pi = 2 * acos(0.0);
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
struct keplerian {
    long double a;
	long double e;
	long double w;
	long double om;
	long double incl;
	long double M;
};

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
// 	std::cout << "DEBUG " << (n[0]*e[0]+n[1]*e[1]+n[2]*e[2]) / veclen(e) / veclen(n) << "\n";

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
long double E_t;
long double* kepler_to_cartesian(
	long double a,
	long double e,
	long double w,
	long double om,
	long double incl,
	long double M,
	long double t,
	long double mu,
	long double h) {
    long double M_t = M;
    if (t != 0) {
        M_t += t * sqrt(mu / a / a / a);
    }

    while (M_t < 0) {
        M_t += 2 * pi;
    }
    while (M_t > 2 * pi) {
        M_t -= 2 * pi;
    }
    if (t_counter == 0)  {
        E_t = M_t;
    } else {
        E_t += M_t;
        E_t -= (t_counter) * sqrt(mu / a / a / a);
    }

	while (t_counter <= t) {
		E_t = E_t - (E_t - e * sinl(E_t) - M_t) / (1 - e * cosl(E_t));
		t_counter += 1;
	}
	long double v_t = 2 * arctan2(
		sqrt(1+e)*sinl(E_t/2), sqrt(1-e)*cosl(E_t/2)
	);

	long double rc_t = a * (1 - e * cosl(E_t));

	long double *o_t = new long double[3]{ cosl(v_t), sinl(v_t), 0.0 };
	for (int k = 0; k < 3; k++) {
	    o_t[k] *= rc_t;
	}

	long double *odt_t = new long double[3]{ -sinl(E_t), sqrt(1-e*e)*cosl(E_t), 0 };
	for (int k = 0; k < 3; k++) {
		odt_t[k] *= sqrt(mu * a) / rc_t;
	}

	long double* x = new long double[3] {
		o_t[0] * (cosl(w)*cosl(om)-sinl(w)*cosl(incl)*sinl(om)) - o_t[1] * (sinl(w)*cosl(om) + cosl(w)*cosl(incl)*sinl(om)),
		o_t[0] * (cosl(w)*sinl(om)+sinl(w)*cosl(incl)*cosl(om)) + o_t[1] * (cosl(w)*cosl(incl)*cosl(om) - sinl(w)*sinl(om)),
		o_t[0] * (sinl(w)*sin(incl))                            + o_t[1] * (cosl(w)*sinl(incl))
	};
    delete[] odt_t;
	return x;
}

int main() {
    std::cout.setf(std::ios::fixed);  // вывод в фиксированном формате
	std::cout.precision(50);      // вывод до 6 знака после точки, включительно
    long double *x = new long double[6] {
        20056834.34185719862580299377441406250000000000000000000000, 133706925.45214439928531646728515625000000000000000000000000, 57994342.94798526912927627563476562500000000000000000000000, -29.94808362540784100851165572265699665877036750316620, 3.86421640889615978196880430317605714662931859493256, 1.67567000626789382167463465789225551816343795508146
    };
    long double mu = 1.98847e30 * 6.67430e-20;
    struct keplerian kepl = cartesian_to_kepler(x, mu);
    long double t = 0;
    long double T = 31536000.;// / 16;
    long double h = 1 * 50000;
    while (t < T) {
        long double *r = kepler_to_cartesian(
            kepl.a,
        	kepl.e,
        	kepl.w,
            kepl.om,
            kepl.incl,
        	kepl.M,
            t,
            mu,
			h
        );
        std::cout << t << " ";
        for (int i = 0; i < 3; i++) {
            std::cout << r[i]<< " ";
        }
        std::cout << std::endl;
        t += h;
        delete[] r;
    }

    return 0;
}
