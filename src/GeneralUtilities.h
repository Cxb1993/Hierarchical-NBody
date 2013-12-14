#ifndef GeneralUtilities_h
#define GeneralUtilities_h

#include <cmath>
#include <vector>
#include <complex>

#include <random>
#include <ctime>

typedef std::complex<double> Complex;
typedef std::vector<std::complex<double>> ComplexVec;
typedef std::vector<std::vector<std::complex<double>>> ComplexMat;

inline void operator+=(ComplexVec &v1, const ComplexVec& v2) {
	for (size_t i = 0; i < v1.size(); ++i) {
		v1[i] += v2[i];
	}
}

inline double randf() {
	return (double)rand() / RAND_MAX;
}

inline void PrintComplexVec(const ComplexVec& vec, bool newlines) {
	if (newlines)
		for (int i = 0; i < vec.size(); i++)
			printf("[%d]: (%06.4f, %06.4f)\n", i, real(vec[i]), imag(vec[i]));
	else
		for (int i = 0; i < vec.size(); i++)
			printf("(%06.4f, %06.4f) ", real(vec[i]), imag(vec[i]));
}

inline double MaxAbsError(const std::vector<double>& approx, const std::vector<double>& exact) {
    double error = 0;
    for (int i = 0; i < approx.size(); i++) {
    	double rel = fabs(approx[i] - exact[i]);
        if (rel > error)
            error = rel;
    }
    return error;
}

inline double MaxRelError(const std::vector<double>& approx, const std::vector<double>& exact) {
    double error = 0;
    for (int i = 0; i < approx.size(); i++) {
    	double rel = fabs(approx[i] - exact[i]) / fabs(exact[i]);
        if (rel > error)
            error = rel;
    }
    return error;
}

inline double AvgAbsError(const std::vector<double>& approx, const std::vector<double>& exact) {
    double error = 0;
    for (int i = 0; i < approx.size(); i++) {
    	error += fabs(approx[i] - exact[i]);
    }
    return error / approx.size();
}

inline double AvgRelError(const std::vector<double>& approx, const std::vector<double>& exact) {
    double error = 0;
    for (int i = 0; i < approx.size(); i++) {
    	error += fabs(approx[i] - exact[i]) / fabs(exact[i]);
        
    }
    return error / approx.size();
}

inline int set_bit(const int n, const int index, const int target) {
	if (target == 1) return n | (1 << index);
	else             return n & ~(1 << index);
}

inline int get_bit(const int n, const int index) {
	return ((n & (1 << index)) != 0) ? 1 : 0;
}

inline int interleave(const int x, const int y, const int level) {
	if (x == 0 && y == 0) {
		return 0;
	} else {
		int rval = 0;
		for (int i = 0; i < level; i++) {
			rval = set_bit(rval, (level - i) * 2 - 1, get_bit(x, level - i - 1));
			rval = set_bit(rval, (level - i) * 2 - 2, get_bit(y, level - i - 1));
		}
		return rval;
	}
}

inline Complex uninterleave(const int n, const int level) {
	if (n == 0) {
		return Complex(0,0);
	} else {
		int x = 0, y = 0;
		for (int i = 0; i < level; i++) {
			x = set_bit(x, level - i - 1, get_bit(n, (level - i - 1) * 2 + 1));
			y = set_bit(y, level - i - 1, get_bit(n, (level - i - 1) * 2));
		}
		return Complex(x, y);
	}
}


#endif