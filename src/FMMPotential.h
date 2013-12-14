#ifndef Potential_h
#define Potential_h

#include "FMMBox.h"

class Potential {

public:

	/// The Truncation number
	int degree;

	/// Constructor
	Potential(const int degree) : degree(degree) 
	{

	}

	/// Directly evaluate the potential
	inline double DirectEvaluate(const Complex& y, const Complex& x) 
	{
		return real(log(y - x));
	}

	/// Matrix-vector multiplication between translation matrix and vector of expansion coefficients
	inline ComplexVec ApplyTranslation(const ComplexMat& matrix, const ComplexVec& coeff) 
	{
		ComplexVec product(degree, Complex(0,0));
		for (int i = 0; i < degree; i++)
			for (int j = 0; j < degree; j++)
				product[i] += coeff[j] * matrix[i][j];
		return product;
	}

	/// Apply Multipole-to-Local translation
	inline ComplexVec MultipoleToLocal(const Complex& from, const Complex& to, const ComplexVec& MultipoleCoeff) 
	{
		Complex t = to - from;
		ComplexMat M2L(degree, ComplexVec(degree, Complex(0,0)));
		M2L[0][0] = log(t);
		M2L[1][0] = 1.0 / t;
		for (int i = 2; i < degree; i++)
			M2L[i][0] = -M2L[i-1][0] * (double)(i - 1) / (t * (double)i);
		for (int j = 1; j < degree; j++)
			M2L[0][j] = 1.0 / pow(t, j);
		for (int i = 1; i < degree; i++)
			for (int j = 1; j < degree; j++)
				M2L[i][j] = M2L[i-1][j] * (double)(i + j - 1) / (-t * (double)i);
		return ApplyTranslation(M2L, MultipoleCoeff);
	}

	/// Apply Multipole-to-Multipole translation
	inline ComplexVec MultipoleToMultipole(const Complex& from, const Complex& to, const ComplexVec& MultipoleCoeff) 
	{
		Complex t = to - from;
		ComplexMat M2M(degree, ComplexVec(degree, Complex(0,0)));
		for (int i = 0; i < degree; i++)
			M2M[i][i] = Complex(1, 0);
		M2M[1][0] = t;
		for (int i = 2; i < degree; i++)
			M2M[i][0] = -M2M[i-1][0] * (double)(i-1) * t / (double)i;
		for (int i = 1; i < degree; i++)
			for (int j = i-1; j >= 1; j--)
				M2M[i][j] = -M2M[i][j+1] * t * (double)j / (double)(i - j);
		return ApplyTranslation(M2M, MultipoleCoeff);
	}

	/// Apply Local-to-Local translation
	inline ComplexVec LocalToLocal(const Complex& from, const Complex& to, const ComplexVec& LocalCoeff) 
	{
		Complex t = to - from;
		ComplexMat L2L(degree, ComplexVec(degree, Complex(0,0)));
		for (int i = 0; i < degree; i++)
			L2L[i][i] = Complex(1, 0);
		for (int j = 1; j < degree; j++)
			L2L[0][j] = L2L[0][j-1] * t;
		for (int i = 1; i < degree; i++)
			for (int j = i+1; j < degree; j++)
				L2L[i][j] = L2L[i-1][j] * (double)(j - i + 1) / (t * (double)i);
		return ApplyTranslation(L2L, LocalCoeff);
	}

	/// Get local expansion coefficients
	inline ComplexVec GetLocalCoeffs(const Complex& y, const Complex& x_star) 
	{
		ComplexVec coeffs(degree, Complex(0,0));
		for (int i = 0; i < degree; i++)
			coeffs[i] = pow(y - x_star, i);
		return coeffs;
	}

	/// Get multipole expansion coefficients
	inline ComplexVec GetMultipoleCoeffs(const Complex& x_i, const Complex& x_star) 
	{
		ComplexVec coeffs(degree, Complex(0,0));
		coeffs[0] = Complex(1, 0);
		for (int i = 1; i < degree; i++)
			coeffs[i] = -pow(x_i - x_star, i) / (double)i;
		return coeffs; 
	}

};

#endif