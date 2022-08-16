#pragma once

#include "ap.h"
#include "linalg.h"
#include "solvers.h"
#include "specialfunctions.h"
#include "fasttransforms.h"
#include <cmath>
#include <omp.h>
#include "stdafx.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "interpolation.h"
#include "Complex_1D_Array.h"

#ifndef ViscoWaveEngine_H
#define ViscoWaveEngine_H

class ViscoWaveEngine
{
public:
	static alglib::real_1d_array bessel_J1_root() {
		alglib::real_1d_array value;
		value.setlength(6);
		value[0] = 0;
		value[1] = 3.83170597020751;
		value[2] = 7.01558666981562;
		value[3] = 10.1734681350627;
		value[4] = 13.3236919363142;
		value[5] = 16.4706300508776;
		return value;
	}

	static alglib::real_1d_array gaussian_nodes() {
		alglib::real_1d_array value;
		value.setlength(6);
		value[0] = -0.932469514203152;
		value[1] = -0.661209386466265;
		value[2] = -0.238619186083197;
		value[3] = 0.238619186083197;
		value[4] = 0.661209386466265;
		value[5] = 0.932469514203152;
		return value;
	}

	static alglib::real_1d_array gaussian_weight() {
		alglib::real_1d_array value;
		value.setlength(6);
		value[0] = 0.17132449237917;
		value[1] = 0.360761573048139;
		value[2] = 0.467913934572691;
		value[3] = 0.467913934572691;
		value[4] = 0.360761573048139;
		value[5] = 0.17132449237917;
		return value;
	}

	static alglib::complex Compute_Power(alglib::complex a, double b) {
		double abs, argument;
		alglib::complex temp, power;

		abs = sqrt(pow(a.x, 2) + pow(a.y, 2));
		argument = atan2(a.y, a.x);

		temp.x = cos(argument * b);
		temp.y = sin(argument * b);

		power = pow(abs, b) * temp;
		return power;
	};

	static alglib::complex Compute_Exp(alglib::complex a)
	{
		alglib::complex temp, exp_value;

		temp.x = cos(a.y);
		temp.y = sin(a.y);

		exp_value = exp(a.x) * temp;
		return exp_value;
	};

	static alglib::complex_1d_array Matrix_Mult_Vec(
		alglib::complex_2d_array matrix,
		alglib::complex_1d_array vector,
		int n)
	{
		Complex_1D_Array empty_arr(n);
		alglib::complex_1d_array multipled = empty_arr.array_1d;

		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				multipled[i] = multipled[i] + matrix(i, j) * vector[j];
			}
		}
		return multipled;
	};
};

#endif