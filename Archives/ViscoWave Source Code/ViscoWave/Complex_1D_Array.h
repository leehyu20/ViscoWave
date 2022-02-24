#pragma once

#include "ap.h"

#ifndef Complex_1D_Array_H
#define Complex_1D_Array_H

class Complex_1D_Array {
private:
	int size;

public:
	alglib::complex_1d_array array_1d;

	Complex_1D_Array(int n);
	void Reset();
	~Complex_1D_Array();
};

#endif