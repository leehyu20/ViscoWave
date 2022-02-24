#include "Complex_1D_Array.h"

Complex_1D_Array::Complex_1D_Array(int n) {
	size = n;
	array_1d.setlength(size);
	Reset();
}

void Complex_1D_Array::Reset() {
	for (int i = 0; i < size; i++) {
		array_1d[i].x = 0;
		array_1d[i].y = 0;
	}
}

Complex_1D_Array::~Complex_1D_Array()
{
}

