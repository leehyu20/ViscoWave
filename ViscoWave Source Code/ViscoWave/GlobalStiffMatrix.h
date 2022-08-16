#pragma once
#include "ap.h"
#include "Complex_2D_Array.h"

#ifndef GlobalStiffMatrix_H
#define GlobalStiffMatrix_H

struct LayerPavementStructure {
	double shaer_mod;
	double possion;
	double density;
	double thickness;
	double damping;
};

class GlobalStiffMatrix {
private:
	int number_pavement_layers;
	int number_prony_elements;
	int number_ve_layer;

	void stiff2node_visco(
		LayerPavementStructure layer_structure,
		alglib::real_1d_array final_prony_values,
		alglib::complex laplace_var,
		double gaussian_point,
		int layer_number);
	void stiff2node_elastic(
		LayerPavementStructure layer_structure,
		alglib::complex laplace_var,
		double gaussian_point,
		int layer_number);
	void stiff1node_elastic(
		LayerPavementStructure layer_structure,
		alglib::complex laplace_var,
		double gaussian_point,
		int layer_number);
	void stiff1node_elastic_nohalf(
		LayerPavementStructure layer_structure,
		alglib::complex laplace_var,
		double gaussian_point,
		int layer_number);

public:
	Complex_2D_Array* matrix;

	GlobalStiffMatrix(
		int num_pavement_layers,
		int num_prony_elements,
		int num_VE_layer);
	void CalculateMatrix(
		alglib::real_1d_array  pavement_structure,
		alglib::real_1d_array  final_prony_values,
		alglib::complex  laplace_var,
		double gaussian_point);
	~GlobalStiffMatrix();
};
#endif