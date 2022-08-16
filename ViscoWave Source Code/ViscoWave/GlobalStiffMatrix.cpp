#include "GlobalStiffMatrix.h"

#include <cmath>
#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "Complex_2D_Array.h"
#include "ViscoWaveEngine.h"

GlobalStiffMatrix::GlobalStiffMatrix(
	int num_pavement_layers,
	int num_prony_elements,
	int num_VE_layer)
{
	number_pavement_layers = num_pavement_layers;
	number_prony_elements = num_prony_elements;
	number_ve_layer = num_VE_layer;
	matrix = new Complex_2D_Array(2 * number_pavement_layers, 2 * number_pavement_layers);
};

void GlobalStiffMatrix::CalculateMatrix(
	alglib::real_1d_array  pavement_structure,
	alglib::real_1d_array  final_prony_values,
	alglib::complex  laplace_var,
	double gaussian_point)
{
	matrix->Reset();
	struct LayerPavementStructure layer_pavement_structure;
	for (int layer_number = 0; layer_number < number_pavement_layers; layer_number++)
	{
		layer_pavement_structure = {
			pavement_structure[5 * layer_number] ,
			pavement_structure[5 * layer_number + 1] ,
			pavement_structure[5 * layer_number + 2] ,
			pavement_structure[5 * layer_number + 3] ,
			pavement_structure[5 * layer_number + 4] ,
		};

		if (layer_number >= 0 && layer_number < number_pavement_layers - 1) {
			if (layer_pavement_structure.shaer_mod == 0) {
				stiff2node_visco(layer_pavement_structure, final_prony_values, laplace_var, gaussian_point, layer_number);
			}
			else {
				stiff2node_elastic(layer_pavement_structure, laplace_var, gaussian_point, layer_number);
			}
		}
		else // (layer_number == (num_pavement_layers - 1))
		{
			if (layer_pavement_structure.thickness == 0) {
				stiff1node_elastic(layer_pavement_structure, laplace_var, gaussian_point, layer_number);
			}
			else {
				stiff1node_elastic_nohalf(layer_pavement_structure, laplace_var, gaussian_point, layer_number);
			}
		}
	}
}

GlobalStiffMatrix::~GlobalStiffMatrix()
{
}

void GlobalStiffMatrix::stiff2node_visco(
	LayerPavementStructure layer_structure,
	alglib::real_1d_array final_prony_values,
	alglib::complex Laplace_Var,
	double gaussian_point,
	int Layer_Number)
{
	alglib::complex Prony_Coeff_E, Prony_Coeff_Rho, Complex_One;
	alglib::complex Hankel_Var, Complex_Poisson, Complex_Density, Complex_Thickness, Complex_Damping;
	alglib::complex Complex_Shear_Mod;
	alglib::complex Complex_Lame;
	alglib::complex C_Wave_Vel;
	alglib::complex S_Wave_Vel;
	alglib::complex Function_F, ehf;
	alglib::complex Function_G, ehg;
	alglib::complex myK, myK1, myQ1, myQ2, myQ3, myQ4, mydenom;
	alglib::complex Ctemp;

	int Num_Pavt_Layers2 = 2 * number_pavement_layers;

	alglib::complex_2d_array c2;

	Complex_One.x = 1;
	Complex_One.y = 0;

	Hankel_Var.x = gaussian_point;
	Hankel_Var.y = 0;

	Complex_Poisson.x = layer_structure.possion;
	Complex_Poisson.y = 0;

	Complex_Density.x = layer_structure.density;
	Complex_Density.y = 0;

	Complex_Thickness.x = layer_structure.thickness;
	Complex_Thickness.y = 0;

	Complex_Damping.x = layer_structure.damping;
	Complex_Damping.y = 0;

	for (int i = 0; i < number_prony_elements; i++) {
		if (i == 0) {
			Complex_Shear_Mod.x = final_prony_values[Layer_Number];
			Complex_Shear_Mod.y = 0;
		}
		else {
			Prony_Coeff_E.x = final_prony_values[(number_ve_layer + 1) * i + Layer_Number];
			Prony_Coeff_E.y = 0;
			Prony_Coeff_Rho.x = final_prony_values[(number_ve_layer + 1) * i + number_ve_layer];
			Prony_Coeff_Rho.y = 0;
			Complex_Shear_Mod += Prony_Coeff_E * Laplace_Var / (Laplace_Var + Complex_One / Prony_Coeff_Rho);
		}
	}
	Complex_Shear_Mod *= (1 + Complex_Damping * Laplace_Var);

	Complex_Lame = 2 * Complex_Shear_Mod * Complex_Poisson / (1 - 2 * Complex_Poisson);
	C_Wave_Vel = ViscoWaveEngine::Compute_Power(((Complex_Lame + 2 * Complex_Shear_Mod) / Complex_Density), 0.5);
	S_Wave_Vel = ViscoWaveEngine::Compute_Power((Complex_Shear_Mod / Complex_Density), 0.5);

	alglib::complex Hankel_Var_Square = Hankel_Var * Hankel_Var;
	alglib::complex Laplace_Var_Square = Laplace_Var * Laplace_Var;
	Function_F = ViscoWaveEngine::Compute_Power((Hankel_Var_Square + (Laplace_Var_Square / (C_Wave_Vel * C_Wave_Vel))), 0.5);
	Function_G = ViscoWaveEngine::Compute_Power((Hankel_Var_Square + (Laplace_Var_Square / (S_Wave_Vel * S_Wave_Vel))), 0.5);

	ehf = ViscoWaveEngine::Compute_Exp(-Complex_Thickness * Function_F);
	ehg = ViscoWaveEngine::Compute_Exp(-Complex_Thickness * Function_G);

	alglib::complex Function_F_Square = Function_F * Function_F;
	alglib::complex Function_G_Square = Function_G * Function_G;
	alglib::complex Function_F_Mul_Function_G = Function_F * Function_G;

	myK = Hankel_Var_Square + Function_G_Square;
	myK1 = Hankel_Var_Square - Function_G_Square;

	alglib::complex ehf_square = ehf * ehf;
	alglib::complex ehg_square = ehg * ehg;
	myQ1 = 1 - ehf_square;
	myQ2 = 1 - ehg_square;
	myQ3 = 1 + ehf_square;
	myQ4 = 1 + ehg_square;
	mydenom = myQ1 * myQ2 * (Function_F_Square * Function_G_Square + Hankel_Var_Square * Hankel_Var_Square)
		- 2 * (myQ3 * myQ4 - 4 * ehf * ehg) * Function_F_Mul_Function_G * Hankel_Var_Square;

	c2.setlength(4, 4);
	c2(0, 0) = Complex_Shear_Mod / mydenom * (-Function_F * myK1 * (myQ1 * myQ4 * Function_F_Mul_Function_G - myQ2 * myQ3 * Hankel_Var_Square));
	c2(0, 1) = Complex_Shear_Mod / mydenom * (-Hankel_Var * (-myQ1 * myQ2 * (2 * Function_F_Square * Function_G_Square + Hankel_Var_Square * myK)
		+ (myQ3 * myQ4 - 4 * ehf * ehg) * Function_F_Mul_Function_G * (Function_G_Square + 3 * Hankel_Var_Square)));
	c2(0, 2) = Complex_Shear_Mod / mydenom * (2 * Function_F * myK1 * (myQ1 * ehg * Function_F_Mul_Function_G - ehf * myQ2 * Hankel_Var_Square));
	c2(0, 3) = Complex_Shear_Mod / mydenom * (-2 * (ehf - ehg) * (1 - ehf * ehg) * Function_F_Mul_Function_G * Hankel_Var * myK1);

	c2(1, 0) = c2(0, 1);
	c2(1, 1) = Complex_Shear_Mod / mydenom * (-Function_G * myK1 * (myQ2 * myQ3 * Function_F_Mul_Function_G - myQ1 * myQ4 * Hankel_Var_Square));
	c2(1, 2) = Complex_Shear_Mod / mydenom * (2 * (ehf - ehg) * (1 - ehf * ehg) * Function_F_Mul_Function_G * Hankel_Var * (Hankel_Var_Square - Function_G_Square));
	c2(1, 3) = Complex_Shear_Mod / mydenom * (-2 * Function_G * myK1 * (-ehf * myQ2 * Function_F_Mul_Function_G + myQ1 * ehg * Hankel_Var_Square));

	c2(2, 0) = c2(0, 2);
	c2(2, 1) = c2(1, 2);
	c2(2, 2) = Complex_Shear_Mod / mydenom * (-Function_F * myK1 * (myQ1 * myQ4 * Function_F_Mul_Function_G - myQ2 * myQ3 * Hankel_Var_Square));
	c2(2, 3) = Complex_Shear_Mod / mydenom * (Hankel_Var * (-myQ1 * myQ2 * (2 * Function_F_Square * Function_G_Square + Hankel_Var_Square * myK)
		+ (myQ3 * myQ4 - 4 * ehf * ehg) * Function_F_Mul_Function_G * (Function_G_Square + 3 * Hankel_Var_Square)));

	c2(3, 0) = c2(0, 3);
	c2(3, 1) = c2(1, 3);
	c2(3, 2) = c2(2, 3);
	c2(3, 3) = Complex_Shear_Mod / mydenom * (-Function_G * myK1 * (myQ2 * myQ3 * Function_F_Mul_Function_G - myQ1 * myQ4 * Hankel_Var_Square));

	//Update Stiffness Matrix
	if (mydenom.x != 0 || mydenom.y != 0)
	{
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				matrix->array_2d(2 * Layer_Number + i, 2 * Layer_Number + j) += c2(i, j);
			}
		}
	}
}

void GlobalStiffMatrix::stiff2node_elastic(
	LayerPavementStructure layer_structure,
	alglib::complex laplace_var,
	double gaussian_point,
	int layer_number)
{
	alglib::complex myK, myK1, myQ1, myQ2, myQ3, myQ4, mydenom;
	alglib::complex Hankel_Var, Complex_Poisson, Complex_Density, Complex_Thickness, Complex_Damping;
	alglib::complex Complex_Shear_Mod;
	alglib::complex Complex_Lame;
	alglib::complex C_Wave_Vel;
	alglib::complex S_Wave_Vel;
	alglib::complex Function_F, ehf;
	alglib::complex Function_G, ehg;
	alglib::complex Ctemp;

	alglib::complex_2d_array c2;

	Hankel_Var.x = gaussian_point;
	Hankel_Var.y = 0;

	Complex_Poisson.x = layer_structure.possion;
	Complex_Poisson.y = 0;

	Complex_Density.x = layer_structure.density;
	Complex_Density.y = 0;

	Complex_Thickness.x = layer_structure.thickness;
	Complex_Thickness.y = 0;

	Complex_Damping.x = layer_structure.damping;
	Complex_Damping.y = 0;

	Complex_Shear_Mod.x = layer_structure.shaer_mod;
	Complex_Shear_Mod.y = 0;
	Complex_Shear_Mod = Complex_Shear_Mod * (1 + Complex_Damping * laplace_var);

	Complex_Lame = 2 * Complex_Shear_Mod * Complex_Poisson / (1 - 2 * Complex_Poisson);
	C_Wave_Vel = ViscoWaveEngine::Compute_Power(((Complex_Lame + 2 * Complex_Shear_Mod) / Complex_Density), 0.5);
	S_Wave_Vel = ViscoWaveEngine::Compute_Power((Complex_Shear_Mod / Complex_Density), 0.5);

	alglib::complex Hankel_Var_Square = Hankel_Var * Hankel_Var;
	Function_F = ViscoWaveEngine::Compute_Power((Hankel_Var * Hankel_Var + ((laplace_var * laplace_var) / (C_Wave_Vel * C_Wave_Vel))), 0.5);
	Function_G = ViscoWaveEngine::Compute_Power((Hankel_Var_Square + ((laplace_var * laplace_var) / (S_Wave_Vel * S_Wave_Vel))), 0.5);

	ehf = ViscoWaveEngine::Compute_Exp(-Complex_Thickness * Function_F);
	ehg = ViscoWaveEngine::Compute_Exp(-Complex_Thickness * Function_G);

	alglib::complex Function_F_Square = Function_F * Function_F;
	alglib::complex Function_G_Square = Function_G * Function_G;
	alglib::complex Function_F_Mul_Function_G = Function_F * Function_G;

	myK = Hankel_Var_Square + Function_G_Square;
	myK1 = Hankel_Var_Square - Function_G_Square;

	alglib::complex ehf_square = ehf * ehf;
	alglib::complex ehg_square = ehg * ehg;
	myQ1 = 1 - ehf_square;
	myQ2 = 1 - ehg_square;
	myQ3 = 1 + ehf_square;
	myQ4 = 1 + ehg_square;
	mydenom = myQ1 * myQ2 * (Function_F_Square * Function_G_Square + Hankel_Var_Square * Hankel_Var_Square) - 2 * (myQ3 * myQ4 - 4 * ehf * ehg) * Function_F_Mul_Function_G * Hankel_Var_Square;

	c2.setlength(4, 4);
	c2(0, 0) = Complex_Shear_Mod / mydenom * (-Function_F * myK1 * (myQ1 * myQ4 * Function_F_Mul_Function_G - myQ2 * myQ3 * Hankel_Var_Square));
	c2(0, 1) = Complex_Shear_Mod / mydenom * (-Hankel_Var * (-myQ1 * myQ2 * (2 * Function_F_Square * Function_G_Square + Hankel_Var_Square * myK)
		+ (myQ3 * myQ4 - 4 * ehf * ehg) * Function_F_Mul_Function_G * (Function_G_Square + 3 * Hankel_Var_Square)));
	c2(0, 2) = Complex_Shear_Mod / mydenom * (2 * Function_F * myK1 * (myQ1 * ehg * Function_F_Mul_Function_G - ehf * myQ2 * Hankel_Var_Square));
	c2(0, 3) = Complex_Shear_Mod / mydenom * (-2 * (ehf - ehg) * (1 - ehf * ehg) * Function_F_Mul_Function_G * Hankel_Var * myK1);

	c2(1, 0) = c2(0, 1);
	c2(1, 1) = Complex_Shear_Mod / mydenom * (-Function_G * myK1 * (myQ2 * myQ3 * Function_F_Mul_Function_G - myQ1 * myQ4 * Hankel_Var_Square));
	c2(1, 2) = Complex_Shear_Mod / mydenom * (2 * (ehf - ehg) * (1 - ehf * ehg) * Function_F_Mul_Function_G * Hankel_Var * (Hankel_Var_Square - Function_G_Square));
	c2(1, 3) = Complex_Shear_Mod / mydenom * (-2 * Function_G * myK1 * (-ehf * myQ2 * Function_F_Mul_Function_G + myQ1 * ehg * Hankel_Var_Square));

	c2(2, 0) = c2(0, 2);
	c2(2, 1) = c2(1, 2);
	c2(2, 2) = Complex_Shear_Mod / mydenom * (-Function_F * myK1 * (myQ1 * myQ4 * Function_F_Mul_Function_G - myQ2 * myQ3 * Hankel_Var_Square));
	c2(2, 3) = Complex_Shear_Mod / mydenom * (Hankel_Var * (-myQ1 * myQ2 * (2 * Function_F_Square * Function_G_Square + Hankel_Var_Square * myK)
		+ (myQ3 * myQ4 - 4 * ehf * ehg) * Function_F_Mul_Function_G * (Function_G_Square + 3 * Hankel_Var_Square)));

	c2(3, 0) = c2(0, 3);
	c2(3, 1) = c2(1, 3);
	c2(3, 2) = c2(2, 3);
	c2(3, 3) = Complex_Shear_Mod / mydenom * (-Function_G * myK1 * (myQ2 * myQ3 * Function_F_Mul_Function_G - myQ1 * myQ4 * Hankel_Var_Square));

	if (mydenom.x != 0 || mydenom.y != 0)
	{
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				matrix->array_2d(2 * layer_number + i, 2 * layer_number + j) += c2(i, j);
			}
		}
	}
}

void GlobalStiffMatrix::stiff1node_elastic(
	LayerPavementStructure layer_structure,
	alglib::complex laplace_var,
	double gaussian_point,
	int layer_number)
{
	alglib::complex myK, mydenom;
	alglib::complex Hankel_Var, Complex_Poisson, Complex_Density, Complex_Thickness, Complex_Damping;
	alglib::complex Complex_Shear_Mod;
	alglib::complex Complex_Lame;
	alglib::complex C_Wave_Vel;
	alglib::complex S_Wave_Vel;
	alglib::complex Function_F;
	alglib::complex Function_G;
	alglib::complex Ctemp;

	alglib::complex_2d_array c2;

	Hankel_Var.x = gaussian_point;
	Hankel_Var.y = 0;

	Complex_Poisson.x = layer_structure.possion;
	Complex_Poisson.y = 0;

	Complex_Density.x = layer_structure.density;
	Complex_Density.y = 0;

	Complex_Thickness.x = layer_structure.thickness;
	Complex_Thickness.y = 0;

	Complex_Damping.x = layer_structure.damping;
	Complex_Damping.y = 0;

	Complex_Shear_Mod.x = layer_structure.shaer_mod;
	Complex_Shear_Mod.y = 0;
	Complex_Shear_Mod = Complex_Shear_Mod * (1 + Complex_Damping * laplace_var);

	Complex_Lame = 2 * Complex_Shear_Mod * Complex_Poisson / (1 - 2 * Complex_Poisson);
	C_Wave_Vel = ViscoWaveEngine::Compute_Power(((Complex_Lame + 2 * Complex_Shear_Mod) / Complex_Density), 0.5);
	S_Wave_Vel = ViscoWaveEngine::Compute_Power((Complex_Shear_Mod / Complex_Density), 0.5);

	alglib::complex Hankel_Var_Square = Hankel_Var * Hankel_Var;
	Function_F = ViscoWaveEngine::Compute_Power((Hankel_Var_Square + ((laplace_var * laplace_var) / (C_Wave_Vel * C_Wave_Vel))), 0.5);
	Function_G = ViscoWaveEngine::Compute_Power((Hankel_Var_Square + ((laplace_var * laplace_var) / (S_Wave_Vel * S_Wave_Vel))), 0.5);

	alglib::complex Function_G_Square = Function_G * Function_G;
	alglib::complex Function_F_Mul_Function_G = Function_F * Function_G;

	myK = (Hankel_Var_Square + Function_G_Square);
	mydenom = (Hankel_Var_Square - Function_F_Mul_Function_G);

	c2.setlength(2, 2);
	c2(0, 0) = Complex_Shear_Mod / mydenom * (2 * Hankel_Var_Square * Function_F - Function_F * myK);
	c2(0, 1) = Complex_Shear_Mod / mydenom * (-2 * Hankel_Var * Function_F_Mul_Function_G + Hankel_Var * myK);
	c2(1, 0) = c2(0, 1);
	c2(1, 1) = Complex_Shear_Mod / mydenom * (2 * Hankel_Var_Square * Function_G - Function_G * myK);

	if (mydenom.x != 0 || mydenom.y != 0)
	{
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < 2; j++) {
				matrix->array_2d(2 * layer_number + i, 2 * layer_number + j) += c2(i, j);
			}
		}
	}
}

void GlobalStiffMatrix::stiff1node_elastic_nohalf(
	LayerPavementStructure layer_structure,
	alglib::complex laplace_var,
	double gaussian_point,
	int layer_number)
{
	alglib::complex myK, myK1, myQ1, myQ2, myQ3, myQ4, mydenom;
	alglib::complex Hankel_Var, Complex_Poisson, Complex_Density, Complex_Thickness, Complex_Damping;
	alglib::complex Complex_Shear_Mod;
	alglib::complex Complex_Lame;
	alglib::complex C_Wave_Vel;
	alglib::complex S_Wave_Vel;
	alglib::complex Function_F, ehf;
	alglib::complex Function_G, ehg;
	alglib::complex Ctemp;

	alglib::complex_2d_array c2;

	Hankel_Var.x = gaussian_point;
	Hankel_Var.y = 0;

	Complex_Poisson.x = layer_structure.possion;
	Complex_Poisson.y = 0;

	Complex_Density.x = layer_structure.density;
	Complex_Density.y = 0;

	Complex_Thickness.x = layer_structure.thickness;
	Complex_Thickness.y = 0;

	Complex_Damping.x = layer_structure.damping;
	Complex_Damping.y = 0;

	Complex_Shear_Mod.x = layer_structure.shaer_mod;
	Complex_Shear_Mod.y = 0;
	Complex_Shear_Mod = Complex_Shear_Mod * (1 + Complex_Damping * laplace_var);

	Complex_Lame = 2 * Complex_Shear_Mod * Complex_Poisson / (1 - 2 * Complex_Poisson);
	C_Wave_Vel = ViscoWaveEngine::Compute_Power(((Complex_Lame + 2 * Complex_Shear_Mod) / Complex_Density), 0.5);
	S_Wave_Vel = ViscoWaveEngine::Compute_Power((Complex_Shear_Mod / Complex_Density), 0.5);

	alglib::complex Hankel_Var_Square = Hankel_Var * Hankel_Var;
	Function_F = ViscoWaveEngine::Compute_Power((Hankel_Var_Square + ((laplace_var * laplace_var) / (C_Wave_Vel * C_Wave_Vel))), 0.5);
	Function_G = ViscoWaveEngine::Compute_Power((Hankel_Var_Square + ((laplace_var * laplace_var) / (S_Wave_Vel * S_Wave_Vel))), 0.5);

	ehf = ViscoWaveEngine::Compute_Exp(-Complex_Thickness * Function_F);
	ehg = ViscoWaveEngine::Compute_Exp(-Complex_Thickness * Function_G);

	alglib::complex Function_F_Square = Function_F * Function_F;
	alglib::complex Function_G_Square = Function_G * Function_G;
	alglib::complex Function_F_Mul_Function_G = Function_F * Function_G;

	myK = (Hankel_Var_Square + Function_G_Square);
	myK1 = (Hankel_Var_Square - Function_G_Square);

	alglib::complex ehf_square = ehf * ehf;
	alglib::complex ehg_square = ehg * ehg;
	myQ1 = 1 - ehf_square;
	myQ2 = 1 - ehg_square;
	myQ3 = 1 + ehf_square;
	myQ4 = 1 + ehg_square;
	mydenom = myQ1 * myQ2 * (Function_F_Square * Function_G_Square + Hankel_Var_Square * Hankel_Var_Square) - 2 * (myQ3 * myQ4 - 4 * ehf * ehg) * Function_F_Mul_Function_G * Hankel_Var_Square;

	c2.setlength(4, 4);
	c2(0, 0) = Complex_Shear_Mod / mydenom * (-Function_F * myK1 * (myQ1 * myQ4 * Function_F_Mul_Function_G - myQ2 * myQ3 * Hankel_Var_Square));
	c2(0, 1) = Complex_Shear_Mod / mydenom * (-Hankel_Var * (-myQ1 * myQ2 * (2 * Function_F_Square * Function_G_Square + Hankel_Var_Square * myK) + (myQ3 * myQ4 - 4 * ehf * ehg) * Function_F_Mul_Function_G * (Function_G_Square + 3 * Hankel_Var_Square)));

	c2(1, 0) = c2(0, 1);
	c2(1, 1) = Complex_Shear_Mod / mydenom * (-Function_G * myK1 * (myQ2 * myQ3 * Function_F_Mul_Function_G - myQ1 * myQ4 * Hankel_Var_Square));

	if (mydenom.x != 0 || mydenom.y != 0)
	{
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < 2; j++) {
				matrix->array_2d(2 * layer_number + i, 2 * layer_number + j) += c2(i, j);
			}
		}
	}
}