// ViscoWave.cpp -- Interface for VBA ViscoWave functions.
// Built through Visual Studio 2019 (v142) with ISO C++ 14 Standard.

#include "ViscoWave.h"
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

#include "ViscoWaveEngine.h"
#include "GlobalStiffMatrix.h"

int ViscoWave(double* displacement,
	double* Sigmoid,
	double* Pavement,
	double Load_Pressure,
	double Load_Radius,
	double* Sensor_Location,
	double* Time,
	double* Timehistory,
	double dt,
	int num_prony_elements,
	int Num_Pavt_Layers,
	int Num_Sensors,
	int Num_Time,
	int Num_VE_Layer)
{
	alglib::real_1d_array final_prony_values, sigmoid_coeffs, pavement_structure, defl_sensor_location, output_time, fwd_load_history, exact_VW_time_vector, exact_VW_defl_vector;
	alglib::real_1d_array bessel_J1_root, gaussian_nodes, gaussian_weight, temporary_displ, temporary_displ1;
	alglib::real_2d_array output_disp_matrix, exact_VW_defl_matrix;
	alglib::ae_int_t prec_dec_digit, i, j, m, mm, mmm, jj, nn, jk, info, total_num_exact_time, num_of_subregions, num_of_J1_roots;

	double a, b, gaussian_point, fixed_talbot_r, fixed_talbot_sumf, fixed_talbot_theta, y, spline_bound_l, spline_bound_r;
	alglib::complex imaginary_i, fixed_talbot_s, Defl_fixed_talbot_s, ftemp, Defl_fixed_talbot_theta, fixed_talbot_gamma, detGSM, detGSM22;
	alglib::spline1dinterpolant mysplinestruc;
	alglib::lsfitreport rep;

	int num_points_for_prony_fitting;
	alglib::real_1d_array time_for_prony_fitting;
	alglib::real_1d_array relaxation_times, sigmoid_modulus, fitted_prony_series;
	alglib::real_2d_array matrix_for_prony_fit;

	// Define contstants

	// First 6 roots of Bessel Function J1
	bessel_J1_root = ViscoWaveEngine::bessel_J1_root();

	// 6 point Gaussian nodes and weights
	gaussian_nodes = ViscoWaveEngine::gaussian_nodes();

	gaussian_weight = ViscoWaveEngine::gaussian_weight();

	// Times for exact VW evaluation (i.e., for impulse response)
	// Time interval for FWD time history is hard coded to 0.0002 second.

	if (Time[Num_Time - 1] <= 0.071) // For FWD Loading, set tmax <= 0.07 sec
	{
		total_num_exact_time = 29;
		exact_VW_defl_vector.setlength(29);
		exact_VW_time_vector.setlength(29);
		exact_VW_time_vector[0] = 0.0000000001;
		exact_VW_time_vector[1] = 0.0002;
		for (i = 2; i < 6; i++) {
			exact_VW_time_vector[i] = exact_VW_time_vector[i - 1] + 0.0002;
		}
		for (i = 6; i < 10; i++) {
			exact_VW_time_vector[i] = exact_VW_time_vector[i - 1] + 0.0005;
		}
		for (i = 10; i < 15; i++) {
			exact_VW_time_vector[i] = exact_VW_time_vector[i - 1] + 0.001;
		}
		for (i = 15; i < 21; i++) {
			exact_VW_time_vector[i] = exact_VW_time_vector[i - 1] + 0.002;
		}
		for (i = 21; i < 25; i++) {
			exact_VW_time_vector[i] = exact_VW_time_vector[i - 1] + 0.005;
		}
		for (i = 25; i < 27; i++) {
			exact_VW_time_vector[i] = exact_VW_time_vector[i - 1] + 0.01;
		}
		for (i = 27; i < 28; i++) {
			exact_VW_time_vector[i] = exact_VW_time_vector[i - 1] + 0.0098;
		}
		for (i = 28; i < 29; i++) {
			exact_VW_time_vector[i] = exact_VW_time_vector[i - 1] + 0.0002;
		}
	}
	else {
		// Need to notify user for error
	}

	// Initialize Exact VW output vector
	exact_VW_defl_matrix.setlength(total_num_exact_time, Num_Sensors);
	for (i = 0; i < total_num_exact_time; i++) {
		for (j = 0; j < Num_Sensors; j++) {
			exact_VW_defl_matrix(i, j) = 0;
		}
	}

	//Initialize Variables
	output_disp_matrix.setlength(Num_Time, Num_Sensors);
	for (i = 0; i < Num_Time; i++) {
		for (j = 0; j < Num_Sensors; j++) {
			output_disp_matrix(i, j) = 0;
		}
	}

	sigmoid_coeffs.setlength(4 * num_prony_elements);
	for (i = 0; i < 4 * num_prony_elements; i++) {
		sigmoid_coeffs[i] = Sigmoid[i];
	}

	
	pavement_structure.setlength(5 * Num_Pavt_Layers);
	for (i = 0; i < 5 * Num_Pavt_Layers; i++) {
		pavement_structure[i] = Pavement[i];
	}

	defl_sensor_location.setlength(Num_Sensors);
	for (i = 0; i < Num_Sensors; i++) {
		defl_sensor_location[i] = Sensor_Location[i];
	}

	output_time.setlength(Num_Time);
	for (i = 0; i < Num_Time; i++) {
		output_time[i] = Time[i];
	}

	fwd_load_history.setlength(Num_Time);
	for (i = 0; i < Num_Time; i++) {
		fwd_load_history[i] = Timehistory[i];
	}

	//Initialize Intermediate Variables and Calculate Prony Coefficients

	temporary_displ.setlength(Num_Time);
	temporary_displ1.setlength(2 * Num_Time - 1);

	imaginary_i.x = 0;
	imaginary_i.y = 1;

	prec_dec_digit = 16;

	num_points_for_prony_fitting = 671;
	num_prony_elements = 15;

	time_for_prony_fitting.setlength(num_points_for_prony_fitting);
	sigmoid_modulus.setlength(num_points_for_prony_fitting);
	fitted_prony_series.setlength(num_prony_elements);
	relaxation_times.setlength(num_prony_elements);
	matrix_for_prony_fit.setlength(num_points_for_prony_fitting, num_prony_elements);
	final_prony_values.setlength((Num_VE_Layer + 1) * num_prony_elements);

	for (i = 0; i < num_prony_elements; i++) {
		if (i == 0) {
			relaxation_times[i] = 0;
		}
		else {
			relaxation_times[i] = pow(10.0, 8 - i);
		}
		final_prony_values[(Num_VE_Layer + 1) * i + Num_VE_Layer] = relaxation_times[i];
	}

	for (j = 0; j < Num_VE_Layer; j++) {
		for (i = 0; i < num_points_for_prony_fitting; i++) {
			time_for_prony_fitting[i] = pow(10, -6.7 + i * 0.02);
			sigmoid_modulus[i] = pow(10, sigmoid_coeffs[0 + 4 * j] + sigmoid_coeffs[1 + 4 * j] / (1 + exp(sigmoid_coeffs[2 + 4 * j] + sigmoid_coeffs[3 + 4 * j] * log10(time_for_prony_fitting[i]))));

			for (jj = 0; jj < num_prony_elements; jj++) {
				if (jj == 0) {
					matrix_for_prony_fit(i, jj) = 1.0;
				}
				else {
					matrix_for_prony_fit(i, jj) = exp(-time_for_prony_fitting[i] / relaxation_times[jj]);
				}
			}
		}
		alglib::lsfitlinear(sigmoid_modulus, matrix_for_prony_fit, num_points_for_prony_fitting, num_prony_elements, info, fitted_prony_series, rep);

		for (i = 0; i < num_prony_elements; i++) {
			final_prony_values[(Num_VE_Layer + 1) * i + j] = fitted_prony_series[i];
		}
	}


//Start looping for solution
#pragma omp parallel for private(jj, fixed_talbot_r, nn, fixed_talbot_theta, fixed_talbot_s, Defl_fixed_talbot_s, fixed_talbot_sumf, Defl_fixed_talbot_theta, ftemp, fixed_talbot_gamma, y, jk, m, mm, a, b, mmm, gaussian_point, num_of_subregions, num_of_J1_roots) firstprivate(info) shared(dt, Load_Radius, Load_Pressure, bessel_J1_root, gaussian_nodes, exact_VW_defl_matrix, exact_VW_time_vector, total_num_exact_time, gaussian_weight, defl_sensor_location, output_disp_matrix, Num_Time, prec_dec_digit, imaginary_i, Num_Pavt_Layers, pavement_structure, final_prony_values, num_prony_elements, output_time, Num_Sensors, Num_VE_Layer)
	for (jj = 0; jj < total_num_exact_time; jj++) {

		if (exact_VW_time_vector[jj] < 0.001) {
			num_of_J1_roots = 5;
		}
		else if ((exact_VW_time_vector[jj] >= 0.001) && (exact_VW_time_vector[jj] < 0.01)) {
			num_of_J1_roots = 3;
		}
		else {
			num_of_J1_roots = 1;
		}

		GlobalStiffMatrix stiff_matrix_calculator(Num_Pavt_Layers, num_prony_elements, Num_VE_Layer);
		Complex_1D_Array* global_load_vector_initializer = new Complex_1D_Array(2 * Num_Pavt_Layers);
		Complex_1D_Array* complexdisp_vector_initializer = new Complex_1D_Array(2 * Num_Pavt_Layers);

		for (m = 0; m < num_of_J1_roots; m++) {
			if (m < 1) {
				num_of_subregions = 10;
			}
			else if (m >= 1 && m < 3) {
				num_of_subregions = 5;
			}
			else {
				num_of_subregions = 3;
			}

			for (mm = 1; mm < num_of_subregions + 1; mm++) {
				// subregions for each range
				a = bessel_J1_root[m] + (mm - 1) * (bessel_J1_root[m + 1] - bessel_J1_root[m]) / num_of_subregions;
				b = bessel_J1_root[m] + mm * (bessel_J1_root[m + 1] - bessel_J1_root[m]) / num_of_subregions;

				for (mmm = 0; mmm < 6; mmm++) {
					//6 global_load_vectoroint Gaussian Quadrature
					gaussian_point = (b - a) / 2 * gaussian_nodes[mmm] + (b + a) / 2;

					global_load_vector_initializer->Reset();
					global_load_vector_initializer->array_1d[1].x = (Load_Pressure * Load_Radius * alglib::besselj1(gaussian_point * Load_Radius) / gaussian_point * dt);
					global_load_vector_initializer->array_1d[1].y = 0;

					fixed_talbot_r = 2 * prec_dec_digit / (5 * exact_VW_time_vector[jj]);

					// Cutoff time
					if (fixed_talbot_r < gaussian_point)
					{
						continue;
					}

					fixed_talbot_s.x = fixed_talbot_r;
					fixed_talbot_s.y = 0;

					stiff_matrix_calculator.CalculateMatrix(pavement_structure, final_prony_values, fixed_talbot_s, gaussian_point);

					complexdisp_vector_initializer->array_1d = global_load_vector_initializer->array_1d;

					info = 0;
					alglib::cmatrixsolvefast(stiff_matrix_calculator.matrix->array_2d, 2 * Num_Pavt_Layers, complexdisp_vector_initializer->array_1d, info);
					if (info > 0) {
						Defl_fixed_talbot_s = complexdisp_vector_initializer->array_1d[1];
					}
					else {
						continue;
					}

					fixed_talbot_sumf = 0;
					ftemp.x = 0;
					ftemp.y = 0;

					for (nn = 1; nn < prec_dec_digit; nn++) {
						fixed_talbot_theta = (nn * alglib::pi() / prec_dec_digit);
						fixed_talbot_s = fixed_talbot_r * fixed_talbot_theta * ((1 / tan(fixed_talbot_theta)) + imaginary_i);

						stiff_matrix_calculator.CalculateMatrix(pavement_structure, final_prony_values, fixed_talbot_s, gaussian_point);

						complexdisp_vector_initializer->array_1d = global_load_vector_initializer->array_1d;

						info = 0;
						alglib::cmatrixsolvefast(stiff_matrix_calculator.matrix->array_2d, 2 * Num_Pavt_Layers, complexdisp_vector_initializer->array_1d, info);
						if (info > 0) {
							Defl_fixed_talbot_theta = complexdisp_vector_initializer->array_1d[1];
						}
						else {
							continue;
						}

						fixed_talbot_gamma = fixed_talbot_theta + (fixed_talbot_theta * (1 / tan(fixed_talbot_theta)) - 1) * (1 / tan(fixed_talbot_theta));
						ftemp = ViscoWaveEngine::Compute_Exp(exact_VW_time_vector[jj] * fixed_talbot_r * fixed_talbot_theta * ((1 / tan(fixed_talbot_theta)) + imaginary_i)) * Defl_fixed_talbot_theta * (1 + imaginary_i * fixed_talbot_gamma);
						fixed_talbot_sumf = fixed_talbot_sumf + ftemp.x;
					}

					y = fixed_talbot_r / prec_dec_digit * (0.5 * Defl_fixed_talbot_s.x * exp(fixed_talbot_r * exact_VW_time_vector[jj]) + fixed_talbot_sumf);

					for (jk = 0; jk < Num_Sensors; jk++) {
						exact_VW_defl_matrix(jj, jk) = exact_VW_defl_matrix(jj, jk) + ((b - a) / 2 * gaussian_weight[mmm] * gaussian_point * y * alglib::besselj0(gaussian_point * defl_sensor_location[jk]));
					}
				}
			}
		}
	}
#pragma omp barrier

//Spline Interpolation
#pragma omp parallel for private(i, j, spline_bound_l, spline_bound_r, mysplinestruc) firstprivate(exact_VW_defl_vector, temporary_displ) shared(Num_Sensors, exact_VW_defl_matrix, exact_VW_time_vector, total_num_exact_time, Num_Time, output_disp_matrix, output_time) //MP 1.31.18
	for (i = 0; i < Num_Sensors; i++) {
		spline_bound_l = (exact_VW_defl_matrix(1, i) - exact_VW_defl_matrix(0, i)) / (exact_VW_time_vector[1] - exact_VW_time_vector[0]);
		spline_bound_r = (exact_VW_defl_matrix(total_num_exact_time - 1, i) - exact_VW_defl_matrix(total_num_exact_time - 2, i)) / (exact_VW_time_vector[total_num_exact_time - 1] - exact_VW_time_vector[total_num_exact_time - 2]);

		for (j = 0; j < total_num_exact_time; j++) {
			exact_VW_defl_vector[j] = exact_VW_defl_matrix(j, i);
		}

		alglib::spline1dconvcubic(exact_VW_time_vector, exact_VW_defl_vector, total_num_exact_time, 1, spline_bound_l, 1, spline_bound_r, output_time, Num_Time, temporary_displ);
		for (j = 0; j < Num_Time; j++) {
			output_disp_matrix(j, i) = temporary_displ[j];
		}
	}
#pragma omp barrier
//Multi Processing

//Convolution
	if (total_num_exact_time == 29) {
		// For FWD loading only
		for (i = 0; i < Num_Sensors; i++) {
			
			//Initialize temp file
			for (j = 0; j < Num_Time; j++) {
				temporary_displ[j] = 0;
			}

			//Initialize temp1 file
			for (j = 0; j < (2 * Num_Time - 1); j++) {
				temporary_displ1[j] = 0;
			}

			for (j = 0; j < Num_Time; j++) {
				temporary_displ[j] = output_disp_matrix(j, i);
			}

			alglib::convr1d(temporary_displ, Num_Time, fwd_load_history, Num_Time, temporary_displ1);

			for (j = 0; j < Num_Time; j++) {
				output_disp_matrix(j, i) = temporary_displ1[j];
			}
		}
	}

//Export Result in vector form
#pragma omp parallel for private(i, j) shared(Num_Time, Num_Sensors, output_disp_matrix, displacement) //MP 12.16.2016

	for (i = 0; i < Num_Sensors; i++) {
		for (j = 0; j < Num_Time; j++) {
			displacement[i * Num_Time + j] = output_disp_matrix(j, i);
		}
	}
#pragma omp barrier //Multi Processing

	return 0;
}