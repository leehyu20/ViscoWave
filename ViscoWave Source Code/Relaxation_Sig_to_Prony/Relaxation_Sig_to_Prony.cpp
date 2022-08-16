// Relaxation_Sig_to_Prony.cpp -- Fit Relaxation Prony Series to Sigmoid Function.
 
#include "Relaxation_Sig_to_Prony.h"
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


int Relaxation_Sig_to_Prony(double* Prony, double* Sigmoid, int Num_Sigmoid)
{
	
	alglib::real_1d_array final_prony_values, sigmoid_coeffs;
	int num_prony_elements;
	int i, j, m, jj;
	alglib::lsfitreport rep;
	alglib::ae_int_t num_points_for_prony_fitting, info;
	alglib::real_1d_array time_for_prony_fitting;
	alglib::real_1d_array relaxation_times, sigmoid_modulus, fitted_prony_series;
	alglib::real_2d_array matrix_for_prony_fit;


	sigmoid_coeffs.setlength(4*Num_Sigmoid);

	for( i = 0; i < 4*Num_Sigmoid; i++){
		sigmoid_coeffs[i]= Sigmoid[i];
	}


	num_points_for_prony_fitting = 671;
	num_prony_elements = 15;
	
	time_for_prony_fitting.setlength(num_points_for_prony_fitting);
	sigmoid_modulus.setlength(num_points_for_prony_fitting);
	fitted_prony_series.setlength(num_prony_elements);
	relaxation_times.setlength(num_prony_elements);
	matrix_for_prony_fit.setlength(num_points_for_prony_fitting, num_prony_elements);
	final_prony_values.setlength((Num_Sigmoid +1)* num_prony_elements);

	for (i=0; i< num_prony_elements; i++){
		if(i==0){
			relaxation_times[i]=0;
		}
		else{
			relaxation_times[i]=pow(10.0, 8-i);
		}
		final_prony_values[(Num_Sigmoid +1)*i+ Num_Sigmoid]= relaxation_times[i];
	}

	for (j=0; j< Num_Sigmoid; j++){
	
		for (i=0; i< num_points_for_prony_fitting; i++){
			time_for_prony_fitting[i]=pow(10,-6.7+i*0.02);
			sigmoid_modulus[i]=pow(10, sigmoid_coeffs[0+4*j]+ sigmoid_coeffs[1+4*j]/(1+exp(sigmoid_coeffs[2+4*j]+ sigmoid_coeffs[3+4*j]*log10(time_for_prony_fitting[i]))));

			for(jj=0; jj< num_prony_elements; jj++){
				if(jj==0){
					matrix_for_prony_fit(i,jj)=1.0;
				}
				else{
					matrix_for_prony_fit(i,jj)=exp(-time_for_prony_fitting[i]/ relaxation_times[jj]);
				}
			}
		}
		alglib::lsfitlinear(sigmoid_modulus, matrix_for_prony_fit, num_points_for_prony_fitting, num_prony_elements, info, fitted_prony_series, rep);

		for (i=0; i< num_prony_elements; i++){
			final_prony_values[(Num_Sigmoid +1)*i+j]= fitted_prony_series[i];
		}
	}

	//Export Result in vector form

	m=0;
	for( i = 0; i < Num_Sigmoid +1; i++){
		for( j = 0; j < num_prony_elements; j++){

			Prony[m]= final_prony_values(m);
			m=m+1;
		}
	}

	return 0;

}

