// ViscoWave.h

// Returns ViscoWave Deflections
extern __declspec(dllexport) int ViscoWave(double* displacement, double* Sigmoid, double* Pavement,
					double Load_Pressure, double Load_Radius, double* Sensor_Location, double* Time,
					double* Timehistory, double dt, int num_prony_elements, int Num_Pavt_Layers, int Num_Sensors, 
					int Num_Time, int Num_VE_Layer);