#include <algorithm>

Info << nl << "Recalculating synthetic body force...\n" << endl;

// find time-index in runTime that best approximates that of pidTime_vector  through interpolation. Returns upper index of the pidTime_vector
pidTimeInterpIndexHigher =
	std::upper_bound( pidTime_vector.begin(), pidTime_vector.end(), runTime.value() ) - pidTime_vector.begin();


// Lower index of the time vector pidTime_vector is simply the upper-index - 1:
pidTimeInterpIndexLower = pidTimeInterpIndexHigher - 1;


// Reference pressure, P0 (or Pref) is obtained through the interpolation of the pressure values at pidP0(higher index) and pidP0(lower index)
pidTempMeasPref = pidP0[pidTimeInterpIndexLower] + ( pidP0[pidTimeInterpIndexHigher] - pidP0[pidTimeInterpIndexLower] )
	* ( runTime.value()                          - pidTime_vector[pidTimeInterpIndexLower] )
	/ ( pidTime_vector[pidTimeInterpIndexHigher] - pidTime_vector[pidTimeInterpIndexLower] );


for (int i = 1; i < (pidNum.value()); ++i) {
	// Pressure at the i-th nudhing point is obtained through the interpolation of the pressure values at pidP(higher time-index) and pidP(lower time-index)
	pidTempMeasP = pidP_matrix[i-1][pidTimeInterpIndexLower] + ( pidP_matrix[i-1][pidTimeInterpIndexHigher] - pidP_matrix[i-1][pidTimeInterpIndexLower] )
		* ( runTime.value()                                   - pidTime_vector[pidTimeInterpIndexLower] )
		/ ( pidTime_vector[pidTimeInterpIndexHigher]            - pidTime_vector[pidTimeInterpIndexLower] );


	// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	// Error at the i-th nudging point: 
	// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	err = ( ( pidTempMeasP - pidTempMeasPref )/rho.value() - (p[pidPositions_inMesh[i]] - p[pidPositions_inMesh[0]]) );

	// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	// Derivative error:
	// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    err_deriv = ( err - err_prev[i] ) / pidDT;

	// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	// Integral error:
	// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    err_int[i] = err_int[i] + ( err * pidDT );

	// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	// Source term at the i-th nudging point is PID controller:
	// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	pidEp[pidPositions_inMesh[i]] = pidKp.value() * ( err + err_int[i]/pidTi.value() + err_deriv*pidTd.value() );

	// Display info to user:
	Info << "Assimilation point #" << i << endl;
	Info << "Measured   (P - Pref)/rho = " << pidTempMeasP/rho.value() << " m2/s2" << endl;
	Info << "Simulation (P - Pref)/rho = " << p[pidPositions_inMesh[i]] << " m2/s2" << endl;
	Info << "Synthetic source term     = " << pidEp[pidPositions_inMesh[i]] << " m3/s2\n" << endl;

	// Update erro at previous time-step
	err_prev[i] = err;

}
