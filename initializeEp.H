#include <algorithm>

Info << nl << "Initializing error values for synthetic source term...\n" << endl;

double pidTempMeasPref = 0; // initialize the time-interpolated value of the reference pressure, P0 (or Pref)
double pidTempMeasP = 0; // initialize the time-interpolated value of the pressure of i-th nudging point. This varible will be individually updated for every nudging point.

double pidDT = 1;  // time-step for the pid controller

int pidTimeInterpIndexLower = 0; // lower input time-index for the interpolation
int pidTimeInterpIndexHigher = 0; // upper input time-index for the interpolation

double err = 0;         // Error of the PID controller. This varible will be individually updated for every nudging point.
double err_deriv = 0;   // Derivative error of the PID controller. This varible will be individually updated for every nudging point.

// For all nudging points, pidNum, define a "error at previous time-step" and "integral error". Vectors are "number of nudging points", pidNum, long.
std::vector<double> err_prev(pidNum.value());
std::vector<double> err_int(pidNum.value());

// Assign zero as initial solution
for (int i = 1; i < (pidNum.value()); ++i) {
    err_prev[i] = 0; // error at previous time step
    err_int[i] = 0;  // integral error
}
