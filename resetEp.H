#include <algorithm>

Info << nl << "Resetting previous and integral error values for synthetic source term...\n" << endl;

for (int i = 1; i < (pidNum.value()); ++i) {
    err_prev[i] = 0; // error at time step t-1 is zeroed 
    err_int[i] = 0;  // integral error is zeroed
}
