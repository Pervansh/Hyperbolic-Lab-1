#pragma once

#include <iostream>

#include "AdvectionEq1DSolvers.h"

template <typename T>
struct Advection1dTaskFileData {
    // Left boundry of computation domain
    T a;
    // Right boundry of computation domain
    T b;
    // Spatial step
    T dx;
    // Time step
    T dt;
    // Time of simulation end
    T tEnd;
    // Initial data array (u0[i] = u0_{i}, where i - cell number), i = 0, ..., N - 1
    const T* u0;
    // Number of elements in u0 array (i.e. a number of cells). required: N >= 2!
    int N;

    int rkMethodId;
    int fluxId;
    int interpolationMethodId;
};

void advection1dTaskRead(std::istream& input) {

}
