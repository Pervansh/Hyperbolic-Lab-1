#include <iostream>
#include <memory>

#include "AdvectonEq1DSolvers.h"

int main() {
    std::cout << "Hello wrld!" << std::endl;
    
    double a = 1.;
    auto linFlux { [&a](double u) {
        return a * u;
    } };

    std::unique_ptr<double[]> u0 = std::make_unique<double[]>(10);
    
    TaskData<double, decltype(linFlux)> taskData(0., 1., 0.1, 0.1, 1., linFlux, u0.get(), 10);
    
    
    uniformMinmodRecRkMethod<double, SSPRK3>(
        taskData,
        rusanovFlux<double>,
        std::cout
    );

    ///*
    uniformMinmodRecRkMethod<double, HeunsMethodRK>(
        taskData,
        rusanovFlux<double>,
        std::cout
    );

    UniformSimpleLaxFriedrichsFlux<double> laxFriedrichsFlux(taskData.dx, taskData.dt);
    uniformMinmodRecRkMethod<double, SSPRK3>(
        taskData,
        laxFriedrichsFlux,
        std::cout
    );

    uniformMinmodRecRkMethod<double, HeunsMethodRK>(
        taskData,
        laxFriedrichsFlux,
        std::cout
    );
    //*/
    
    return 0;
}