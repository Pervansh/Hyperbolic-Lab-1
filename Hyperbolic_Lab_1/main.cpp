#include <iostream>
#include <memory>

#include "AdvectonEq1DSolvers.h"

int main() {
    std::cout << "Hello world!" << std::endl;
    
    double a = 1.;
    auto linFlux { [&a](double u) {
        return a * u;
    } };

    std::unique_ptr<double[]> u0 = std::make_unique<double[]>(10);
    
    TaskData<double, decltype(linFlux)> taskData(0., 1., 0.1, 0.1, 1., linFlux, u0.get(), 10);

    float aF = 1.;
    auto linFluxF{ [&aF](float u) {
        return aF * u;
    } };

    std::unique_ptr<float[]> u0F = std::make_unique<float[]>(10);

    TaskData<float, decltype(linFluxF)> taskDataF(0.f, 1.f, 0.1f, 0.1f, 1.f, linFluxF, u0F.get(), 10);
    
    
    uniformMinmodRecRkMethod<float, SSPRK3>(
        taskDataF,
        rusanovFlux<float>,
        std::cout
    );

    /*
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
    */
    
    return 0;
}