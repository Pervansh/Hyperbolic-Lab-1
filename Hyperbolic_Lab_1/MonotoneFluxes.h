#pragma once

/*
// Rusanov monotone flux struct
template <typename T>
struct RusanovFlux {
    static inline T flux(T fl, T fr, T ul, T ur) {
        return T(0.5f) * (fl + fr - std::max(std::fabs(ur), std::fabs(ul)) * (ur - ul));
    }
};
*/

// Rusanov monotone flux
template <typename T>
inline T rusanovFlux(T fl, T fr, T ul, T ur) {
    return T(0.5f) * (fl + fr - std::max(std::fabs(ur), std::fabs(ul)) * (ur - ul));
}

// Lax-Friedrichs monotone flux for uniform grid
template <typename T>
class UniformSimpleLaxFriedrichsFlux {
    T _deltaX;
    T _deltaT;
    T _halfS;

public:
    UniformSimpleLaxFriedrichsFlux(T deltaX, T deltaT)
        : _deltaX(deltaX), _deltaT(deltaT), _halfS(T(0.5f)* deltaX / deltaT)
    {}

    inline T operator()(T fl, T fr, T ul, T ur) {
        return T(0.5f) * (fl + fr - _halfS * (ur - ul));
    }
};
