#pragma once
#include "vars.hpp"

//Функция f

const double DSum = DX + DY + DZ; //dx + dy + dz

const double DSum_PI_PI = DSum * PI * PI;

inline double F(double x, double y, double z){
    
    return DSum_PI_PI * sin(PI * x) * sin(PI * y) * sin(PI * z);
}