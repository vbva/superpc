#pragma once
#include "vars.hpp"

//функции для рассчета координат x y z при заданном n

const unsigned NY_NX = NY * NX; 


void Count_XYZ(unsigned n, int* arr){
    
    arr[2] = n / NY_NX;
    arr[1] = (n - arr[2] * NY_NX) / NX;
    arr[0] = n - arr[2] * NY_NX - arr[1] * NX; 
    return;
}


//функция для рассчета координат n при заданном x y z

unsigned Count_n(unsigned x, unsigned y, unsigned z){
    
    return z * NY_NX + y * NX + x;
}


