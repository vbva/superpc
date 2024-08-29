#pragma once

//Глобальные константы

const unsigned NX = 122, NY = 62, NZ = 22; //Значение количества точек по осям X, Y, Z

const double HX = 1 / double(NX - 1), HY = 1 / double(NY - 1), HZ = 1 / double(NZ - 1); //Значение шага по осям X, Y, Z

const double DX = 0.25, DY = 0.15, DZ = 0.1; //Значение dx, dy, dz в тензоре D

const double PI = 3.14159265359; //число Пи

const double HT = 0.45 / (DX/HX/HX + DY/HY/HY + DZ/HZ/HZ); //Значение шага по T

const unsigned NT = unsigned(1 / HT); //Значение количества шагов по T

const unsigned N = NX * NY * NZ; // общее количество элеметов в массиве (кубе)