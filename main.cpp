// https://github.com/raketanamarse/super_pc/

#include <iostream>
#include <cmath>
#include <fstream> 
#include <cstring>
#include "l.hpp"
#include "f.hpp"
#include "vars.hpp"
#include "coords.hpp"
#include <omp.h>
#include <ctime> // заголовочный файл с прототипом функции clock()
using namespace std;

#define file_name "../out.txt"



void write_to_file(double *massiv){
    try{
        ofstream MyFile(file_name);
        int *xyz = new int[3];
        // Write to the file
        MyFile << NX << ' ' << NY << ' ' << NZ <<"\n"; // NX NY NZ
            
        for(int n = 0; n < N; n++){
            
            Count_XYZ(n, xyz);
            MyFile << massiv[n] << ' ' << xyz[0] << ' ' << xyz[1] << ' ' << xyz[2] <<"\n"; // U x y z
        }
        
        // Close the file
        MyFile.close();
    }
    catch(int cod_error){cout << "error write to " << file_name << "with error code" << cod_error;}
}



int main(int, char**) {
    
    cout << "start prog\n";
    unsigned int start_time =  clock(); // начальное время
    

    double *cube = new double[N];//текущий куб
    double *next_cube = new double[N]; //следующий куб

    //инициализация нулями cube БЫСТРОЕ
    memset(cube, '0', sizeof(double)*N); 
    memset(next_cube, '0', sizeof(double)*N);

    // #pragma omp parallel for
    for (unsigned t = 1; t < NT; ++t){

        double z = HZ;
        for (unsigned k = NY_NX; k < N - NY_NX; k += NY_NX){

            double y = HY; 
            for (unsigned j = NX; j < NY_NX - NX; j += NX){

               double x = HX;
               for (unsigned i = 1; i < NX - 1; i += 1){

                    unsigned n = i + j + k;
                    next_cube[n] = cube[n] + HT * (F(x, y, z) + DX * LX(cube[n - 1], cube[n], cube[n + 1]) + DY * LY(cube[n - NX], cube[n], cube[n + NX]) + DZ * LZ(cube[n - NY_NX], cube[n], cube[n + NY_NX]));
                    x += HX;
                }
                y += HY;
            }
            z += HZ;
        }
        auto p = cube;
        cube = next_cube;
        next_cube = p;
    }

    unsigned int end_time = clock(); // конечное время
    unsigned int search_time = end_time - start_time; // искомое время
    // cout << "runtime = " << search_time << " mks" << endl;
    // cout << "runtime = " << search_time/1000 << " ms" << endl;
    cout << "runtime = " << search_time/1000000 << " s" << endl;
    // cout << "runtime = " << (double)search_time/60000000 << " min" << endl;
    cout << "NX = " << NX << endl;
    cout << "NY = " << NY << endl;
    cout << "NZ = " << NZ << endl;
    cout << "NT = " << NT << endl;
    cout << "Number of elements = " << N * NT << endl;

    write_to_file(cube);
    //cout << HT <<endl;
}


//https://github.com/raketanamarse/super_pc/

// сборка и запуск mpi выполняется через команды 
// mpic++ -O3 main.cpp -o main
// mpiexec -n 1 main
// вместо 1 можно указать любое количество процессоров/потоков

// #include <iostream>
// #include <cmath>
// #include <fstream> 
// #include <cstring>
// #include "l.hpp"
// #include "f.hpp"
// #include "vars.hpp"
// #include "mpi.h"
// #include <ctime> // заголовочный файл с прототипом функции clock()
// using namespace std;


// const int TAG = 0;

// int main(int argc, char** argv) {
    
//     int start_time = clock(); // конечное время
//     int size, id;
// 	MPI_Init(&argc, &argv);
// 	MPI_Comm_size(MPI_COMM_WORLD, &size);
// 	MPI_Comm_rank(MPI_COMM_WORLD, &id);

//     const int s = (int)sqrt(size);

//     const int nx = (NX - 2) / s + 2;
//     const int ny = (NY - 2) / s + 2;

//     const int N = nx * ny * NZ; 
//     const int nx_ny = nx * ny;

//     int id_top = id + s;
//     int id_down = id - s;
//     int id_ost_s = id % s;
//     int N_nx_ny = N - nx_ny; // N - nx_ny

//     //инициализация соседей
//     const bool top = (0 <= (id_top) && (id_top) < size);
//     const bool down = (0 <= (id_down) && (id_down) < size);
//     const bool left = (id_ost_s != 0);
//     const bool right = (((id + 1) % s) != 0);

//     int req_size = 0;
//     if (top == true){

//         req_size += 2;
//     }
//     if (down == true){

//         req_size += 2;
//     }
//     if (left == true){

//         req_size += 2;
//     }
//     if (right == true){

//         req_size += 2;
//     }

//     MPI_Request reqs[req_size];
//     MPI_Status stats[req_size];


//     double *cube = new double[N];//текущий куб
//     double *next_cube = new double[N]; //следующий куб

//     //инициализация нулями cube БЫСТРОЕ
//     memset(cube, '0', sizeof(double)*N); 
//     memset(next_cube, '0', sizeof(double)*N);

//     //массивы для обмена данными
//     double* buffer_top;
//     double* buffer_top_r;
//     double* buffer_down;
//     double* buffer_down_r;
//     double* buffer_right;
//     double* buffer_right_r;
//     double* buffer_left;
//     double* buffer_left_r;

//     for (int t = 1; t < NT; ++t){

//         double z = HZ;
//         for (int k = nx_ny; k < N_nx_ny; k += nx_ny){

//             double y = HY * ((ny - 2) * (id / s) + 1); 
//             for (int j = nx; j < nx_ny - nx; j += nx){

//                double x = HX * ((nx - 2) * (id_ost_s) + 1);
//                for (int i = 1; i < nx - 1; i += 1){

//                     int n = i + j + k;
//                     next_cube[n] = cube[n] + HT * (F(x, y, z) + DX * LX(cube[n - 1], cube[n], cube[n + 1]) + DY * LY(cube[n - nx], cube[n], cube[n + nx]) + DZ * LZ(cube[n - nx_ny], cube[n], cube[n + nx_ny]));
//                     x += HX;
//                 }
//                 y += HY;
//             }
//             z += HZ;
//         }
//         auto p = cube;
//         cube = next_cube;
//         next_cube = p;

//         int req_i = 0;



//         //верхний сосед
//         if (top == true){

//             //подготовка сообщения
//             buffer_top = new double[(nx - 2) * (NZ - 2)];

//             int index = 0;

//             for (int k = nx_ny; k < N_nx_ny; k += nx_ny){

//                 for (int i = 1; i < nx - 1; i += 1){

//                     int n = i + (nx_ny - 2 * nx) + k; //j = nx*ny - 2nx
//                     buffer_top[index] = cube[n];
//                     ++index;
//                 }
//             }

//             //отправка сообщения
//             MPI_Isend(buffer_top, (nx - 2) * (NZ - 2), MPI_DOUBLE, id_top, TAG, MPI_COMM_WORLD, &reqs[req_i]);
//             req_i++;


//             buffer_top_r = new double[(nx - 2) * (NZ - 2)];
//             //принятие сообщения
//             MPI_Irecv(buffer_top_r, (nx - 2) * (NZ - 2), MPI_DOUBLE, id_top, TAG, MPI_COMM_WORLD, &reqs[req_i]);
//             req_i++;            
//         }

//         //нижний сосед
//         if (down == true){

//             //подготовка сообщения
//             buffer_down = new double[(nx - 2) * (NZ - 2)];

//             int index = 0;

//             for (int k = nx_ny; k < N_nx_ny; k += nx_ny){

//                 for (int i = 1; i < nx - 1; i += 1){

//                     int n = i + (nx) + k; //j = nx
//                     buffer_down[index] = cube[n];
//                     index++;
//                 }
//             }

//             //отправка сообщения
//             MPI_Isend(buffer_down, (nx - 2) * (NZ - 2), MPI_DOUBLE, id_down, TAG, MPI_COMM_WORLD, &reqs[req_i]);
//             req_i++;

//             buffer_down_r = new double[(nx - 2) * (NZ - 2)];

//             //принятие сообщения
//             MPI_Irecv(buffer_down_r, (nx - 2) * (NZ - 2), MPI_DOUBLE, id_down, TAG, MPI_COMM_WORLD, &reqs[req_i]);
//             req_i++;
//         }

//         //левый сосед
//         if (left == true){

//             //подготовка сообщения
//             buffer_right = new double[(ny - 2) * (NZ - 2)];

//             int index = 0;

//             for (int k = nx_ny; k < N_nx_ny; k += nx_ny){

//                 for (int j = nx; j < nx_ny - nx; j += nx){

//                     int n = 1 + j + k; //i = 1
//                     buffer_right[index] = cube[n];
//                     index++;
//                 }
//             }

//             //отправка сообщения
//             MPI_Isend(buffer_right, (ny - 2) * (NZ - 2), MPI_DOUBLE, id - 1, TAG, MPI_COMM_WORLD, &reqs[req_i]);
//             req_i++;

//             buffer_right_r = new double[(ny - 2) * (NZ - 2)];

//             //принятие сообщения
//             MPI_Irecv(buffer_right_r, (ny - 2) * (NZ - 2), MPI_DOUBLE, id - 1, TAG, MPI_COMM_WORLD, &reqs[req_i]);
//             req_i++;
//         }

//         //правый сосед
//         if (right == true){

//              //подготовка сообщения
//             buffer_left = new double[(ny - 2) * (NZ - 2)];

//             int index = 0;

//             for (int k = nx_ny; k < N_nx_ny; k += nx_ny){

//                 for (int j = nx; j < nx_ny - nx; j += nx){

//                     int n = (nx - 2) + j + k; //i = nx - 2
//                     buffer_left[index] = cube[n];
//                     index++;
//                 }
//             }

//             //отправка сообщения
//             MPI_Isend(buffer_left, (ny - 2) * (NZ - 2), MPI_DOUBLE, id + 1, TAG, MPI_COMM_WORLD, &reqs[req_i]);
//             req_i++;

//             buffer_left_r = new double[(ny - 2) * (NZ - 2)];

//             //принятие сообщения
//             MPI_Irecv(buffer_left_r, (ny - 2) * (NZ - 2), MPI_DOUBLE, id + 1, TAG, MPI_COMM_WORLD, &reqs[req_i]);
//             req_i++;
//         }


//         //ожидание окончания работы всех прцоессоров над текущим кубом
//         MPI_Waitall(req_size, reqs, stats);

//         //обработка принятых сообщений
//         if (top == true){

//             int index = 0;

//             for (int k = nx_ny; k < N_nx_ny; k += nx_ny){

//                 for (int i = 1; i < nx - 1; i += 1){

//                     int n = i + (nx_ny - nx) + k; //j = nx*ny - nx
//                     cube[n] = buffer_top_r[index];
//                     index++;
//                 }
//             }
//         }

//         if (down == true){

//             int index = 0;

//             for (int k = nx_ny; k < N_nx_ny; k += nx_ny){

//                 for (int i = 1; i < nx - 1; i += 1){

//                     int n = i + 0 + k; //j = 0
//                     cube[n] = buffer_down_r[index];
//                     index++;
//                 }
//             }
//         }

//         if (left == true){

//             int index = 0;

//             for (int k = nx_ny; k < N_nx_ny; k += nx_ny){

//                 for (int j = nx; j < nx_ny - nx; j += nx){

//                     int n = j + k; //i = 0
//                     cube[n] = buffer_right_r[index];
//                     index++;
//                 }
//             }
//         }

//         if (right == true){

//             int index = 0;

//             for (int k = nx_ny; k < N_nx_ny; k += nx_ny){

//                 for (int j = nx; j < nx_ny - nx; j += nx){

//                     int n = (nx - 1) + j + k; //i = nx - 1
//                     cube[n] = buffer_left_r[index];
//                     index++;
//                 }
//             }
//         }

//     }

//     //рассчет времени работы программы на кадом из процессоров
//     int end_time = clock(); // конечное время
//     int search_time = end_time - start_time; // искомое время

//     MPI_Status status;

//     //поиск максимального времени исполнения
//     if(id == 0){ // принять данные по времени

//         int time = search_time;
//         int buffer_time;

//         for (int i = 1; i < size; i++){

//             MPI_Recv(&buffer_time, 1, MPI_INT, i, TAG, MPI_COMM_WORLD, &status);
//             if (buffer_time > time){

//                 time = buffer_time;
//             }         
//         }
//         cout << "runtime = " << time << " mks" << endl;
//         cout << "runtime = " << time/1000 << " ms" << endl;
//         cout << "runtime = " << time/1000000 << " s" << endl;
//         cout << "runtime = " << (double)time/60000000 << " min" << endl;


//     }
//     else{ // отправить даные по времени

//         MPI_Send(&search_time, 1, MPI_INT, 0, TAG, MPI_COMM_WORLD);
//     }

    


//     //вывод результатов в файлы (каждый процессор выводит результат в собственный файл, далее при необходимости файлы объеденяются python-скриптом)
//     string file_name =  "./out_result/out_";
//     file_name += to_string(id);
//     file_name += ".txt";
//     //cout << file_name;
//     try{
//         ofstream MyFile(file_name);
//         // Write to the file
//         MyFile << NX << ' ' << NY << ' ' << NZ <<"\n"; // NX NY NZ

//         for (int k = nx_ny; k < N_nx_ny; k += nx_ny){

//             for (int j = nx; j < nx_ny - nx; j += nx){

//                for (int i = 1; i < nx - 1; i += 1){

//                     int n = i + j + k;
//                     MyFile << cube[n] << ' ' << (id_ost_s) * (nx - 2) + i << ' ' << (id / s) * (ny - 2) + j / nx << ' ' << k / nx_ny <<"\n"; // U i j k
                    
//                 }
//             }
//         }
//         // Close the file
//         MyFile.close();
//     }
//     catch(int cod_error){cout << "error write to " << file_name << "with error code" << cod_error;}

    
//     MPI_Finalize();
// }


