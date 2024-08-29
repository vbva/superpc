cmake_minimum_required(VERSION 3.0.0)
project(super_pc VERSION 0.1.0)

# Добавляем флаги компиляции и линковки для MPI
set(CMAKE_CXX_COMPILE_FLAGS "${CMAKE_CXX_COMPILE_FLAGS} ${MPI_CXX_COMPILE_FLAGS}")
set(CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} ${MPI_CXX_LINK_FLAGS}")
# Указываем компилятор MPI (если используется OpenMPI)
set(CMAKE_CXX_COMPILER mpicxx)


include(CTest)
enable_testing()

find_package(MPI REQUIRED)
add_executable(super_pc main.cpp f.hpp l.hpp vars.hpp coords.hpp)



# Включим MPI компилятор и компоновщик
include_directories(SYSTEM ${MPI_INCLUDE_PATH})


# Линковка с библиотекой MPI
target_link_libraries(super_pc ${MPI_C_LIBRARIES})



set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
