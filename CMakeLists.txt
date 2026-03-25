cmake_minimum_required(VERSION 3.20)
project(Gaus_Parallel C)

set(CMAKE_C_STANDARD 11)

find_package(MPI REQUIRED)

add_executable(Gaus_Parallel main.c gauss_2d.c gauss_1d.c)

target_link_libraries(Gaus_Parallel PUBLIC MPI::MPI_C)

# Явно указываем путь к заголовкам MPI (если нужно)
include_directories("C:/Program Files (x86)/Microsoft SDKs/MPI/Include")