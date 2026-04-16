cmake_minimum_required(VERSION 3.20)
project(Gaus_Parallel C)

set(CMAKE_C_STANDARD 11)

# Поиск MPI
find_package(MPI REQUIRED)

# Создаем исполняемый файл
add_executable(Gaus_Parallel
        main.c
        gauss_2d.c
        gauss_1d.c
        test.c
        gauss_2d_parallel.c
        gauss_1d_parallel.c
        matrix_io.c
        main_mpi.c
)

target_link_libraries(Gaus_Parallel PRIVATE MPI::MPI_C)

target_compile_options(Gaus_Parallel PRIVATE
        # Для Release режима
        $<$<CONFIG:Release>:-O3>
        $<$<CONFIG:RelWithDebInfo>:-O2 -g>
        # Для Debug режима (не рекомендуется для замеров!)
        $<$<CONFIG:Debug>:-O0>
)

# Если используете MS-MPI на Windows, добавьте пути к DLL
if(WIN32)
    set_target_properties(Gaus_Parallel PROPERTIES
            VS_DEBUGGER_WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
    )
endif()