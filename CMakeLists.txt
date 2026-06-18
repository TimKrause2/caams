cmake_minimum_required(VERSION 3.24)
project(caams_demo)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
add_compile_options(-ggdb)

cmake_policy(SET CMP0072 NEW) # Pefer GLVND over legacy GL libraries
find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
find_package(SDL2 REQUIRED)
find_package(Eigen3 REQUIRED)

add_subdirectory(libcaams)

set(CAAMS_LINK_LIBS Eigen3::Eigen OpenGL::GL GLUT::GLUT SDL2::SDL2 caams)

add_subdirectory(apps/chain)
add_subdirectory(apps/mb_gyro)
#add_subdirectory(apps/multi_sp)
add_subdirectory(apps/n21test)
add_subdirectory(apps/p22test)
add_subdirectory(apps/revp_test)
add_subdirectory(apps/sp_sp_test)
add_subdirectory(apps/tennis_racket)
add_subdirectory(apps/uni_test)
add_subdirectory(apps/vehicle)
add_subdirectory(apps/screw_test)
add_subdirectory(apps/planar_test)
