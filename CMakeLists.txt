cmake_minimum_required(VERSION 3.25)
project(Simple-RayTracer-Program CXX)

set(CMAKE_CXX_STANDARD 17)

add_executable(Simple-RayTracer-Program Raytracer.cpp RayUtilities.h)
