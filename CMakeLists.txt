cmake_minimum_required(VERSION 3.5)
PROJECT(LANGEVIN)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++11 -g -O3")
ADD_SUBDIRECTORY(src bin)
