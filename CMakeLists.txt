cmake_minimum_required (VERSION 3.8)

project(lab-10-graphs)

set(CMAKE_CXX_STANDARD 17)

if("${CMAKE_CXX_COMPUTER_ID}" STREQUAL "MSVS")
	add_compile_options(/W4/D_CRT_SECURE_NO_WARNING)
elseif("${CMAKE_CXX_COMPUTER_ID}" STREQUAL "GNU")
	add_compile_options(-Wall -Wextra -Wpedantie)
endif()

add_subdirectory("lab-10-graphs")
