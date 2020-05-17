cmake_minimum_required (VERSION 3.5)

project("TMP" CXX)

set(CMAKE_CXX_STANDARD 17)

if("${CMAKE_CXX_COMPUTER_ID}" STREQUAL "MSVS")
	add_compile_options(/W4/D_CRT_SECURE_NO_WARNING)
elseif("${CMAKE_CXX_COMPUTER_ID}" STREQUAL "GNU")
	add_compile_options(-Wall -Wextra -Wpedantie)
endif()

add_subdirectory(libutils)

add_subdirectory(lab-1-encapsulation)

add_subdirectory(lab-2-overload)

add_subdirectory(lab-3-polymorphism)

add_subdirectory(lab-4-inheritance)

add_subdirectory(lab-6-sort_algorythms)

add_subdirectory(lab-8-data_structures)

add_subdirectory(lab-10-graphs)