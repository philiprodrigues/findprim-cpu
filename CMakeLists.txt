cmake_minimum_required(VERSION 3.6)
project(findprim-cpu CXX)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx2 -pthread")

include_directories($ENV{BOOST_INC})
link_directories($ENV{BOOST_LIB})

add_library(findprim SHARED process_avx2.cpp process_naive.cpp design_fir.cpp)

add_executable(do_processing do_processing.cxx)
set_property(TARGET do_processing PROPERTY CXX_STANDARD 17)
target_link_libraries(do_processing boost_program_options findprim)

add_executable(check_output check_output.cpp)
set_property(TARGET check_output PROPERTY CXX_STANDARD 17)
target_link_libraries(check_output findprim)

add_executable(run_naive run_naive.cpp)
set_property(TARGET run_naive PROPERTY CXX_STANDARD 17)
target_link_libraries(run_naive findprim)

add_executable(print_taps print_taps.cpp)
set_property(TARGET print_taps PROPERTY CXX_STANDARD 17)
target_link_libraries(print_taps findprim)

