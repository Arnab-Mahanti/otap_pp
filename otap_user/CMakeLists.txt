cmake_minimum_required(VERSION 3.16)
project(otap_user)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include_directories(../libotap/src)
include_directories(../libotap/vendor)
include_directories(src)

# add_subdirectory(vendor/matplotplusplus)

add_executable(otap_user src/otap_user.cpp)
target_link_libraries(otap_user PRIVATE libotap)
# target_link_libraries(otap_user PUBLIC matplot)