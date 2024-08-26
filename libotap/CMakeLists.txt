cmake_minimum_required(VERSION 3.16)
project(libotap)

file(GLOB_RECURSE OTAP_SRC CONFIGURE_DEPENDS "${PROJECT_SOURCE_DIR}/src/*.h" "${PROJECT_SOURCE_DIR}/src/*.cpp")

add_subdirectory(vendor/rapidyaml ryml)

add_library(libotap STATIC ${OTAP_SRC})
target_link_libraries(libotap ryml::ryml)
target_include_directories(libotap PRIVATE src)
target_include_directories(libotap PRIVATE vendor)
message(${OTAP_SRC})