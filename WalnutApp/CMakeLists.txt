cmake_minimum_required(VERSION 3.16)

project(walnut_app)
set(CMAKE_CXX_STANDARD 17)

include_directories(../libotap/src)
include_directories(../libotap/vendor)
include_directories(../Walnut/src)
include_directories(../Walnut/vendor)


file(GLOB_RECURSE WALNUT_APP_SRC CONFIGURE_DEPENDS "${PROJECT_SOURCE_DIR}/src/*.h" "${PROJECT_SOURCE_DIR}/src/*.cpp" "${PROJECT_SOURCE_DIR}/vendor/*.h" "${PROJECT_SOURCE_DIR}/vendor/*.cpp")

add_executable(walnutapp
${WALNUT_APP_SRC}
)
add_subdirectory(vendor/nfd)

target_link_libraries(walnutapp PRIVATE walnut libotap nfd)