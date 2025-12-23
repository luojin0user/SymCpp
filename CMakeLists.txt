cmake_minimum_required(VERSION 3.10)
project(SymCpp LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 设置可执行文件输出目录
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_CXX_FLAGS "-O0")


# find_package(symengine CONFIG)
find_package(OpenMP REQUIRED)
find_package(Eigen3 REQUIRED)
# 添加可执行文件 
add_executable(SymCpp main.cpp)

# 设置源文件 
file(GLOB_RECURSE SRC_FILES
    ${PROJECT_SOURCE_DIR}/src/*.cpp
    ${PROJECT_SOURCE_DIR}/src/Boundarys/*.cpp
    ${PROJECT_SOURCE_DIR}/src/Coil/*.cpp
    ${PROJECT_SOURCE_DIR}/src/Regions/*.cpp
    ${PROJECT_SOURCE_DIR}/main.cpp
    # debug.cpp
)
target_sources(SymCpp PRIVATE ${SRC_FILES})


# 头文件
target_include_directories(SymCpp PRIVATE
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/include/Media
    ${PROJECT_SOURCE_DIR}/include/Regions
    ${PROJECT_SOURCE_DIR}/include/Coil
    ${PROJECT_SOURCE_DIR}/include/Boundarys
)

target_link_libraries(SymCpp PRIVATE symengine gmp OpenMP::OpenMP_CXX Eigen3::Eigen)


target_compile_options(SymCpp PRIVATE
    $<$<CONFIG:Release>:-O3>
    $<$<CONFIG:Release>:-march=native>
    $<$<CONFIG:Release>:-ffast-math>
)

