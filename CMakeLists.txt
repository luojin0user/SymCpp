cmake_minimum_required(VERSION 3.10)
project(SymCpp LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 设置可执行文件输出目录
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

find_package(symengine CONFIG)

# 设置源文件 
set(SOURCES     
    main.cpp ) 
# 添加可执行文件 
add_executable(SymCpp ${SOURCES})

target_link_libraries(SymCpp PRIVATE symengine gmp)
