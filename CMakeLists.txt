cmake_minimum_required(VERSION 3.10)

# 设置项目名称
project(matrix_calculator)

# 设置C语言标准
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

# 添加包含目录
include_directories(${PROJECT_SOURCE_DIR}/inc)

# 添加源文件
set(SOURCES 
    src/main.c
    src/algebra.c
)

# 生成可执行文件
add_executable(matrix_calculator ${SOURCES})