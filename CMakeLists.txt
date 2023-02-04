cmake_minimum_required(VERSION 3.16)

project("PinAstar")
message("Start build PinAstar")

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(SRC
    ${CMAKE_SOURCE_DIR}/PinAstar/ABP_AStarDecode.cpp
    ${CMAKE_SOURCE_DIR}/PinAstar/AStarDecode_Adaptive.cpp
    ${CMAKE_SOURCE_DIR}/PinAstar/AStarDecode_DM_III.cpp
    ${CMAKE_SOURCE_DIR}/PinAstar/AStarDecode_repeatition.cpp
    ${CMAKE_SOURCE_DIR}/PinAstar/AStarDecode_RSS_II.cpp
    ${CMAKE_SOURCE_DIR}/PinAstar/AStarDecode_RSS_OSC_II.cpp
    ${CMAKE_SOURCE_DIR}/PinAstar/AstarDecode_section_RSS.cpp
    ${CMAKE_SOURCE_DIR}/PinAstar/AstarDecode_segments.cpp
    ${CMAKE_SOURCE_DIR}/PinAstar/SPA_AStarDecoder.cpp
    ${CMAKE_SOURCE_DIR}/PinAstar/Test.cpp
    ${CMAKE_SOURCE_DIR}/AStarDecode_CBC.cpp
    ${CMAKE_SOURCE_DIR}/AStarDecode_DM_I.cpp
    ${CMAKE_SOURCE_DIR}/AStarDecode_DM_II.cpp
    ${CMAKE_SOURCE_DIR}/AStarDecode_ReduceCOM.cpp
    ${CMAKE_SOURCE_DIR}/AStarDecode_RSS_OSC.cpp
    ${CMAKE_SOURCE_DIR}/AStarDecode_RSS.cpp
    ${CMAKE_SOURCE_DIR}/AStarDecode.cpp
    ${CMAKE_SOURCE_DIR}/AWGN.cpp
    ${CMAKE_SOURCE_DIR}/DefineClass.cpp
    ${CMAKE_SOURCE_DIR}/FileOperation.cpp
    ${CMAKE_SOURCE_DIR}/LDPC_OPER.cpp
    ${CMAKE_SOURCE_DIR}/LinearBlockCodes.cpp
    ${CMAKE_SOURCE_DIR}/MatrixOperation.cpp
    ${CMAKE_SOURCE_DIR}/POLAR_OPER.cpp
    ${CMAKE_SOURCE_DIR}/Source.cpp
    ${CMAKE_SOURCE_DIR}/TURBO_OPER.cpp
)

add_executable(PinAstar ${SRC})

add_compile_definitions(_HAS_STD_BYTE=0)

add_compile_options(-g)