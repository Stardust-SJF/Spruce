cmake_minimum_required(VERSION 3.16)
project(GraphIndex)

FIND_PACKAGE( OpenMP REQUIRED)
if(OPENMP_FOUND)
    message("OPENMP FOUND")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()
set(CMAKE_CXX_STANDARD 20)

#file(GLOB GAPBS "third_party/gapbs/src/*")
add_library(spruce STATIC src/index_algorithms.cpp src/memory_analysis.cpp src/succinct_algorithms.cpp src/spruce_transaction.h src/spruce_transaction.cpp)

add_executable(Spruce main.cpp
        src/memory_analysis.h src/memory_analysis.cpp
        src/header.h
        src/index_algorithms.cpp src/index_algorithms.h
        src/succinct_algorithms.cpp src/succinct_algorithms.h
        test/demo.cpp test/demo.h
#        graph_algorithms.cpp graph_algorithms.h
#        third_party/gapbs/src/bfs.h
#        third_party/gapbs/src/bfs.cc
#        third_party/gapbs/src/platform_atomics.cpp third_party/gapbs/src/util.cpp third_party/gapbs/src/benchmark.cpp third_party/gapbs/src/pr_spmv.h
#        third_party/gapbs/src/benchmark.cpp third_party/gapbs/src/pr_spmv.cc third_party/gapbs/src/tc.h
#        third_party/gapbs/src/tc.cc
#        third_party/gapbs/src/sssp.h
#        third_party/gapbs/src/sssp.cc
#        third_party/gapbs/src/cc_sv.h
#        third_party/gapbs/src/cc_sv.cc
#        third_party/gapbs/src/platform_atomics.cpp block_optimization.cpp block_optimization.h
        src/spruce_transaction.h src/spruce_transaction.cpp)
target_link_libraries(Spruce
#        divsufsort
#        divsufsort64
#        numa
#        event
#        event_pthreads
        pthread
        junction
        turf
#        /usr/local/lib/libtbb.so.12.7
        )
