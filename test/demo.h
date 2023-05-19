//
// Created by sjf on 8/17/2022.
//

#ifndef GRAPHINDEX_TEST_H
#define GRAPHINDEX_TEST_H

#include "graph_algorithms.h"
#include "memory_analysis.h"
#include "Comparisons/teseo-master/include/teseo.hpp"

void Concurrent_rw_test();

bool ParallelRead(BVGraphTreeParallel::TopBlock* root);

bool ParallelInsert(BVGraphTreeParallel::TopBlock* root);

bool ParallelDelete(BVGraphTreeParallel::TopBlock* root);

bool BlockCompressionTest(std::string file_path);

bool TeseoTest(std::string file_path);

#endif //GRAPHINDEX_TEST_H
