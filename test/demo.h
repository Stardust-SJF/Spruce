//
// Created by sjf on 8/17/2022.
//

#ifndef GRAPHINDEX_TEST_H
#define GRAPHINDEX_TEST_H

//#include "graph_algorithms.h"
#include "../src/memory_analysis.h"
#include "../src/index_algorithms.h"

void Concurrent_rw_test();

bool ParallelRead(SpruceUniverse &spruce);

bool ParallelInsert(SpruceUniverse &spruce);

bool ParallelDelete(SpruceUniverse &spruce);


#endif //GRAPHINDEX_TEST_H
