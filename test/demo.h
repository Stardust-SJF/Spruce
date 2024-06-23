//
// Created by sjf on 8/17/2022.
//

#ifndef GRAPHINDEX_TEST_H
#define GRAPHINDEX_TEST_H

//#include "graph_algorithms.h"
#include "../src/memory_analysis.h"
#include "../src/index_algorithms.h"
#include "../src/spruce_transaction.h"

void Concurrent_rw_test();

bool ParallelRead(SpruceTransVer &spruce);

bool ParallelInsert(SpruceTransVer &spruce);

bool ParallelDelete(SpruceTransVer &spruce);


#endif //GRAPHINDEX_TEST_H
