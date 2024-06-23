//
// Created by sjf on 3/15/2022.
//

#ifndef GRAPHINDEX_GRAPH_PREPROCESS_H
#define GRAPHINDEX_GRAPH_PREPROCESS_H

#include "header.h"
#include <inttypes.h>

// g500 definition
typedef struct packed_edge {
    uint32_t v0_low;
    uint32_t v1_low;
    uint32_t high; /* v1 in high half, v0 in low half */
} packed_edge;

// Convert graph format to undirected format used for metis etc.
int ConvertToMetisFormat(const std::string &graph_path);

bool ConvertToUndirectedGraph(std::string inputFileName);

// Convert graph format to ligra type (adj list) https://github.com/jshun/ligra
int ConvertToLigraFormat(const std::string &graph_path);

int ConvertToGFEFormat(const std::string &graph_path);

int ConvertG500toEdgeList(char* input_file);

int ConvertEdgeListToG500(char* input_file);

int RemovetopthousandDataset(const std::string &graph_path);

int AnalyzeDataset(const std::string &graph_path);

#endif //GRAPHINDEX_GRAPH_PREPROCESS_H
