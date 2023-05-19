//
// Created by sjf on 8/17/2022.
//

#include "demo.h"
//#include "third_party/gapbs/src/bfs.h"
#define TEST_NUM 20000
#define AVG_DEGREE 128
#define BARRIER_NUM 1

std::atomic<int> barrier_flag;
std::atomic<long int> space_count_temp = 0;

bool ParallelRead(SpruceUniverse &spruce) {
    // randomly generate edges
    srand((int)time(NULL));
    std::vector<std::pair<uint32_t, uint32_t>> edges;
    uint32_t to_node_id, from_node_id;
    for (int i = 0; i < AVG_DEGREE * TEST_NUM; i++) {
        to_node_id = rand()%TEST_NUM;
        from_node_id = rand()%TEST_NUM;
        edges.push_back(std::make_pair(from_node_id, to_node_id));
    }
    barrier_flag++;
    while (barrier_flag < BARRIER_NUM);
    double start = omp_get_wtime();
#pragma omp parallel for num_threads(10)
    for (int i = 0; i < edges.size(); i++) {
        std::vector<SpruceUniverse::WeightedOutEdge> neighbours;
        SpruceUniverse::get_neighbours(spruce, edges[i].first, neighbours);
    }
    double stop = omp_get_wtime();
    printf("Time consumption for reading %ld vertices:%lf s\n",edges.size(),
           double(stop - start));
    return true;
}

bool ParallelInsert(SpruceUniverse &spruce) {
    srand((int)time(NULL));
    std::vector<std::pair<uint32_t, uint32_t>> edges;
    uint32_t to_node_id, from_node_id;
    for (int i = 0; i < AVG_DEGREE * TEST_NUM; i++) {
        to_node_id = rand()%TEST_NUM;
        from_node_id = rand()%TEST_NUM;
        edges.push_back(std::make_pair(from_node_id, to_node_id));
    }
    barrier_flag++;
    while (barrier_flag < BARRIER_NUM);
    double start = omp_get_wtime();
#pragma omp parallel for num_threads(10)
    for (int i = 0; i < edges.size(); i++) {
        SpruceUniverse::InsertEdge(spruce, {edges[i].first, edges[i].second, 0.01});
    }
    double stop = omp_get_wtime();
    printf("Time consumption for inserting %ld edges:%lf s\n",edges.size(),
           double(stop - start));
    return true;
}

bool ParallelDelete(SpruceUniverse &spruce) {
    srand((int)time(NULL));
    std::vector<std::pair<uint32_t, uint32_t>> edges;
    uint32_t to_node_id, from_node_id;
    for (int i = 0; i < AVG_DEGREE * TEST_NUM; i++) {
        to_node_id = rand()%TEST_NUM;
        from_node_id = rand()%TEST_NUM;
        edges.push_back(std::make_pair(from_node_id, to_node_id));
    }
    barrier_flag++;
    while (barrier_flag < BARRIER_NUM);
    double start = omp_get_wtime();
#pragma omp parallel for num_threads(10)
    for (int i = 0; i < edges.size(); i++) {
        SpruceUniverse::DeleteEdge(spruce, edges[i].first, edges[i].second);
    }
    double stop = omp_get_wtime();
    printf("Time consumption for deleting %ld edges:%lf s\n",edges.size(),
           double(stop - start));
    return true;
}









