//
// Created by sjf on 8/17/2022.
//

#include "demo.h"
//#include "third_party/gapbs/src/bfs.h"
#define TEST_NUM 2000000
#define AVG_DEGREE 128

std::atomic<int> barrier_flag;
std::atomic<long int> space_count_temp = 0;

void Concurrent_rw_test() {
    barrier_flag = 0;
    SpruceUniverse spruce;
    memset(root, 0, sizeof(BVGraphTreeParallel::TopBlock));
    std::thread insert_edges(ParallelInsert, spruce);
    std::thread read_neighbours(ParallelRead, spruce);
    std::thread delete_edges(ParallelDelete, spruce);
    insert_edges.join();
    read_neighbours.join();
    delete_edges.join();
}

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
    while (barrier_flag < 3);
    double start = omp_get_wtime();
#pragma omp parallel for num_threads(10)
    for (int i = 0; i < edges.size(); i++) {
        std::vector<uint64_t> neighbours;
        SpruceUniverse::get_neighbours(root, &edges[i].first, neighbours);
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
    while (barrier_flag < 3);
    double start = omp_get_wtime();
#pragma omp parallel for num_threads(10)
    for (int i = 0; i < edges.size(); i++) {
        SpruceUniverse::InsertEdge(root, &edges[i].first, &edges[i].second);
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
    while (barrier_flag < 3);
    double start = omp_get_wtime();
#pragma omp parallel for num_threads(10)
    for (int i = 0; i < edges.size(); i++) {
        SpruceUniverse::DeleteEdge(root, &edges[i].first, &edges[i].second);
    }
    double stop = omp_get_wtime();
    printf("Time consumption for deleting %ld edges:%lf s\n",edges.size(),
           double(stop - start));
    return true;
}









