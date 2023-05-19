//
// Created by sjf on 3/22/2022.
//

#ifndef GRAPHINDEX_INDEX_ALGORITHMS_H
#define GRAPHINDEX_INDEX_ALGORITHMS_H

#include "header.h"
#include "memory_analysis.h"
#include "succinct_algorithms.h"

#define UNLOCKED 0
#define READ_LOCKED 2
#define WRITE_LOCKED 1
#define RESTART_THRESHOLD 10
#define BLOCK_SIZE 64

extern std::atomic<long int> space_count_temp;

bool BuildSubgraph(const std::string &subgraph_path, const std::string &graph_path, int nodes_num, int part_num);

class SpruceUniverse{
public:
    typedef struct _top_block{
        uint8_t bitmap_8kb[1<<13];
        std::atomic<uint64_t> ptr_to_children[1<<16];
        std::atomic<uint8_t> mtx[1<<10];
        std::atomic<uint8_t> obsolete_flag;
    } TopBlock;

    typedef struct _middle_block{
        uint8_t bitmap_8kb[1<<13];
        std::atomic<uint64_t> ptr_to_children[1<<10];
        std::atomic<uint8_t> mtx[1<<10];
        std::atomic<uint8_t> obsolete_flag;
    } MiddleBlock;

    typedef struct _weighted_edge {
        uint64_t src;
        uint64_t des;
        double weight;
    } WeightedEdge;

    typedef struct _weighted_out_edge {
        uint64_t des;
        double weight;
    } WeightedOutEdge;

    typedef struct _weighted_edge_4B {
        uint32_t src;
        uint32_t des;
        float weight;
    } WeightedEdge4B;

    typedef struct _weighted_out_edge_4B {
        uint32_t des;
        float weight;
    } WeightedOutEdge4B;


//    //May not need;
//    typedef uint64_t* PtrToPtrBlock;

    // 32B/block

    typedef struct _ptr_block {
        std::atomic<uint8_t> obsolete_flag;
        std::atomic<uint64_t> ptr_to_buffer[64];
        std::atomic<uint8_t> buffer_locks[64];
    } PtrBlock;

    typedef struct _adj_subsequent_block_one {
        std::atomic<uint8_t> type;
        std::atomic<uint8_t> obsolete_flag;
        std::atomic<uint8_t> fb_flag;
        std::atomic<uint32_t> timestamp;
        std::atomic<uint64_t> bitvector_64;
        WeightedOutEdge adj_vertex[4];
    }AdjSubsequentBlockOne;

    typedef struct _adj_subsequent_block_two {
        std::atomic<uint8_t> type;
        std::atomic<uint8_t> obsolete_flag;
        std::atomic<uint8_t> fb_flag;
        std::atomic<uint32_t> timestamp;
        std::atomic<uint64_t> bitvector_64;
        WeightedOutEdge adj_vertex[8];
    }AdjSubsequentBlockTwo;

    typedef struct _adj_subsequent_block_three {
        std::atomic<uint8_t> type;
        std::atomic<uint8_t> obsolete_flag;
        std::atomic<uint8_t> fb_flag;
        std::atomic<uint32_t> timestamp;
        std::atomic<uint64_t> bitvector_64;
        WeightedOutEdge adj_vertex[16];
    }AdjSubsequentBlockThree;

    typedef struct _adj_subsequent_block_four {
        std::atomic<uint8_t> type;
        std::atomic<uint8_t> obsolete_flag;
        std::atomic<uint8_t> fb_flag;
        std::atomic<uint32_t> timestamp;
        std::atomic<uint64_t> bitvector_64;
        WeightedOutEdge adj_vertex[32];
    }AdjSubsequentBlockFour;

    typedef struct _adj_subsequent_block_five {
        std::atomic<uint8_t> type;
        std::atomic<uint8_t> obsolete_flag;
        std::atomic<uint8_t> fb_flag;
        std::atomic<uint32_t> timestamp;
        std::atomic<uint64_t> bitvector_64;
        WeightedOutEdge adj_vertex[64];
        std::atomic<uint64_t> next_block;
    }AdjSubsequentBlockFive;

    typedef struct _adj_subsequent_block_one_4B {
        std::atomic<uint8_t> type;
        std::atomic<uint8_t> obsolete_flag;
        std::atomic<uint8_t> fb_flag;
        std::atomic<uint32_t> timestamp;
        std::atomic<uint64_t> bitvector_64;
        WeightedOutEdge4B adj_vertex[4];
    }AdjSubsequentBlockOne4B;

    typedef struct _adj_subsequent_block_two_4B {
        std::atomic<uint8_t> type;
        std::atomic<uint8_t> obsolete_flag;
        std::atomic<uint8_t> fb_flag;
        std::atomic<uint32_t> timestamp;
        std::atomic<uint64_t> bitvector_64;
        WeightedOutEdge4B adj_vertex[8];
    }AdjSubsequentBlockTwo4B;

    typedef struct _adj_subsequent_block_three_4B {
        std::atomic<uint8_t> type;
        std::atomic<uint8_t> obsolete_flag;
        std::atomic<uint8_t> fb_flag;
        std::atomic<uint32_t> timestamp;
        std::atomic<uint64_t> bitvector_64;
        WeightedOutEdge4B adj_vertex[16];
    }AdjSubsequentBlockThree4B;

    typedef struct _adj_subsequent_block_four_4B {
        std::atomic<uint8_t> type;
        std::atomic<uint8_t> fb_flag;
        std::atomic<uint8_t> obsolete_flag;
        std::atomic<uint32_t> timestamp;
        std::atomic<uint64_t> bitvector_64;
        WeightedOutEdge4B adj_vertex[32];
    }AdjSubsequentBlockFour4B;

    typedef struct _adj_subsequent_block_five_4B {
        std::atomic<uint8_t> type;
        std::atomic<uint8_t> obsolete_flag;
        std::atomic<uint8_t> fb_flag;
        std::atomic<uint32_t> timestamp;
        std::atomic<uint64_t> bitvector_64;
        WeightedOutEdge4B adj_vertex[64];
        std::atomic<uint64_t> next_block;
    }AdjSubsequentBlockFive4B;

    inline static bool CompareOutEdges (WeightedOutEdge a, WeightedOutEdge b) {
        return a.des < b.des;
    }

    inline static bool CompareOutEdges4B (WeightedOutEdge4B a, WeightedOutEdge4B b) {
        return a.des < b.des;
    }

    typedef junction::ConcurrentMap_Linear<turf::u32, TopBlock*> ConcurrentMap;
    std::atomic<uint8_t> fb_flag; //4-B flag
    ConcurrentMap spruce_hash;
    TopBlock* top_block;

    static void ReadGraphToVector(const std::string &graph_path, std::vector<WeightedEdge> &edges, bool undirected_flag = 0, bool weight_flag = 0);

    static void InsertEdgeFromVector(SpruceUniverse &spruce, std::vector<WeightedEdge> &edges, bool shuffle_flag = 0);

    static void ClearStatistics();

    static void PrintStatistics();

    static uint64_t GetDegree(SpruceUniverse &spruce, uint64_t from_node_id);

    inline static TopBlock* CreateTopBlock() {
        auto root = (SpruceUniverse::TopBlock*) malloc(sizeof(SpruceUniverse::TopBlock));
        memset(root, 0, sizeof(SpruceUniverse::TopBlock));
        return root;
    }

    static bool BuildTree(SpruceUniverse &spruce, const std::string &graph_path, bool undirected_flag = 0, bool weight_flag = 0);

    static bool InsertEdge(SpruceUniverse &spruce, WeightedEdge edge);

    static bool get_neighbours(SpruceUniverse &spruce, uint64_t from_node_id,
                               std::vector<WeightedOutEdge> &neighbours /*tbb::concurrent_vector<int> &neighbours*/);

    static bool get_neighbours_exclusively(SpruceUniverse &spruce, uint64_t from_node_id,
            /*tbb::concurrent_vector<int> &neighbours*/std::vector<WeightedOutEdge> &neighbours);

    static bool DeleteEdge(SpruceUniverse &spruce, uint64_t from_node_id, uint64_t to_node_id);

    static bool get_neighbours_sorted(SpruceUniverse &spruce, uint64_t from_node_id, std::vector<WeightedOutEdge> &neighbours);

    static bool UpdateEdge(SpruceUniverse &spruce, WeightedEdge edge);

    SpruceUniverse() {
        fb_flag = 0;
        top_block = CreateTopBlock();
        spruce_hash.assign(0, top_block);
    };

    inline static WeightedOutEdge convert_to_8b (WeightedOutEdge4B e) {
        return {e.des, e.weight};
    }
};

#endif //GRAPHINDEX_INDEX_ALGORITHMS_H
