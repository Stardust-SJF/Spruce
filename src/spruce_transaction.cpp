//
// Created by Stardust on 2023/6/25.
//

#include "spruce_transaction.h"

extern std::atomic<int> type1, type2, type3, type4, type5, type64, edge_array_num;
extern std::atomic<int> middle_block_num_, dense_ptr_block_num_, sparse_ptr_block_num, second_block_num_, third_block_num_, fourth_block_num_, adj_head_num_, adj_subsequent_block_num_;
extern std::atomic<double> total_insert_time;
std::atomic<uint64_t> global_timestamp;
extern std::atomic<int> compact_4, compact_8, compact_16, compact_32;

void SpruceTransVer::ReadGraphToVector(const std::string &graph_path, std::vector<WeightedEdge> &edges, bool undirected_flag,
                                  bool weight_flag) {
    edges.clear();
    //Open file
    std::ifstream initial_graph(graph_path, std::ios::in);
    if (!initial_graph) {
        std::cout << "File Error!! " << std::endl;
        return;
    }

    std::string temp_line;
    int cnt = 0;
    std::string hyperlink = "Hyperlink";
    std::string twitter_2010 = "twitter-2010";
    std::string graph500 = "graph500";
    if (graph_path.find(twitter_2010) != std::string::npos || graph_path.find("uniform")!= std::string::npos || graph_path.find(hyperlink) != std::string::npos || graph_path.find(graph500) != std::string::npos) {

    } else {
        while (cnt++ < 4) {
            getline(initial_graph, temp_line);
            std::cout << temp_line << std::endl;
        }
    }

    // Read Graph Data (Normal)
    uint32_t from_node_id, to_node_id;

    cnt = 0;
    SpruceTransVer::WeightedEdge buffer_edge;
    srand(time(NULL));
    // randomly generate weights
    if (!undirected_flag) {
        while (initial_graph >> from_node_id) {
            initial_graph >> to_node_id;
            buffer_edge = {from_node_id, to_node_id, (float)rand()/RAND_MAX};
            edges.push_back(buffer_edge);
        }
    }
    else {
        while (initial_graph >> from_node_id) {
            initial_graph >> to_node_id;
            buffer_edge = {from_node_id, to_node_id, (float)rand()/RAND_MAX};
            edges.emplace_back(buffer_edge);
            buffer_edge = {to_node_id, from_node_id, (float)rand()/RAND_MAX};
            edges.emplace_back(buffer_edge);
        }
    }

    initial_graph.close();
}

void SpruceTransVer::InsertEdgeFromVector(SpruceTransVer &spruce, std::vector<WeightedEdge> &edges,
                                          bool shuffle_flag) {
    if (shuffle_flag) {
        std::shuffle(edges.begin(), edges.end(), std::mt19937(std::random_device()()));
    }
    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);
#pragma omp parallel for schedule(guided) num_threads(NUM_THREADS)
    for (int i = 0; i < edges.size(); i++) {
        SpruceTransVer::InsertEdge(spruce, edges[i]);
    }
    clock_gettime(CLOCK_MONOTONIC, &end);
    total_insert_time += (end.tv_sec - start.tv_sec) * 1000000000 + end.tv_nsec - start.tv_nsec;
}

void SpruceTransVer::ClearStatistics() {
    middle_block_num_ = dense_ptr_block_num_ = sparse_ptr_block_num = adj_subsequent_block_num_ = adj_head_num_ = 0;
    type1 = type2 = type3 = type4 = type5 = type64 = edge_array_num = 0;
    total_insert_time = 0;
}

void SpruceTransVer::PrintStatistics() {
    std::cout << "Num of middle blocks: " << middle_block_num_ << std::endl;
    std::cout << "Num of sparse pointer blocks: " << sparse_ptr_block_num << std::endl;
    std::cout << "Num of dense pointer blocks: " << dense_ptr_block_num_ << std::endl;
    std::cout << "Num of adj head blocks: " << adj_head_num_ << std::endl;

    std::cout << "Num of adj type 1 subsequent blocks: " << type1 << std::endl;
    std::cout << "Num of adj type 2 subsequent blocks: " << type2 << std::endl;
    std::cout << "Num of adj type 3 subsequent blocks: " << type3 << std::endl;
    std::cout << "Num of adj type 4 subsequent blocks: " << type4 << std::endl;
    std::cout << "Num of adj type 5 subsequent blocks: " << type5 << std::endl;
    std::cout << "Num of edge arrays: " << edge_array_num << std::endl;
    std::cout << "Num of adj type 64 edge arrays: " << type64 << std::endl;

    std::cout << "Real Insertion time:" << total_insert_time / 1000000000 << "s" << std::endl;

    //actually size may be little smaller; not exactly data structure space
    long int analytical_space =
            sizeof(SpruceTransVer::TopBlock) + sizeof(SpruceTransVer::MiddleBlock) * middle_block_num_ +
            sizeof(uint64_t) * sparse_ptr_block_num * (64 + 1) + sizeof(uint64_t) * dense_ptr_block_num_ * 32 +
            //            sizeof(SpruceTransVer::AdjHeadBlock) * adj_head_num_ +
            sizeof(SpruceTransVer::AdjSubsequentBlockOne) * type1 +
            sizeof(SpruceTransVer::AdjSubsequentBlockTwo) * type2 + sizeof(AdjSubsequentBlockThree) * type3 +
            sizeof(AdjSubsequentBlockFour) * type4 + (sizeof(SpruceTransVer::AdjSubsequentBlockFive)) * type5
            +sizeof(WeightedOutEdge) * (64) * type64 + edge_array_num * 2 * sizeof(uint32_t);
    analytical_space = analytical_space / 1024;
    std::cout << "Analytical space: " << analytical_space << "KB" << std::endl;
}

uint64_t SpruceTransVer::GetDegree(SpruceTransVer &spruce, const uint64_t from_node_id) {
    uint64_t degree = 0;
    //do not use locks;
    SpruceTransVer::TopBlock* root;
    auto hash_index = (uint32_t)(from_node_id >> 32);
    auto from_node_id_low = (uint16_t)from_node_id;
    auto from_node_id_high = (uint16_t)(from_node_id >> 16);
    if (hash_index || spruce.fb_flag) {
        // if existed in hash table
        root = spruce.spruce_hash.get(hash_index);
        // if not existed
        if (!root) {
            root = SpruceTransVer::CreateTopBlock();
            spruce.spruce_hash.assign(hash_index,root);
        }
    }
    else {
        root = spruce.top_block;
    }
    uint64_t temp_ptr;
    if (!get_bit(&(root->bitmap_8kb), (uint32_t) from_node_id_high)) {
        return degree;
    }
    temp_ptr = root->ptr_to_children[from_node_id_high].load();
    auto* middle_block_ptr = (SpruceTransVer::MiddleBlock*)temp_ptr;
    if (!middle_block_ptr) {
        return degree;
    }
    if (!get_bit(&(middle_block_ptr->bitmap_8kb), (uint32_t) from_node_id_low)) {
        return degree;
    }
    auto ptr_block_index = from_node_id_low / 64;
    auto auxiliary_ptr = reinterpret_cast<uint64_t*>(&middle_block_ptr->bitmap_8kb);
    uint64_t auxiliary_64 = auxiliary_ptr[ptr_block_index];
    uint64_t auxiliary_64_rev = __builtin_bswap64(auxiliary_64);
    auto ptr_num = __builtin_popcountl(auxiliary_64);
    auto index_in_64 = (uint32_t) (from_node_id_low % 64);

    SpruceTransVer::AdjSubsequentBlockFive4B* bottom_head_block;


    SpruceTransVer::PtrBlock* ptr_block = (SpruceTransVer::PtrBlock*)middle_block_ptr->ptr_to_children[ptr_block_index].load();
    bottom_head_block = (SpruceTransVer::AdjSubsequentBlockFive4B*)(ptr_block->ptr_to_buffer[index_in_64].load());


    if (!bottom_head_block) {
        return degree;
    }

    SpruceTransVer::AdjSubsequentBlockFive4B* local_head_block_ptr = bottom_head_block;
    uint32_t get_index = 0;
    uint64_t temp_bitvector;

    degree += __builtin_popcountl(~(local_head_block_ptr->bitvector_64));
    auto type = local_head_block_ptr->type.load();
    if(type == 5) {

        if (!(local_head_block_ptr->fb_flag_log_size >> 15)) {
            temp_ptr = local_head_block_ptr->next_block.load();
            if (!temp_ptr) {
                // only happens when block was deleted
                return degree;
            } else {
                auto edge_array = (uint32_t*) temp_ptr;
                uint32_t block_size, delete_num;
                if (!edge_array) {
                    block_size = 0;
                    delete_num = 0;
                } else {
                    block_size = edge_array[0];
                    delete_num = edge_array[1];
                }
                degree = degree + block_size - delete_num;
            }
        }
        else {
            //8B case
            auto local_head_block_ptr_8b = (SpruceTransVer::AdjSubsequentBlockFive*)local_head_block_ptr;
            temp_ptr = local_head_block_ptr->next_block.load();
            if (!temp_ptr) {
                // only happens when block was deleted
                return degree;
            } else {
                auto edge_array = (uint64_t*) temp_ptr;
                uint64_t block_size, delete_num;
                if (!edge_array) {
                    block_size = 0;
                    delete_num = 0;
                } else {
                    block_size = edge_array[0];
                    delete_num = edge_array[1];
                }
                degree = degree + block_size - delete_num;
            }
        }
    }
    return degree;
}

bool SpruceTransVer::BuildTree(SpruceTransVer &spruce, const std::string &graph_path, bool undirected_flag, bool weight_flag) {
    auto root = spruce.top_block;
    memset(root, 0, sizeof(SpruceTransVer::TopBlock));

    SpruceTransVer::ClearStatistics();
    auto g = SpruceTransVer::CreateTopBlock();
    std::vector<WeightedEdge> edges;
    SpruceTransVer::ReadGraphToVector(graph_path, edges, undirected_flag);
    SpruceTransVer::InsertEdgeFromVector(spruce, edges, 1);
    SpruceTransVer::PrintStatistics();
    return true;
}

bool SpruceTransVer::InsertEdge(SpruceTransVer &spruce, SpruceTransVer::WeightedEdge edge) {
    // check the size of graph
    SpruceTransVer::TopBlock * root;
    auto from_node_id = edge.src;
    auto to_node_id = edge.des;
    auto weight = edge.weight;
    auto hash_index = (uint32_t)(from_node_id >> 32);
    auto from_node_id_low = (uint16_t)from_node_id;
    auto from_node_id_high = (uint16_t)(from_node_id >> 16);




    //set up lock status
    int unlocked = UNLOCKED, read_locked = READ_LOCKED, write_locked = WRITE_LOCKED;
    uint8_t unlocked_m = UNLOCKED;
    uint8_t write_locked_m = WRITE_LOCKED;

    //need lock: in case for deletion
    SpruceTransVer::MiddleBlock* middle_block_ptr;
    restart:
    if (hash_index || spruce.fb_flag) {
        // if existed in hash table
        root = spruce.spruce_hash.get(hash_index);
        // if not existed
        if (!root) {
            root = SpruceTransVer::CreateTopBlock();
            spruce.spruce_hash.assign(hash_index,root);
            spruce.fb_flag.store(0x01);
        }
    }
    else {
        root = spruce.top_block;
    }

    restart_middle:
    if (root->obsolete_flag) {
        goto restart;
    }
    if (!get_bit(&(root->bitmap_8kb), (uint32_t) from_node_id_high)) {
        //Noted that when thread get lock, corresponding middle block may be created, thereby recheck is needed;
        //try to get write lock (spin lock)
        while(!root->mtx[from_node_id_high/64].compare_exchange_strong(unlocked_m, write_locked_m)) {
            unlocked_m = UNLOCKED;
        }
        if (!get_bit(&(root->bitmap_8kb), (uint32_t) from_node_id_high)) {
            // Need to malloc new middle block
            middle_block_ptr = (SpruceTransVer::MiddleBlock*) malloc(sizeof(SpruceTransVer::MiddleBlock));
            memset(middle_block_ptr, 0, sizeof(SpruceTransVer::MiddleBlock));
            root->ptr_to_children[from_node_id_high].store((uint64_t)middle_block_ptr);
            set_bit(&(root->bitmap_8kb), (uint32_t) from_node_id_high);
            //Analysis
            middle_block_num_++;
        }
        else {
            //atomically load
            uint64_t temp = root->ptr_to_children[from_node_id_high].load();
            middle_block_ptr = (SpruceTransVer::MiddleBlock*)temp;
        }
        root->mtx[from_node_id_high/64]--/*.store(UNLOCKED)*/;
    }
    else {
        //read atomically
        uint64_t temp = root->ptr_to_children[from_node_id_high].load();
        middle_block_ptr = (SpruceTransVer::MiddleBlock*)temp;
    }

    // Get to middle block
    if (!middle_block_ptr) {
        goto restart_middle;
    }
    restart_bottom:
    if (middle_block_ptr->obsolete_flag) {
        goto restart_middle;
    }
    //now we need to get corresponding uint64_t and check the number of 1
    auto ptr_block_index = from_node_id_low / 64;
    auto auxiliary_ptr = reinterpret_cast<uint64_t*>(&middle_block_ptr->bitmap_8kb);
    uint64_t auxiliary_64 = auxiliary_ptr[ptr_block_index];
    uint64_t auxiliary_64_rev = __builtin_bswap64(auxiliary_64);
//    auto ptr_num = __builtin_popcountl(auxiliary_64);
    auto index_in_64 = (uint32_t)(from_node_id_low % 64);

    //Decide the ptr type due to the number of 1

    SpruceTransVer::AdjSubsequentBlockOne* bottom_head_block;

    int lock_flag = 0;
    unlocked = UNLOCKED;

    if (!get_bit(&auxiliary_64, index_in_64)) {
        // does not exist, lock;
        while(!middle_block_ptr->mtx[ptr_block_index].compare_exchange_strong(unlocked_m, write_locked_m)) {
            unlocked_m = UNLOCKED;
        }
        lock_flag = 1;
    }




    //Get bottom head block
    // change set bit sequence for parallel: firstly edit pointer, then set bit
    // reget values!!!!
    auxiliary_64 = auxiliary_ptr[ptr_block_index];
    auxiliary_64_rev = __builtin_bswap64(auxiliary_64);
    SpruceTransVer::PtrBlock* ptr_block;
    // recheck
    if (!get_bit(&auxiliary_64, index_in_64)) {
        //bottom block does not exist, malloc a new block
        bottom_head_block = (SpruceTransVer::AdjSubsequentBlockOne*) malloc(sizeof(SpruceTransVer::AdjSubsequentBlockOne));
        memset(bottom_head_block, 0, sizeof(SpruceTransVer::AdjSubsequentBlockOne));
        bottom_head_block->bitvector_64 = UINT64_MAX;
        bottom_head_block->type.store(1);

        //Edit middle block bitmap and ptr block
        uint64_t temp = middle_block_ptr->ptr_to_children[ptr_block_index].load();
        ptr_block = (SpruceTransVer::PtrBlock*)temp;
        if(!ptr_block) {
            // + 1 for obsolete flag
            auto new_ptr_block = (SpruceTransVer::PtrBlock*) malloc(sizeof(SpruceTransVer::PtrBlock)) ;
            memset(new_ptr_block, 0, sizeof(SpruceTransVer::PtrBlock));
            new_ptr_block->ptr_to_buffer[index_in_64].store(reinterpret_cast<unsigned long>(bottom_head_block));
            middle_block_ptr->ptr_to_children[ptr_block_index].store((uint64_t)new_ptr_block);
            ptr_block = new_ptr_block;
            sparse_ptr_block_num++;
        }
        else {
            if (ptr_block->obsolete_flag) {
                //obsoleted
                if (lock_flag) {
                    middle_block_ptr->mtx[ptr_block_index]--;
                }
                goto restart_bottom;
            }
            ptr_block->ptr_to_buffer[index_in_64].store(reinterpret_cast<unsigned long>(bottom_head_block));
        }


        set_bit(&(middle_block_ptr->bitmap_8kb), (uint32_t) from_node_id_low);
        //Analysis
//        adj_head_num_++;
        type1++;
    } else {
        //corresponding block exists
        uint64_t temp = middle_block_ptr->ptr_to_children[ptr_block_index].load();
        ptr_block = (SpruceTransVer::PtrBlock*)temp;
//        bottom_head_block = (SpruceTransVer::AdjSubsequentBlockFive*)(ptr_block->ptr_to_buffer[index_in_64].load());
    }



    unlocked = UNLOCKED;
    if (lock_flag) {
        middle_block_ptr->mtx[ptr_block_index]--;
    }


    uint32_t insert_index;
    uint8_t type;

    while(!ptr_block->buffer_locks[index_in_64].compare_exchange_strong(unlocked_m, write_locked_m)) {
        unlocked_m = UNLOCKED;
    }

    auto local_head_block_ptr = (SpruceTransVer::AdjSubsequentBlockFive4B*)ptr_block->ptr_to_buffer[index_in_64].load();
    if (!local_head_block_ptr) {
        ptr_block->buffer_locks[index_in_64]--;
        goto restart_bottom;
    }
    if (local_head_block_ptr->obsolete_flag) {
        ptr_block->buffer_locks[index_in_64]--;
        goto restart_bottom;
    }
    //Edit timestamp
    local_head_block_ptr->timestamp = SpruceTransVer::get_global_timestamp();
    //head insertion

    //Edit for SpruceTransVer

    if ((local_head_block_ptr->fb_flag_log_size>>15)^spruce.fb_flag) {
        // 4B->8B case
        WeightedOutEdge out_edge = {edge.des, edge.weight, local_head_block_ptr->timestamp};
        SpruceTransVer::AdjSubsequentBlockFive* local_head_block_ptr_8b;
        if (local_head_block_ptr->bitvector_64 != 0) {
            //not full yet
            uint64_t temp_bv_rev = __builtin_bswap64(local_head_block_ptr->bitvector_64.load());
            insert_index = __builtin_clzl(temp_bv_rev);
            type = local_head_block_ptr->type;
            clear_bit(&local_head_block_ptr->bitvector_64, insert_index);
            uint64_t temp;
            if (insert_index > ((1 << (type - 1)) * 4 - 1)) {
                //need to malloc new space;
                switch (type) {
                    case 1: {
                        auto new_block = (SpruceTransVer::AdjSubsequentBlockTwo*) malloc(
                                sizeof(AdjSubsequentBlockTwo));
                        memset(new_block, 0, sizeof(AdjSubsequentBlockTwo));
                        temp = (uint64_t)local_head_block_ptr;

                        //copy values
//                        memcpy(new_block,
//                               (SpruceTransVer::AdjSubsequentBlockOne*)(temp),
//                               sizeof(SpruceTransVer::AdjSubsequentBlockOne));
                        new_block->type = local_head_block_ptr->type + 1;
                        new_block->fb_flag_log_size = local_head_block_ptr->fb_flag_log_size.load() & 0x7FFF;
                        new_block->bitvector_64 = local_head_block_ptr->bitvector_64.load();
                        new_block->timestamp = local_head_block_ptr->timestamp.load();
                        for (int c = 0; c < ((1 << (type - 1)) * 4); c++) {
                            new_block->adj_vertex[c] = {local_head_block_ptr->adj_vertex->des, local_head_block_ptr->adj_vertex->weight};
                        }

                        new_block->adj_vertex[insert_index] = out_edge;
                        local_head_block_ptr->obsolete_flag = 1;
                        // add previous block to delete queue
                        ptr_block->ptr_to_buffer[index_in_64].store(reinterpret_cast<uint64_t>(new_block));
                        free((SpruceTransVer::AdjSubsequentBlockOne4B*)(temp));
                        local_head_block_ptr_8b = (SpruceTransVer::AdjSubsequentBlockFive*)new_block;
                        //analysis
                        type2++;
                        type1--;

                        break;
                    }
                    case 2: {
                        auto new_block = (SpruceTransVer::AdjSubsequentBlockThree*) malloc(
                                sizeof(AdjSubsequentBlockThree));
                        memset(new_block, 0, sizeof(AdjSubsequentBlockThree));
                        temp = (uint64_t)local_head_block_ptr;
                        //copy values
                        new_block->type = local_head_block_ptr->type + 1;
                        new_block->fb_flag_log_size = local_head_block_ptr->fb_flag_log_size.load() & 0x7FFF;
                        new_block->bitvector_64 = local_head_block_ptr->bitvector_64.load();
                        new_block->timestamp = local_head_block_ptr->timestamp.load();
                        for (int c = 0; c < ((1 << (type - 1)) * 4); c++) {
                            new_block->adj_vertex[c] = {local_head_block_ptr->adj_vertex->des, local_head_block_ptr->adj_vertex->weight};
                        }
                        new_block->adj_vertex[insert_index] = out_edge;
                        new_block->type.store(3);
                        local_head_block_ptr->obsolete_flag = 1;
                        // add previous block to delete queue
                        ptr_block->ptr_to_buffer[index_in_64].store(reinterpret_cast<uint64_t>(new_block));
                        free((SpruceTransVer::AdjSubsequentBlockTwo*)temp);
                        //analysis
                        type3++;
                        type2--;

                        break;
                    }
                    case 3: {
                        auto new_block = (SpruceTransVer::AdjSubsequentBlockFour*) malloc(
                                sizeof(AdjSubsequentBlockFour));
                        memset(new_block, 0, sizeof(AdjSubsequentBlockFour));
                        temp = (uint64_t)local_head_block_ptr;
                        //copy values
                        new_block->type = local_head_block_ptr->type + 1;
                        new_block->fb_flag_log_size = local_head_block_ptr->fb_flag_log_size.load() & 0x7FFF;
                        new_block->bitvector_64 = local_head_block_ptr->bitvector_64.load();
                        new_block->timestamp = local_head_block_ptr->timestamp.load();
                        for (int c = 0; c < ((1 << (type - 1)) * 4); c++) {
                            new_block->adj_vertex[c] = {local_head_block_ptr->adj_vertex->des, local_head_block_ptr->adj_vertex->weight};
                        }
                        new_block->adj_vertex[insert_index] = out_edge;
                        local_head_block_ptr->obsolete_flag = 1;
                        // add previous block to delete queue
                        ptr_block->ptr_to_buffer[index_in_64].store(reinterpret_cast<uint64_t>(new_block));
                        free((SpruceTransVer::AdjSubsequentBlockThree*)temp);
                        //analysis
                        type4++;
                        type3--;

                        break;

                    }
                    case 4: {
                        auto new_block = (SpruceTransVer::AdjSubsequentBlockFive*) malloc(
                                sizeof(AdjSubsequentBlockFive));
                        temp = (uint64_t)local_head_block_ptr;
                        memset(new_block, 0, sizeof(AdjSubsequentBlockFive));
                        // new_block->bitvector_64 = UINT64_MAX; //no subsequent block
                        //copy values
                        new_block->type = local_head_block_ptr->type + 1;
                        new_block->fb_flag_log_size = local_head_block_ptr->fb_flag_log_size.load() & 0x7FFF;
                        new_block->bitvector_64 = local_head_block_ptr->bitvector_64.load();
                        new_block->timestamp = local_head_block_ptr->timestamp.load();
                        for (int c = 0; c < ((1 << (type - 1)) * 4); c++) {
                            new_block->adj_vertex[c] = {local_head_block_ptr->adj_vertex->des, local_head_block_ptr->adj_vertex->weight};
                        }
                        new_block->adj_vertex[insert_index] = out_edge;
                        local_head_block_ptr->obsolete_flag = 1;
                        // add previous block to delete queue
                        ptr_block->ptr_to_buffer[index_in_64].store(reinterpret_cast<uint64_t>(new_block));
                        free((SpruceTransVer::AdjSubsequentBlockFour*)temp);
                        //analysis
                        type5++;
                        type4--;

                        break;
                    }
                    default: {
                        std::cout << "Impossible type!" << std::endl;
                        break;
                    }

                }
            } else {
                // just change type
                switch (type) {
                    case 1: {
                        auto new_block = (SpruceTransVer::AdjSubsequentBlockOne*) malloc(
                                sizeof(AdjSubsequentBlockOne));
                        memset(new_block, 0, sizeof(AdjSubsequentBlockOne));
                        temp = (uint64_t)local_head_block_ptr;

                        //copy values
//                        memcpy(new_block,
//                               (SpruceTransVer::AdjSubsequentBlockOne*)(temp),
//                               sizeof(SpruceTransVer::AdjSubsequentBlockOne));
                        new_block->type = local_head_block_ptr->type.load();
                        new_block->fb_flag_log_size = local_head_block_ptr->fb_flag_log_size.load() & 0x7FFF;
                        new_block->bitvector_64 = local_head_block_ptr->bitvector_64.load();
                        new_block->timestamp = local_head_block_ptr->timestamp.load();
                        for (int c = 0; c < ((1 << (type - 1)) * 4); c++) {
                            new_block->adj_vertex[c] = {local_head_block_ptr->adj_vertex->des, local_head_block_ptr->adj_vertex->weight};
                        }

                        new_block->adj_vertex[insert_index] = out_edge;
                        local_head_block_ptr->obsolete_flag = 1;
                        // add previous block to delete queue
                        ptr_block->ptr_to_buffer[index_in_64].store(reinterpret_cast<uint64_t>(new_block));
                        free((SpruceTransVer::AdjSubsequentBlockOne4B*)(temp));
                        local_head_block_ptr_8b = (SpruceTransVer::AdjSubsequentBlockFive*)new_block;
                        break;
                    }
                    case 2: {
                        auto new_block = (SpruceTransVer::AdjSubsequentBlockTwo*) malloc(
                                sizeof(AdjSubsequentBlockTwo));
                        memset(new_block, 0, sizeof(AdjSubsequentBlockTwo));
                        temp = (uint64_t)local_head_block_ptr;
                        //copy values
                        new_block->type = local_head_block_ptr->type.load();
                        new_block->fb_flag_log_size = local_head_block_ptr->fb_flag_log_size.load() & 0x7FFF;
                        new_block->bitvector_64 = local_head_block_ptr->bitvector_64.load();
                        new_block->timestamp = local_head_block_ptr->timestamp.load();
                        for (int c = 0; c < ((1 << (type - 1)) * 4); c++) {
                            new_block->adj_vertex[c] = {local_head_block_ptr->adj_vertex->des, local_head_block_ptr->adj_vertex->weight};
                        }
                        new_block->adj_vertex[insert_index] = out_edge;
                        new_block->type.store(3);
                        local_head_block_ptr->obsolete_flag = 1;
                        // add previous block to delete queue
                        ptr_block->ptr_to_buffer[index_in_64].store(reinterpret_cast<uint64_t>(new_block));
                        free((SpruceTransVer::AdjSubsequentBlockTwo*)temp);
                        break;
                    }
                    case 3: {
                        auto new_block = (SpruceTransVer::AdjSubsequentBlockThree*) malloc(
                                sizeof(AdjSubsequentBlockThree));
                        memset(new_block, 0, sizeof(AdjSubsequentBlockThree));
                        temp = (uint64_t)local_head_block_ptr;
                        //copy values
                        new_block->type = local_head_block_ptr->type.load();
                        new_block->fb_flag_log_size = local_head_block_ptr->fb_flag_log_size.load() & 0x7FFF;
                        new_block->bitvector_64 = local_head_block_ptr->bitvector_64.load();
                        new_block->timestamp = local_head_block_ptr->timestamp.load();
                        for (int c = 0; c < ((1 << (type - 1)) * 4); c++) {
                            new_block->adj_vertex[c] = {local_head_block_ptr->adj_vertex->des, local_head_block_ptr->adj_vertex->weight};
                        }
                        new_block->adj_vertex[insert_index] = out_edge;
                        local_head_block_ptr->obsolete_flag = 1;
                        // add previous block to delete queue
                        ptr_block->ptr_to_buffer[index_in_64].store(reinterpret_cast<uint64_t>(new_block));
                        free((SpruceTransVer::AdjSubsequentBlockThree*)temp);
                        //analysis
                        break;

                    }
                    case 4: {
                        auto new_block = (SpruceTransVer::AdjSubsequentBlockFour*) malloc(
                                sizeof(AdjSubsequentBlockFour));
                        temp = (uint64_t)local_head_block_ptr;
                        memset(new_block, 0, sizeof(AdjSubsequentBlockFour));
                        //copy values
                        new_block->type = local_head_block_ptr->type.load();
                        new_block->fb_flag_log_size = local_head_block_ptr->fb_flag_log_size.load() & 0x7FFF;
                        new_block->bitvector_64 = local_head_block_ptr->bitvector_64.load();
                        new_block->timestamp = local_head_block_ptr->timestamp.load();
                        for (int c = 0; c < ((1 << (type - 1)) * 4); c++) {
                            new_block->adj_vertex[c] = {local_head_block_ptr->adj_vertex->des, local_head_block_ptr->adj_vertex->weight};
                        }
                        new_block->adj_vertex[insert_index] = out_edge;
                        local_head_block_ptr->obsolete_flag = 1;
                        // add previous block to delete queue
                        ptr_block->ptr_to_buffer[index_in_64].store(reinterpret_cast<uint64_t>(new_block));
                        free((SpruceTransVer::AdjSubsequentBlockFour*)temp);
                        //analysis
                        break;
                    }
                    case 5: {
                        auto new_block = (SpruceTransVer::AdjSubsequentBlockFive*) malloc(
                                sizeof(AdjSubsequentBlockFive));
                        temp = (uint64_t)local_head_block_ptr;
                        memset(new_block, 0, sizeof(AdjSubsequentBlockFive));
                        //copy values
                        new_block->type = local_head_block_ptr->type.load();
                        new_block->fb_flag_log_size = local_head_block_ptr->fb_flag_log_size.load() & 0x7FFF;
                        new_block->bitvector_64 = local_head_block_ptr->bitvector_64.load();
                        new_block->timestamp = local_head_block_ptr->timestamp.load();
                        for (int c = 0; c < ((1 << (type - 1)) * 4); c++) {
                            new_block->adj_vertex[c] = {local_head_block_ptr->adj_vertex->des, local_head_block_ptr->adj_vertex->weight};
                        }
                        new_block->adj_vertex[insert_index] = out_edge;
                        local_head_block_ptr->obsolete_flag = 1;

                        //also change the edge array from 4B->8B
                        uint32_t old_block_size, old_delete_num;
                        uint64_t temp1 = (uint64_t)local_head_block_ptr;
                        uint64_t temp2 = local_head_block_ptr->next_block.load();
                        auto old_edge_array = (uint32_t*)temp2;
                        if (old_edge_array) {
                            // does not exist, set initial size
                            old_block_size = old_edge_array[0];
                            old_delete_num = old_edge_array[1];
                            auto old_edges = (SpruceTransVer::WeightedOutEdge4B*) (old_edge_array + 2);
                            int64_t old_delete_64 = 0;
                            if (old_delete_num > 64) {
                                old_delete_64 = old_delete_num / 64;
                            }

                            // + 1 for block_size and delete_num (equals to weighted out edge)
                            // resize block according to delete_num
                            auto new_block_size = old_block_size + 1;
                            auto new_edge_array = (uint64_t*) malloc(
                                    sizeof(SpruceTransVer::WeightedOutEdge) * (new_block_size));
                            memset(new_edge_array, 0, sizeof(SpruceTransVer::WeightedOutEdge) * (new_block_size));

                            //then use merge sort to place spaces in new block
                            uint64_t k = 0;
                            uint64_t start1 = 0, start2 = 0;
                            uint64_t end1 = 64, end2 = old_block_size;

                            auto new_edges = (SpruceTransVer::WeightedOutEdge*) (new_edge_array + 2);

                            while (start2 < end2) {
                                if (old_edges[start2].des == UINT32_MAX) {
                                    // skip invalid data
                                    new_edges[k++] = convert_to_8b(old_edges[start2++]);
                                    new_edges[k++].des = UINT64_MAX;
                                    continue;
                                }
                                new_edges[k++] = convert_to_8b(old_edges[start2++]);
                            }

                            //reset subsequent block status
                            new_edge_array[0] = new_block_size -
                                                1; // remember: block_size does not include first 2 uint32_t for information
                            new_edge_array[1] = old_delete_num;  // copy delete_num
                            new_block->next_block.store(reinterpret_cast<uint64_t>(new_edge_array));
                            //do not execute free function in new;
                            // Noted that when execute on large datasets, use free function to avoid memory exceeded
                            // when execute parallel test, comment it to avoid errors
                            free(old_edge_array);
                        }
                        // add previous block to delete queue
                        ptr_block->ptr_to_buffer[index_in_64].store(reinterpret_cast<uint64_t>(new_block));
                        free((SpruceTransVer::AdjSubsequentBlockFour*)temp);
                        break;
                    }
                    default: {
                        std::cout << "Impossible type!" << std::endl;
                        break;
                    }
                }
            }
        } else {
            type = 5;
            // while ends, no spaces
            // Firstly check the block size,we assume that we only check invalid label in deletion
            uint32_t old_block_size, old_delete_num;
            uint64_t temp1 = (uint64_t)local_head_block_ptr;
            uint64_t temp2 = local_head_block_ptr->next_block.load();
            auto old_edge_array = (uint32_t*)temp2;
            if (!old_edge_array) {
                // does not exist, set initial size

                old_block_size = 0;
                old_delete_num = 0;
                // analysis
                edge_array_num++;
            }
            else {
                old_block_size = old_edge_array[0];
                old_delete_num = old_edge_array[1];
            }

            auto old_edges = (SpruceTransVer::WeightedOutEdge4B*)(old_edge_array + 2);

            int64_t old_delete_64 = 0;
            if (old_delete_num > 64) {
                old_delete_64 = old_delete_num/64;
            }

            // + 1 for block_size and delete_num (equals to weighted out edge)
            // resize block according to delete_num
            auto new_block_size = (64 + 1 + old_block_size - old_delete_64 * 64);
            auto new_edge_array = (uint64_t*) malloc(
                    sizeof(SpruceTransVer::WeightedOutEdge) * (new_block_size));
            memset(new_edge_array, 0, sizeof(SpruceTransVer::WeightedOutEdge) * (new_block_size));
//

            //sort the vertex using bubble sort
            for (int i = 0; i < 64; i++) {
                bool flag = false;
                for (int j = 0; j < 64 - i - 1; j++) {
                    if (local_head_block_ptr->adj_vertex[j].des > local_head_block_ptr->adj_vertex[j + 1].des) {
                        auto temp = local_head_block_ptr->adj_vertex[j];
                        local_head_block_ptr->adj_vertex[j] = local_head_block_ptr->adj_vertex[j + 1];
                        local_head_block_ptr->adj_vertex[j + 1] = temp;
                        flag = true;
                    }
                }
                if (!flag) {
                    break;
                }
            }

            //then use merge sort to place spaces in new block
            uint64_t k = 0;
            uint64_t start1 = 0, start2 = 0;
            uint64_t end1 = 64, end2 = old_block_size;

            auto new_edges = (SpruceTransVer::WeightedOutEdge*)(new_edge_array + 2);

            while (start1 < end1 && start2 < end2) {
                if (old_edges[start2].des == UINT32_MAX) {
                    // skip invalid data
                    start2++;
                    continue;
                }
                new_edges[k++] = local_head_block_ptr->adj_vertex[start1].des < old_edges[start2].des ?
                                 convert_to_8b(local_head_block_ptr->adj_vertex[start1++]) : convert_to_8b(old_edges[start2++]);
            }
            while (start1 < end1) {
                new_edges[k++] = convert_to_8b(local_head_block_ptr->adj_vertex[start1++]);
            }
            while (start2 < end2) {
                if (old_edges[start2].des == UINT32_MAX) {
                    // skip invalid data
                    start2++;
                    continue;
                }
                new_edges[k++] = convert_to_8b(old_edges[start2++]);
            }
            while (k < new_block_size - 1) {
                // shift invalid values with size < 64
                new_edges[k++].des = UINT64_MAX;
            }

            auto new_block = (SpruceTransVer::AdjSubsequentBlockFive*) malloc(
                    sizeof(AdjSubsequentBlockFive));
            auto temp = (uint64_t)local_head_block_ptr;
            memset(new_block, 0, sizeof(AdjSubsequentBlockFive));
            //copy values
            new_block->type = local_head_block_ptr->type.load();
            new_block->fb_flag_log_size = local_head_block_ptr->fb_flag_log_size.load() & 0x7FFF;
            new_block->bitvector_64.store(UINT64_MAX);
            new_block->timestamp = local_head_block_ptr->timestamp.load();
            local_head_block_ptr->obsolete_flag = 1;
            // add previous block to delete queue
            ptr_block->ptr_to_buffer[index_in_64].store(reinterpret_cast<uint64_t>(new_block));
            free((SpruceTransVer::AdjSubsequentBlockFour*)temp);

            //reset subsequent block status
            new_edge_array[0] = new_block_size - 1; // remember: block_size does not include first 2 uint32_t for information
            new_edge_array[1] = old_delete_num % 64;  // copy delete_num
            local_head_block_ptr->next_block.store(reinterpret_cast<uint64_t>(new_edge_array));

            // Generally you need to use your own garbage collecting function in your system instead of simply free it.
            if(old_edge_array) free(old_edge_array);

            //Then insert new edge to subsequent block
            new_block->adj_vertex[0] = out_edge;
            clear_bit(&new_block->bitvector_64, 0);

            //analysis
            type64++;

        }
    }
    else if (!spruce.fb_flag) {
        //4B case
        WeightedOutEdge4B out_edge = {(uint32_t)edge.des, (uint32_t)local_head_block_ptr->timestamp, weight};
        if (local_head_block_ptr->bitvector_64 != 0) {
            //not full yet
            uint64_t temp_bv_rev = __builtin_bswap64(local_head_block_ptr->bitvector_64.load());
            insert_index = __builtin_clzl(temp_bv_rev);
            type = local_head_block_ptr->type;
            clear_bit(&local_head_block_ptr->bitvector_64, insert_index);
            uint64_t temp;
            if (insert_index > ((1 << (type - 1)) * 4 - 1)) {
                //need to malloc new space;
                switch (type) {
                    case 1: {
                        auto new_block = (SpruceTransVer::AdjSubsequentBlockTwo4B*) malloc(
                                sizeof(AdjSubsequentBlockTwo4B));
                        memset(new_block, 0, sizeof(AdjSubsequentBlockTwo4B));
                        temp = (uint64_t)local_head_block_ptr;
                        memcpy(new_block,
                               (SpruceTransVer::AdjSubsequentBlockOne4B*)(temp),
                               sizeof(SpruceTransVer::AdjSubsequentBlockOne4B));
                        new_block->adj_vertex[insert_index] = out_edge;
                        new_block->type.store(2);
                        local_head_block_ptr->obsolete_flag = 1;
                        // add previous block to delete queue
                        ptr_block->ptr_to_buffer[index_in_64].store(reinterpret_cast<uint64_t>(new_block));
//                        free((SpruceTransVer::AdjSubsequentBlockOne4B*)(temp));
                        local_head_block_ptr = (SpruceTransVer::AdjSubsequentBlockFive4B*)new_block;
                        //analysis
                        type2++;
                        type1--;

                        break;
                    }
                    case 2: {
                        auto new_block = (SpruceTransVer::AdjSubsequentBlockThree4B*) malloc(
                                sizeof(AdjSubsequentBlockThree4B));
                        memset(new_block, 0, sizeof(AdjSubsequentBlockThree4B));
                        temp = (uint64_t)local_head_block_ptr;
                        memcpy(new_block,
                               (SpruceTransVer::AdjSubsequentBlockTwo4B*)temp,
                               sizeof(SpruceTransVer::AdjSubsequentBlockTwo4B));
                        new_block->adj_vertex[insert_index] = out_edge;
                        new_block->type.store(3);
                        local_head_block_ptr->obsolete_flag = 1;
                        // add previous block to delete queue
                        ptr_block->ptr_to_buffer[index_in_64].store(reinterpret_cast<uint64_t>(new_block));
//                        free((SpruceTransVer::AdjSubsequentBlockTwo4B*)temp);
                        //analysis
                        type3++;
                        type2--;

                        break;
                    }
                    case 3: {
                        auto new_block = (SpruceTransVer::AdjSubsequentBlockFour4B*) malloc(
                                sizeof(AdjSubsequentBlockFour4B));
                        memset(new_block, 0, sizeof(AdjSubsequentBlockFour4B));
                        temp = (uint64_t)local_head_block_ptr;
                        memcpy(new_block,
                               (SpruceTransVer::AdjSubsequentBlockThree4B*)temp,
                               sizeof(SpruceTransVer::AdjSubsequentBlockThree4B));
                        new_block->adj_vertex[insert_index] = out_edge;
                        new_block->type.store(4);
                        local_head_block_ptr->obsolete_flag = 1;
                        // add previous block to delete queue
                        ptr_block->ptr_to_buffer[index_in_64].store(reinterpret_cast<uint64_t>(new_block));
//                        free((SpruceTransVer::AdjSubsequentBlockThree4B*)temp);
                        //analysis
                        type4++;
                        type3--;

                        break;

                    }
                    case 4: {
                        auto new_block = (SpruceTransVer::AdjSubsequentBlockFive4B*) malloc(
                                sizeof(AdjSubsequentBlockFive4B));
                        temp = (uint64_t)local_head_block_ptr;
                        memset(new_block, 0, sizeof(AdjSubsequentBlockFive4B));
                        // new_block->bitvector_64 = UINT64_MAX; //no subsequent block
                        memcpy(new_block,
                               (SpruceTransVer::AdjSubsequentBlockFour4B*)(temp),
                               sizeof(SpruceTransVer::AdjSubsequentBlockFour4B));
                        new_block->adj_vertex[insert_index] = out_edge;
                        new_block->type.store(5);
                        local_head_block_ptr->obsolete_flag = 1;
                        // add previous block to delete queue
                        ptr_block->ptr_to_buffer[index_in_64].store(reinterpret_cast<uint64_t>(new_block));
//                        free((SpruceTransVer::AdjSubsequentBlockFour4B*)temp);
                        //analysis
                        type5++;
                        type4--;

                        break;
                    }
                    default: {
                        std::cout << "Impossible type!" << std::endl;
                        break;
                    }

                }
            } else {
                local_head_block_ptr->adj_vertex[insert_index] = out_edge;
            }
        } else {
            type = 5;
            // while ends, no spaces
            // Firstly check the block size,we assume that we only check invalid label in deletion
            uint32_t old_block_size, old_delete_num;
            uint64_t temp1 = (uint64_t)local_head_block_ptr;
            uint64_t temp2 = local_head_block_ptr->next_block.load();
            auto old_edge_array = (uint32_t*)temp2;
            if (!old_edge_array) {
                // does not exist, set initial size

                old_block_size = 0;
                old_delete_num = 0;
                // analysis
                edge_array_num++;
            }
            else {
                old_block_size = old_edge_array[0];
                old_delete_num = old_edge_array[1];
            }

            auto old_edges = (SpruceTransVer::WeightedOutEdge4B*)(old_edge_array + 2);

            uint32_t old_delete_64 = 0;
            if (old_delete_num > 64) {
                old_delete_64 = old_delete_num/64;
            }

            // + 1 for block_size and delete_num (equals to weighted out edge)
            // resize block according to delete_num
            auto new_block_size = (64 + 1 + old_block_size - old_delete_64 * 64);
            auto new_edge_array = (uint32_t*) malloc(
                    sizeof(SpruceTransVer::WeightedOutEdge4B) * (new_block_size));
            memset(new_edge_array, 0, sizeof(SpruceTransVer::WeightedOutEdge4B) * (new_block_size));
//

            //sort the vertex using bubble sort
            for (int i = 0; i < 64; i++) {
                bool flag = false;
                for (int j = 0; j < 64 - i - 1; j++) {
                    if (local_head_block_ptr->adj_vertex[j].des > local_head_block_ptr->adj_vertex[j + 1].des) {
                        auto temp = local_head_block_ptr->adj_vertex[j];
                        local_head_block_ptr->adj_vertex[j] = local_head_block_ptr->adj_vertex[j + 1];
                        local_head_block_ptr->adj_vertex[j + 1] = temp;
                        flag = true;
                    }
                }
                if (!flag) {
                    break;
                }
            }

            //then use merge sort to place spaces in new block
            uint32_t k = 0;
            uint32_t start1 = 0, start2 = 0;
            uint32_t end1 = 64, end2 = old_block_size;

            auto new_edges = (SpruceTransVer::WeightedOutEdge4B*)(new_edge_array + 2);

            while (start1 < end1 && start2 < end2) {
                if (old_edges[start2].des == UINT32_MAX) {
                    // skip invalid data
                    start2++;
                    continue;
                }
                new_edges[k++] = local_head_block_ptr->adj_vertex[start1].des < old_edges[start2].des ?
                                 local_head_block_ptr->adj_vertex[start1++] : old_edges[start2++];
            }
            while (start1 < end1) {
                new_edges[k++] = local_head_block_ptr->adj_vertex[start1++];
            }
            while (start2 < end2) {
                if (old_edges[start2].des == UINT32_MAX) {
                    // skip invalid data
                    start2++;
                    continue;
                }
                new_edges[k++] = old_edges[start2++];
            }
            while (k < new_block_size - 1) {
                // shift invalid values with size < 64
                new_edges[k++].des = UINT32_MAX;
            }


            //reset head block status
            local_head_block_ptr->bitvector_64.store(UINT64_MAX);

            //reset subsequent block status
            memset(local_head_block_ptr->adj_vertex, 0, sizeof(SpruceTransVer::WeightedOutEdge4B)*64);
            new_edge_array[0] = new_block_size - 1; // remember: block_size does not include first 2 uint32_t for information
            new_edge_array[1] = old_delete_num % 64;  // copy delete_num
            local_head_block_ptr->next_block.store(reinterpret_cast<uint64_t>(new_edge_array));
            // Generally you need to use your own garbage collecting function in your system instead of simply free itrrors
//            if(old_edge_array) free(old_edge_array);

            //Then insert new edge to subsequent block
            local_head_block_ptr->adj_vertex[0] = out_edge;
            clear_bit(&local_head_block_ptr->bitvector_64, 0);

            //analysis
            type64++;

        }
    }
    else {
        //8B case
        WeightedOutEdge out_edge = {edge.des, edge.weight, local_head_block_ptr->timestamp};
        auto local_head_block_ptr_8b = (AdjSubsequentBlockFive*)local_head_block_ptr;
        if (local_head_block_ptr_8b->bitvector_64 != 0) {
            //not full yet
            uint64_t temp_bv_rev = __builtin_bswap64(local_head_block_ptr_8b->bitvector_64.load());
            insert_index = __builtin_clzl(temp_bv_rev);
            type = local_head_block_ptr_8b->type;
            clear_bit(&local_head_block_ptr_8b->bitvector_64, insert_index);
            uint64_t temp;
            if (insert_index > ((1 << (type - 1)) * 4 - 1)) {
                //need to malloc new space;
                switch (type) {
                    case 1: {
                        auto new_block = (SpruceTransVer::AdjSubsequentBlockTwo*) malloc(
                                sizeof(AdjSubsequentBlockTwo));
                        memset(new_block, 0, sizeof(AdjSubsequentBlockTwo));
                        temp = (uint64_t)local_head_block_ptr_8b;
                        memcpy(new_block,
                               (SpruceTransVer::AdjSubsequentBlockOne*)(temp),
                               sizeof(SpruceTransVer::AdjSubsequentBlockOne));
                        new_block->adj_vertex[insert_index] = out_edge;
                        new_block->type.store(2);
                        local_head_block_ptr_8b->obsolete_flag = 1;
                        // add previous block to delete queue
                        ptr_block->ptr_to_buffer[index_in_64].store(reinterpret_cast<uint64_t>(new_block));
                        free((SpruceTransVer::AdjSubsequentBlockOne*)(temp));
                        local_head_block_ptr_8b = (SpruceTransVer::AdjSubsequentBlockFive*)new_block;
                        //analysis
                        type2++;
                        type1--;

                        break;
                    }
                    case 2: {
                        auto new_block = (SpruceTransVer::AdjSubsequentBlockThree*) malloc(
                                sizeof(AdjSubsequentBlockThree));
                        memset(new_block, 0, sizeof(AdjSubsequentBlockThree));
                        temp = (uint64_t)local_head_block_ptr_8b;
                        memcpy(new_block,
                               (SpruceTransVer::AdjSubsequentBlockTwo*)temp,
                               sizeof(SpruceTransVer::AdjSubsequentBlockTwo));
                        new_block->adj_vertex[insert_index] = out_edge;
                        new_block->type.store(3);
                        local_head_block_ptr_8b->obsolete_flag = 1;
                        // add previous block to delete queue
                        ptr_block->ptr_to_buffer[index_in_64].store(reinterpret_cast<uint64_t>(new_block));
                        free((SpruceTransVer::AdjSubsequentBlockTwo*)temp);
                        //analysis
                        type3++;
                        type2--;

                        break;
                    }
                    case 3: {
                        auto new_block = (SpruceTransVer::AdjSubsequentBlockFour*) malloc(
                                sizeof(AdjSubsequentBlockFour));
                        memset(new_block, 0, sizeof(AdjSubsequentBlockFour));
                        temp = (uint64_t)local_head_block_ptr_8b;
                        memcpy(new_block,
                               (SpruceTransVer::AdjSubsequentBlockThree*)temp,
                               sizeof(SpruceTransVer::AdjSubsequentBlockThree));
                        new_block->adj_vertex[insert_index] = out_edge;
                        new_block->type.store(4);
                        local_head_block_ptr_8b->obsolete_flag = 1;
                        // add previous block to delete queue
                        ptr_block->ptr_to_buffer[index_in_64].store(reinterpret_cast<uint64_t>(new_block));
                        free((SpruceTransVer::AdjSubsequentBlockThree*)temp);
                        //analysis
                        type4++;
                        type3--;

                        break;

                    }
                    case 4: {
                        auto new_block = (SpruceTransVer::AdjSubsequentBlockFive*) malloc(
                                sizeof(AdjSubsequentBlockFive));
                        temp = (uint64_t)local_head_block_ptr_8b;
                        memset(new_block, 0, sizeof(AdjSubsequentBlockFive));
                        // new_block->bitvector_64 = UINT64_MAX; //no subsequent block
                        memcpy(new_block,
                               (SpruceTransVer::AdjSubsequentBlockFour*)(temp),
                               sizeof(SpruceTransVer::AdjSubsequentBlockFour));
                        new_block->adj_vertex[insert_index] = out_edge;
                        new_block->type.store(5);
                        local_head_block_ptr_8b->obsolete_flag = 1;
                        // add previous block to delete queue
                        ptr_block->ptr_to_buffer[index_in_64].store(reinterpret_cast<uint64_t>(new_block));
                        free((SpruceTransVer::AdjSubsequentBlockFour*)temp);
                        //analysis
                        type5++;
                        type4--;

                        break;
                    }
                    default: {
                        std::cout << "Impossible type!" << std::endl;
                        break;
                    }

                }
            } else {
                local_head_block_ptr_8b->adj_vertex[insert_index] = out_edge;
            }
        } else {
            type = 5;
            // while ends, no spaces
            // Firstly check the block size,we assume that we only check invalid label in deletion
            uint64_t old_block_size, old_delete_num;
            uint64_t temp1 = (uint64_t)local_head_block_ptr_8b;
            uint64_t temp2 = local_head_block_ptr_8b->next_block.load();
            auto old_edge_array = (uint64_t*)temp2;
            if (!old_edge_array) {
                // does not exist, set initial size

                old_block_size = 0;
                old_delete_num = 0;
                // analysis
                edge_array_num++;
            }
            else {
                old_block_size = old_edge_array[0];
                old_delete_num = old_edge_array[1];
            }

            auto old_edges = (SpruceTransVer::WeightedOutEdge*)(old_edge_array + 2);

            int64_t old_delete_64 = 0;
            if (old_delete_num > 64) {
                old_delete_64 = old_delete_num/64;
            }

            // + 1 for block_size and delete_num (equals to weighted out edge)
            // resize block according to delete_num
            auto new_block_size = (64 + 1 + old_block_size - old_delete_64 * 64);
            auto new_edge_array = (uint64_t*) malloc(
                    sizeof(SpruceTransVer::WeightedOutEdge) * (new_block_size));
            memset(new_edge_array, 0, sizeof(SpruceTransVer::WeightedOutEdge) * (new_block_size));
//

            //sort the vertex using bubble sort
            for (int i = 0; i < 64; i++) {
                bool flag = false;
                for (int j = 0; j < 64 - i - 1; j++) {
                    if (local_head_block_ptr_8b->adj_vertex[j].des > local_head_block_ptr_8b->adj_vertex[j + 1].des) {
                        auto temp = local_head_block_ptr_8b->adj_vertex[j];
                        local_head_block_ptr_8b->adj_vertex[j] = local_head_block_ptr_8b->adj_vertex[j + 1];
                        local_head_block_ptr_8b->adj_vertex[j + 1] = temp;
                        flag = true;
                    }
                }
                if (!flag) {
                    break;
                }
            }

            //then use merge sort to place spaces in new block
            uint64_t k = 0;
            uint64_t start1 = 0, start2 = 0;
            uint64_t end1 = 64, end2 = old_block_size;

            auto new_edges = (SpruceTransVer::WeightedOutEdge*)(new_edge_array + 2);

            while (start1 < end1 && start2 < end2) {
                if (old_edges[start2].des == UINT64_MAX) {
                    // skip invalid data
                    start2++;
                    continue;
                }
                new_edges[k++] = local_head_block_ptr_8b->adj_vertex[start1].des < old_edges[start2].des ?
                                 local_head_block_ptr_8b->adj_vertex[start1++] : old_edges[start2++];
            }
            while (start1 < end1) {
                new_edges[k++] = local_head_block_ptr_8b->adj_vertex[start1++];
            }
            while (start2 < end2) {
                if (old_edges[start2].des == UINT64_MAX) {
                    // skip invalid data
                    start2++;
                    continue;
                }
                new_edges[k++] = old_edges[start2++];
            }
            while (k < new_block_size - 1) {
                // shift invalid values with size < 64
                new_edges[k++].des = UINT64_MAX;
            }


            //reset head block status
            local_head_block_ptr_8b->bitvector_64.store(UINT64_MAX);

            //reset subsequent block status
            memset(local_head_block_ptr_8b->adj_vertex, 0, sizeof(SpruceTransVer::WeightedOutEdge)*64);
            new_edge_array[0] = new_block_size - 1; // remember: block_size does not include first 2 uint32_t for information
            new_edge_array[1] = old_delete_num % 64;  // copy delete_num
            local_head_block_ptr_8b->next_block.store(reinterpret_cast<uint64_t>(new_edge_array));
            // Generally you need to use your own garbage collecting function in your system instead of simply free it
            if(old_edge_array) free(old_edge_array);

            //Then insert new edge to subsequent block
            local_head_block_ptr_8b->adj_vertex[0] = out_edge;
            clear_bit(&local_head_block_ptr_8b->bitvector_64, 0);

            //analysis
            type64++;

        }
    }
    ptr_block->buffer_locks[index_in_64]--;
    return true;
}

bool SpruceTransVer::get_neighbours(SpruceTransVer &spruce, const uint64_t from_node_id,
                                    std::vector<WeightedOutEdgeSimple> &neighbours) {
    //do not use locks
    SpruceTransVer::TopBlock* root;
    uint32_t restart_num = 0;
    auto hash_index = (uint32_t)(from_node_id >> 32);
    auto from_node_id_low = (uint16_t)from_node_id;
    auto from_node_id_high = (uint16_t)(from_node_id >> 16);
    restart:
    if (hash_index || spruce.fb_flag) {
        // if existed in hash table
        root = spruce.spruce_hash.get(hash_index);
        // if not existed
        if (!root) {
            return false;
        }
    }
    else {
        root = spruce.top_block;
    }
    neighbours.clear();
    uint64_t temp_ptr;
    if (!get_bit(&(root->bitmap_8kb), (uint32_t) from_node_id_high)) {
        return false;
    }
    temp_ptr = root->ptr_to_children[from_node_id_high].load();
    auto* middle_block_ptr = (SpruceTransVer::MiddleBlock*)temp_ptr;
    if (!middle_block_ptr) {
        return false;
    }
    if (!get_bit(&(middle_block_ptr->bitmap_8kb), (uint32_t)from_node_id_low)) {
        return false;
    }
    auto ptr_block_index = from_node_id_low / 64;
    auto auxiliary_ptr = reinterpret_cast<uint64_t*>(&middle_block_ptr->bitmap_8kb);
    uint64_t auxiliary_64 = auxiliary_ptr[ptr_block_index];
    uint64_t auxiliary_64_rev = __builtin_bswap64(auxiliary_64);
    auto ptr_num = __builtin_popcountl(auxiliary_64);
    auto index_in_64 = (uint32_t) (from_node_id_low % 64);

    SpruceTransVer::AdjSubsequentBlockFive4B* bottom_head_block;

    SpruceTransVer::PtrBlock* ptr_block = (SpruceTransVer::PtrBlock*)middle_block_ptr->ptr_to_children[ptr_block_index].load();
    bottom_head_block = (SpruceTransVer::AdjSubsequentBlockFive4B*)(ptr_block->ptr_to_buffer[index_in_64].load());

    if (!bottom_head_block) {
        return false;
    }

    SpruceTransVer::AdjSubsequentBlockFive4B* local_head_block_ptr = bottom_head_block;
    uint32_t get_index = 0;
    uint64_t temp_bitvector;

    recheck_bottom_lock:
    while (ptr_block->buffer_locks[index_in_64].load() != UNLOCKED);
    uint32_t timestamp_before = bottom_head_block->timestamp.load();
    // also recheck
    if (ptr_block->buffer_locks[index_in_64].load() == WRITE_LOCKED) {
        goto recheck_bottom_lock;
    }
    // allow reading from obsoleted block
    if (!(local_head_block_ptr->fb_flag_log_size >> 15)){
        //4B case
        if (local_head_block_ptr->bitvector_64 != UINT64_MAX) {
            //value existed
            uint64_t temp_bv_rev = __builtin_bswap64(local_head_block_ptr->bitvector_64.load());
            auto pre_bv = bottom_head_block->bitvector_64.load();
            auto type = local_head_block_ptr->type.load();
            if (type == 5) {
                //type = 5;
                //first check subsequent block
                for (int i = 0; i < 64; i++) {
                    if (!get_bit(&pre_bv, i)) {
                        neighbours.push_back(convert_to_8b_simple(local_head_block_ptr->adj_vertex[i]));
                    }
                }
                //then check edge array
                temp_ptr = local_head_block_ptr->next_block.load();
                auto edge_array = (uint32_t*)temp_ptr;
                uint32_t block_size, delete_space;
                if (!edge_array) {
                    block_size = 0;
                    delete_space = 0;
                }
                else {
                    block_size = edge_array[0];
                    delete_space = edge_array[1];
                }
                uint32_t new_index = neighbours.size(); // use new index as array index to improve efficiency (instead of push_back)
                neighbours.resize(neighbours.size() + block_size - delete_space);
                auto edges_ptr = (SpruceTransVer::WeightedOutEdge4B*)(edge_array + 2);

                for (int i = 0; i < block_size; i++) {
                    if (edges_ptr[i].des != UINT32_MAX) {
                        neighbours[new_index++] = convert_to_8b_simple(edges_ptr[i]);
                    }
                }


            } else {
                //type = 4;
                for (int i = 0; i < (2 << type) ; i++) {
                    if (!get_bit(&pre_bv, i)) {
                        neighbours.push_back(convert_to_8b_simple(local_head_block_ptr->adj_vertex[i]));
                    }
                }
            }
        } else {
            // check edge array
            if (!local_head_block_ptr->next_block.load()) {
                return false;
            }
            else {
                temp_ptr = local_head_block_ptr->next_block.load();
                auto edge_array = (uint32_t*)temp_ptr;
                uint32_t block_size, delete_space;
                block_size = edge_array[0];
                delete_space = edge_array[1];
                uint32_t new_index = neighbours.size(); // use new index as array index to improve efficiency (instead of push_back)
                neighbours.resize(neighbours.size() + block_size - delete_space);
                auto edges_ptr = (SpruceTransVer::WeightedOutEdge4B*)(edge_array + 2);

                for (int i = 0; i < block_size; i++) {
                    if (edges_ptr[i].des != UINT32_MAX) {
                        neighbours[new_index++] = convert_to_8b_simple(edges_ptr[i]);
                    }
                }
            }
        }
    }
    else {
        // 8B case
        auto local_head_block_ptr_8b = (SpruceTransVer::AdjSubsequentBlockFive*)local_head_block_ptr;
        if (local_head_block_ptr_8b->bitvector_64 != UINT64_MAX) {
            //value existed
            uint64_t temp_bv_rev = __builtin_bswap64(local_head_block_ptr_8b->bitvector_64.load());
            auto pre_bv = bottom_head_block->bitvector_64.load();
            auto type = local_head_block_ptr_8b->type.load();
            if (type == 5) {
                //type = 5;
                //first check subsequent block
                for (int i = 0; i < 64; i++) {
                    if (!get_bit(&pre_bv, i)) {
                        neighbours.push_back({local_head_block_ptr_8b->adj_vertex[i].des, local_head_block_ptr_8b->adj_vertex[i].weight});
                    }
                }
                //then check edge array
                temp_ptr = local_head_block_ptr_8b->next_block.load();
                auto edge_array = (uint64_t*)temp_ptr;
                uint64_t block_size, delete_space;
                if (!edge_array) {
                    block_size = 0;
                    delete_space = 0;
                }
                else {
                    block_size = edge_array[0];
                    delete_space = edge_array[1];
                }
                uint64_t new_index = neighbours.size(); // use new index as array index to improve efficiency (instead of push_back)
                neighbours.resize(neighbours.size() + block_size - delete_space);
                auto edges_ptr = (SpruceTransVer::WeightedOutEdge*)(edge_array + 2);

                for (int i = 0; i < block_size; i++) {
                    if (edges_ptr[i].des != UINT64_MAX) {
                        neighbours[new_index++] = {edges_ptr[i].des, edges_ptr[i].weight};
                    }
                }


            } else {
                //type = 4;
                for (int i = 0; i < (2 << type) ; i++) {
                    if (!get_bit(&pre_bv, i)) {
                        neighbours.push_back({local_head_block_ptr_8b->adj_vertex[i].des, local_head_block_ptr_8b->adj_vertex[i].weight});
                    }
                }
            }
        } else {
            // check edge array
            if (!local_head_block_ptr_8b->next_block.load()) {
                return false;
            }
            else {
                temp_ptr = local_head_block_ptr_8b->next_block.load();
                auto edge_array = (uint64_t*)temp_ptr;
                uint64_t block_size, delete_space;
                block_size = edge_array[0];
                delete_space = edge_array[1];
                uint64_t new_index = neighbours.size(); // use new index as array index to improve efficiency (instead of push_back)
                neighbours.resize(neighbours.size() + block_size - delete_space);
                auto edges_ptr = (SpruceTransVer::WeightedOutEdge*)(edge_array + 2);

                for (int i = 0; i < block_size; i++) {
                    if (edges_ptr[i].des != UINT64_MAX) {
                        neighbours[new_index++] = {edges_ptr[i].des, edges_ptr[i].weight};
                    }
                }
            }
        }
    }

    uint32_t timestamp_after = bottom_head_block->timestamp.load();
    if (timestamp_before != timestamp_after) {
        if (restart_num < RESTART_THRESHOLD){
            restart_num++;
            goto restart;
        }
        else {
            // read bottom head block exclusively
            return get_neighbours_exclusively(spruce, from_node_id, neighbours);
        }
    }

    return true;
}

bool SpruceTransVer::get_neighbours_exclusively(SpruceTransVer &spruce, const uint64_t from_node_id,
                                                std::vector<WeightedOutEdgeSimple> &neighbours) {
    //do not use locks
    neighbours.clear();
    SpruceTransVer::TopBlock* root;
    uint32_t restart_num = 0;
    auto hash_index = (uint32_t)(from_node_id >> 32);
    auto from_node_id_low = (uint16_t)from_node_id;
    auto from_node_id_high = (uint16_t)(from_node_id >> 16);
    restart:
    if (hash_index || spruce.fb_flag) {
        // if existed in hash table
        root = spruce.spruce_hash.get(hash_index);
        // if not existed
        if (!root) {
            return false;
        }
    }
    else {
        root = spruce.top_block;
    }

    uint64_t temp_ptr;
    if (!get_bit(&(root->bitmap_8kb), (uint32_t) from_node_id_high)) {
        return false;
    }
    temp_ptr = root->ptr_to_children[from_node_id_high].load();
    auto* middle_block_ptr = (SpruceTransVer::MiddleBlock*)temp_ptr;
    if (!middle_block_ptr) {
        return false;
    }
    if (!get_bit(&(middle_block_ptr->bitmap_8kb), (uint32_t)from_node_id_low)) {
        return false;
    }
    auto ptr_block_index = from_node_id_low / 64;
    auto auxiliary_ptr = reinterpret_cast<uint64_t*>(&middle_block_ptr->bitmap_8kb);
    uint64_t auxiliary_64 = auxiliary_ptr[ptr_block_index];
    uint64_t auxiliary_64_rev = __builtin_bswap64(auxiliary_64);
    auto ptr_num = __builtin_popcountl(auxiliary_64);
    auto index_in_64 = (uint32_t) (from_node_id_low % 64);

    SpruceTransVer::AdjSubsequentBlockFive4B* bottom_head_block;

    SpruceTransVer::PtrBlock* ptr_block = (SpruceTransVer::PtrBlock*)middle_block_ptr->ptr_to_children[ptr_block_index].load();
    bottom_head_block = (SpruceTransVer::AdjSubsequentBlockFive4B*)(ptr_block->ptr_to_buffer[index_in_64].load());

    if (!bottom_head_block) {
        return false;
    }

    SpruceTransVer::AdjSubsequentBlockFive4B* local_head_block_ptr = bottom_head_block;
    uint32_t get_index = 0;
    uint64_t temp_bitvector;
    uint8_t unlocked_m = UNLOCKED;
    uint8_t write_locked_m = WRITE_LOCKED;
    while(!ptr_block->buffer_locks[index_in_64].compare_exchange_strong(unlocked_m, write_locked_m)) {
        unlocked_m = UNLOCKED;
    }

    // allow reading from obsoleted block
    if (!(local_head_block_ptr->fb_flag_log_size >> 15)){
        //4B case
        if (local_head_block_ptr->bitvector_64 != UINT64_MAX) {
            //value existed
            uint64_t temp_bv_rev = __builtin_bswap64(local_head_block_ptr->bitvector_64.load());
            auto pre_bv = bottom_head_block->bitvector_64.load();
            auto type = local_head_block_ptr->type.load();
            if (type == 5) {
                //type = 5;
                //first check subsequent block
                for (int i = 0; i < 64; i++) {
                    if (!get_bit(&pre_bv, i)) {
                        neighbours.push_back(convert_to_8b_simple(local_head_block_ptr->adj_vertex[i]));
                    }
                }
                //then check edge array
                temp_ptr = local_head_block_ptr->next_block.load();
                auto edge_array = (uint32_t*)temp_ptr;
                uint32_t block_size, delete_space;
                if (!edge_array) {
                    block_size = 0;
                    delete_space = 0;
                }
                else {
                    block_size = edge_array[0];
                    delete_space = edge_array[1];
                }
                uint32_t new_index = neighbours.size(); // use new index as array index to improve efficiency (instead of push_back)
                neighbours.resize(neighbours.size() + block_size - delete_space);
                auto edges_ptr = (SpruceTransVer::WeightedOutEdge4B*)(edge_array + 2);

                for (int i = 0; i < block_size; i++) {
                    if (edges_ptr[i].des != UINT32_MAX) {
                        neighbours[new_index++] = convert_to_8b_simple(edges_ptr[i]);
                    }
                }


            } else {
                //type = 4;
                for (int i = 0; i < (2 << type) ; i++) {
                    if (!get_bit(&pre_bv, i)) {
                        neighbours.push_back(convert_to_8b_simple(local_head_block_ptr->adj_vertex[i]));
                    }
                }
            }
        } else {
            // check edge array
            if (!local_head_block_ptr->next_block.load()) {
                return false;
            }
            else {
                temp_ptr = local_head_block_ptr->next_block.load();
                auto edge_array = (uint32_t*)temp_ptr;
                uint32_t block_size, delete_space;
                block_size = edge_array[0];
                delete_space = edge_array[1];
                uint32_t new_index = neighbours.size(); // use new index as array index to improve efficiency (instead of push_back)
                neighbours.resize(neighbours.size() + block_size - delete_space);
                auto edges_ptr = (SpruceTransVer::WeightedOutEdge4B*)(edge_array + 2);

                for (int i = 0; i < block_size; i++) {
                    if (edges_ptr[i].des != UINT32_MAX) {
                        neighbours[new_index++] = convert_to_8b_simple(edges_ptr[i]);
                    }
                }
            }
        }
    }
    else {
        // 8B case
        auto local_head_block_ptr_8b = (SpruceTransVer::AdjSubsequentBlockFive*)local_head_block_ptr;
        if (local_head_block_ptr_8b->bitvector_64 != UINT64_MAX) {
            //value existed
            uint64_t temp_bv_rev = __builtin_bswap64(local_head_block_ptr_8b->bitvector_64.load());
            auto pre_bv = bottom_head_block->bitvector_64.load();
            auto type = local_head_block_ptr_8b->type.load();
            if (type == 5) {
                //type = 5;
                //first check subsequent block
                for (int i = 0; i < 64; i++) {
                    if (!get_bit(&pre_bv, i)) {
                        neighbours.push_back({local_head_block_ptr_8b->adj_vertex[i].des, local_head_block_ptr_8b->adj_vertex[i].weight});
                    }
                }
                //then check edge array
                temp_ptr = local_head_block_ptr_8b->next_block.load();
                auto edge_array = (uint64_t*)temp_ptr;
                uint64_t block_size, delete_space;
                if (!edge_array) {
                    block_size = 0;
                    delete_space = 0;
                }
                else {
                    block_size = edge_array[0];
                    delete_space = edge_array[1];
                }
                uint64_t new_index = neighbours.size(); // use new index as array index to improve efficiency (instead of push_back)
                neighbours.resize(neighbours.size() + block_size - delete_space);
                auto edges_ptr = (SpruceTransVer::WeightedOutEdge*)(edge_array + 2);

                for (int i = 0; i < block_size; i++) {
                    if (edges_ptr[i].des != UINT64_MAX) {
                        neighbours[new_index++] = {edges_ptr[i].des, edges_ptr[i].weight};
                    }
                }


            } else {
                //type = 4;
                for (int i = 0; i < (2 << type) ; i++) {
                    if (!get_bit(&pre_bv, i)) {
                        neighbours.push_back({local_head_block_ptr_8b->adj_vertex[i].des, local_head_block_ptr_8b->adj_vertex[i].weight});
                    }
                }
            }
        } else {
            // check edge array
            if (!local_head_block_ptr_8b->next_block.load()) {
                return false;
            }
            else {
                temp_ptr = local_head_block_ptr_8b->next_block.load();
                auto edge_array = (uint64_t*)temp_ptr;
                uint64_t block_size, delete_space;
                block_size = edge_array[0];
                delete_space = edge_array[1];
                uint64_t new_index = neighbours.size(); // use new index as array index to improve efficiency (instead of push_back)
                neighbours.resize(neighbours.size() + block_size - delete_space);
                auto edges_ptr = (SpruceTransVer::WeightedOutEdge*)(edge_array + 2);

                for (int i = 0; i < block_size; i++) {
                    if (edges_ptr[i].des != UINT64_MAX) {
                        neighbours[new_index++] = {edges_ptr[i].des, edges_ptr[i].weight};
                    }
                }
            }
        }
    }
    ptr_block->buffer_locks[index_in_64]--;
    return true;
}

bool SpruceTransVer::get_neighbours_sorted(SpruceTransVer &spruce, const uint64_t from_node_id,
                                           std::vector<WeightedOutEdgeSimple> &neighbours) {
    //do not use locks
    SpruceTransVer::TopBlock* root;
    uint32_t restart_num = 0;
    auto hash_index = (uint32_t)(from_node_id >> 32);
    auto from_node_id_low = (uint16_t)from_node_id;
    auto from_node_id_high = (uint16_t)(from_node_id >> 16);
    restart:
    if (hash_index || spruce.fb_flag) {
        // if existed in hash table
        root = spruce.spruce_hash.get(hash_index);
        // if not existed
        if (!root) {
            return false;
        }
    }
    else {
        root = spruce.top_block;
    }
    neighbours.clear();
    uint64_t temp_ptr;
    if (!get_bit(&(root->bitmap_8kb), (uint32_t) from_node_id_high)) {
        return false;
    }
    temp_ptr = root->ptr_to_children[from_node_id_high].load();
    auto* middle_block_ptr = (SpruceTransVer::MiddleBlock*)temp_ptr;
    if (!middle_block_ptr) {
        return false;
    }
    if (!get_bit(&(middle_block_ptr->bitmap_8kb), (uint32_t)from_node_id_low)) {
        return false;
    }
    auto ptr_block_index = from_node_id_low / 64;
    auto auxiliary_ptr = reinterpret_cast<uint64_t*>(&middle_block_ptr->bitmap_8kb);
    uint64_t auxiliary_64 = auxiliary_ptr[ptr_block_index];
    uint64_t auxiliary_64_rev = __builtin_bswap64(auxiliary_64);
    auto ptr_num = __builtin_popcountl(auxiliary_64);
    auto index_in_64 = (uint32_t) (from_node_id_low % 64);

    SpruceTransVer::AdjSubsequentBlockFive4B* bottom_head_block;

    SpruceTransVer::PtrBlock* ptr_block = (SpruceTransVer::PtrBlock*)middle_block_ptr->ptr_to_children[ptr_block_index].load();
    bottom_head_block = (SpruceTransVer::AdjSubsequentBlockFive4B*)(ptr_block->ptr_to_buffer[index_in_64].load());

    if (!bottom_head_block) {
        return false;
    }

    SpruceTransVer::AdjSubsequentBlockFive4B* local_head_block_ptr = bottom_head_block;
    uint32_t get_index = 0;
    uint64_t temp_bitvector;

    recheck_bottom_lock:
    while (ptr_block->buffer_locks[index_in_64].load() != UNLOCKED);
    uint32_t timestamp_before = bottom_head_block->timestamp.load();
    // also recheck
    if (ptr_block->buffer_locks[index_in_64].load() == WRITE_LOCKED) {
        goto recheck_bottom_lock;
    }
    // allow reading from obsoleted block
    if (!(local_head_block_ptr->fb_flag_log_size >> 15)){
        //4B case
        if (local_head_block_ptr->bitvector_64 != UINT64_MAX) {
            //value existed
            uint64_t temp_bv_rev = __builtin_bswap64(local_head_block_ptr->bitvector_64.load());
            auto pre_bv = bottom_head_block->bitvector_64.load();
            auto type = local_head_block_ptr->type.load();
            if (type == 5) {
                //type = 5;
                //first check subsequent block
                std::vector<SpruceTransVer::WeightedOutEdgeSimple> temp_neighbours;
                for (int i = 0; i < 64; i++) {
                    if (!get_bit(&pre_bv, i)) {
                        auto new_edge = convert_to_8b(local_head_block_ptr->adj_vertex[i]);
                        temp_neighbours.push_back({new_edge.des, new_edge.weight});
                    }
                }
                sort(temp_neighbours.begin(), temp_neighbours.end(), CompareOutEdgesSimple);
                //then check edge array
                temp_ptr = local_head_block_ptr->next_block.load();
                auto edge_array = (uint32_t*)temp_ptr;
                uint32_t block_size, delete_space;
                if (!edge_array) {
                    block_size = 0;
                    delete_space = 0;
                }
                else {
                    block_size = edge_array[0];
                    delete_space = edge_array[1];
                }
                uint32_t new_index = neighbours.size(); // use new index as array index to improve efficiency (instead of push_back)
                neighbours.resize(temp_neighbours.size() + block_size -delete_space);

                uint32_t start1 = 0, start2 = 0;
                uint32_t end1 = temp_neighbours.size(), end2 = block_size;
                uint32_t  k = 0;
                auto edges_ptr = (SpruceTransVer::WeightedOutEdge4B*)(edge_array + 2);
                while (start1 < end1 && start2 < end2) {
                    if (edges_ptr[start2].des == UINT32_MAX) {
                        // skip invalid data
                        start2++;
                        continue;
                    }
                    neighbours[k++] = temp_neighbours[start1].des < edges_ptr[start2].des ?
                                      temp_neighbours[start1++]: convert_to_8b_simple(edges_ptr[start2++]);
                }
                while (start1 < end1) {
                    neighbours[k++] = {temp_neighbours[start1++].des, temp_neighbours[start1++].weight};
                }
                while (start2 < end2) {
                    if (edges_ptr[start2].des == UINT32_MAX) {
                        // skip invalid data
                        start2++;
                        continue;
                    }
                    neighbours[k++] = convert_to_8b_simple(edges_ptr[start2++]);
                }


            } else {
                //type = 4;
                for (int i = 0; i < (2 << type); i++) {
                    if (!get_bit(&pre_bv, i)) {
                        neighbours.push_back(convert_to_8b_simple(local_head_block_ptr->adj_vertex[i]));
                    }
                }
                sort(neighbours.begin(), neighbours.end(), CompareOutEdgesSimple);
            }
        } else {
            // just check edge array
            temp_ptr = local_head_block_ptr->next_block.load();
            auto edge_array = (uint32_t*)temp_ptr;
            uint32_t block_size, delete_space;
            if (!edge_array) {
                return false;
            }
            else {
                block_size = edge_array[0];
                delete_space = edge_array[1];
            }
            uint32_t new_index = neighbours.size(); // use new index as array index to improve efficiency (instead of push_back)
            neighbours.resize(block_size -delete_space);

            uint32_t start2 = 0;
            uint32_t end2 = block_size;
            uint32_t  k = 0;
            auto edges_ptr = (SpruceTransVer::WeightedOutEdge4B*)(edge_array + 2);
            while (start2 < end2) {
                if (edges_ptr[start2].des == UINT32_MAX) {
                    // skip invalid data
                    start2++;
                    continue;
                }
                neighbours[k++] = convert_to_8b_simple(edges_ptr[start2++]);
            }
        }
    }
    else {
        // 8B case
        auto local_head_block_ptr_8b = (SpruceTransVer::AdjSubsequentBlockFive*)local_head_block_ptr;
        if (local_head_block_ptr_8b->bitvector_64 != UINT64_MAX) {
            //value existed
            uint64_t temp_bv_rev = __builtin_bswap64(local_head_block_ptr_8b->bitvector_64.load());
            auto pre_bv = bottom_head_block->bitvector_64.load();
            auto type = local_head_block_ptr_8b->type.load();
            if (type == 5) {
                //type = 5;
                //first check subsequent block
                std::vector<SpruceTransVer::WeightedOutEdgeSimple> temp_neighbours;
                for (int i = 0; i < 64; i++) {
                    if (!get_bit(&pre_bv, i)) {
                        SpruceTransVer::WeightedOutEdgeSimple new_edge = {local_head_block_ptr_8b->adj_vertex[i].des, local_head_block_ptr_8b->adj_vertex[i].weight};
                        temp_neighbours.push_back(new_edge);
                    }
                }
                sort(temp_neighbours.begin(), temp_neighbours.end(), CompareOutEdgesSimple);
                //then check edge array
                temp_ptr = local_head_block_ptr_8b->next_block.load();
                auto edge_array = (uint64_t*)temp_ptr;
                uint64_t block_size, delete_space;
                if (!edge_array) {
                    block_size = 0;
                    delete_space = 0;
                }
                else {
                    block_size = edge_array[0];
                    delete_space = edge_array[1];
                }
                uint64_t new_index = neighbours.size(); // use new index as array index to improve efficiency (instead of push_back)
                neighbours.resize(temp_neighbours.size() + block_size -delete_space);

                uint64_t start1 = 0, start2 = 0;
                uint64_t end1 = temp_neighbours.size(), end2 = block_size;
                uint64_t  k = 0;
                auto edges_ptr = (SpruceTransVer::WeightedOutEdge*)(edge_array + 2);
                while (start1 < end1 && start2 < end2) {
                    if (edges_ptr[start2].des == UINT64_MAX) {
                        // skip invalid data
                        start2++;
                        continue;
                    }
                    SpruceTransVer::WeightedOutEdgeSimple tt_edge = {edges_ptr[start2].des, edges_ptr[start2].weight};
                    if (temp_neighbours[start1].des < edges_ptr[start2].des){
                        neighbours[k++] = temp_neighbours[start1++];
                    }
                    else {
                        neighbours[k++] = tt_edge;
                        start2++;
                    }
                }
                while (start1 < end1) {
                    neighbours[k++] = temp_neighbours[start1++];
                }
                while (start2 < end2) {
                    if (edges_ptr[start2].des == UINT64_MAX) {
                        // skip invalid data
                        start2++;
                        continue;
                    }
                    SpruceTransVer::WeightedOutEdgeSimple tt_edge = {edges_ptr[start2].des, edges_ptr[start2].weight};
                    neighbours[k++] = tt_edge;
                    start2++;
                }


            } else {
                //type = 4;
                for (int i = 0; i < (2 << type); i++) {
                    if (!get_bit(&pre_bv, i)) {
                        neighbours.push_back({local_head_block_ptr_8b->adj_vertex[i].des, local_head_block_ptr_8b->adj_vertex[i].weight});
                    }
                }
                sort(neighbours.begin(), neighbours.end(), CompareOutEdgesSimple);
            }
        } else {
            // just check edge array
            temp_ptr = local_head_block_ptr_8b->next_block.load();
            auto edge_array = (uint64_t*)temp_ptr;
            uint64_t block_size, delete_space;
            if (!edge_array) {
                return false;
            }
            else {
                block_size = edge_array[0];
                delete_space = edge_array[1];
            }
            uint64_t new_index = neighbours.size(); // use new index as array index to improve efficiency (instead of push_back)
            neighbours.resize(block_size -delete_space);

            uint64_t start2 = 0;
            uint64_t end2 = block_size;
            uint64_t  k = 0;
            auto edges_ptr = (SpruceTransVer::WeightedOutEdge*)(edge_array + 2);
            while (start2 < end2) {
                if (edges_ptr[start2].des == UINT64_MAX) {
                    // skip invalid data
                    start2++;
                    continue;
                }
                neighbours[k++] = {edges_ptr[start2++].des, edges_ptr[start2++].weight};
            }
        }
    }

    uint32_t timestamp_after = bottom_head_block->timestamp.load();
    if (timestamp_before != timestamp_after) {
        goto restart;
    }

    return true;
}

bool
SpruceTransVer::DeleteEdge(SpruceTransVer &spruce, const uint64_t from_node_id, const uint64_t to_node_id) {
    SpruceTransVer::TopBlock * root;
    double temp_weight;
    auto hash_index = (uint32_t)(from_node_id >> 32);
    auto from_node_id_low = (uint16_t)from_node_id;
    auto from_node_id_high = (uint16_t)(from_node_id >> 16);
    restart:
    if (hash_index || spruce.fb_flag) {
        // if existed in hash table
        root = spruce.spruce_hash.get(hash_index);
        // if not existed
        if (!root) {
            root = SpruceTransVer::CreateTopBlock();
            spruce.spruce_hash.assign(hash_index,root);
        }
    }
    else {
        root = spruce.top_block;
    }
    restart_middle:
    if (root->obsolete_flag) {
        goto restart;
    }
    uint64_t temp_ptr;
    if (!get_bit(&(root->bitmap_8kb), (uint32_t) from_node_id_high)) {
        return false;
    }
    temp_ptr = root->ptr_to_children[from_node_id_high].load();
    auto* middle_block_ptr = (SpruceTransVer::MiddleBlock*)temp_ptr;
    if (!middle_block_ptr) {
        return false;
    }
    restart_bottom:
    if (middle_block_ptr->obsolete_flag) {
        goto restart_middle;
    }
    if (!get_bit(&(middle_block_ptr->bitmap_8kb), (uint32_t)from_node_id_low)) {
        return false;
    }

    auto ptr_block_index = from_node_id_low / 64;
    auto auxiliary_ptr = reinterpret_cast<uint64_t*>(&middle_block_ptr->bitmap_8kb);
    uint64_t auxiliary_64 = auxiliary_ptr[ptr_block_index];
    uint64_t auxiliary_64_rev = __builtin_bswap64(auxiliary_64);
    auto ptr_num = __builtin_popcountl(auxiliary_64);
    auto index_in_64 = (uint32_t) (from_node_id_low % 64);


    auto ptr_block = (SpruceTransVer::PtrBlock*)middle_block_ptr->ptr_to_children[ptr_block_index].load();
    auto bottom_head_block = (SpruceTransVer::AdjSubsequentBlockFive4B*)(ptr_block->ptr_to_buffer[index_in_64].load());

    //precheck
    if (!bottom_head_block) {
        return false;
    }
    uint8_t unlocked_m = UNLOCKED, write_locked_m = WRITE_LOCKED;

    uint32_t get_index = 0;
    uint64_t temp_bitvector;
    uint32_t insert_index;
    uint8_t type;


    while(!ptr_block->buffer_locks[index_in_64].compare_exchange_strong(unlocked_m, write_locked_m)) {
        unlocked_m = UNLOCKED;
    }
    //reget
    auto local_head_block_ptr = (SpruceTransVer::AdjSubsequentBlockFive4B*)ptr_block->ptr_to_buffer[index_in_64].load();
    if (!local_head_block_ptr) {
        ptr_block->buffer_locks[index_in_64]--;
        goto restart_bottom;
    }
    if (local_head_block_ptr->obsolete_flag) {
        ptr_block->buffer_locks[index_in_64]--;
        goto restart_bottom;
    }




    auto pre_bv = local_head_block_ptr->bitvector_64.load();
    bool bottom_empty_flag = 0;
    local_head_block_ptr->timestamp += 1;

    if (!(local_head_block_ptr->fb_flag_log_size >> 15)){
        //4B case
        if (local_head_block_ptr->bitvector_64 != UINT64_MAX) {
            //value existed
            uint64_t temp_bv_rev = __builtin_bswap64(local_head_block_ptr->bitvector_64.load());
            auto type = local_head_block_ptr->type.load();
            if (type == 5) {
                //type = 5;
                for (int i = 0; i < 64; i++) {
                    if (!get_bit(&pre_bv, i)) {
                        if (local_head_block_ptr->adj_vertex[i].des == to_node_id) {
//                        local_head_block_ptr->adj_vertex[i].des = UINT32_MAX;    // invalid
                            set_bit(&local_head_block_ptr->bitvector_64, i);
                            if (local_head_block_ptr->bitvector_64.load() == UINT64_MAX) {
                                temp_ptr = local_head_block_ptr->next_block.load();
                                if (!temp_ptr) {
                                    //bottom_head_block empty
                                    bottom_empty_flag = 1;
                                    local_head_block_ptr->obsolete_flag.store(1);
                                }
                                else {
                                    auto next_block_ptr = (uint32_t*)temp_ptr;
                                    if (next_block_ptr[0] == next_block_ptr[1]) {
                                        //delete_num == edge array size
                                        local_head_block_ptr->next_block.store(0);
                                        bottom_empty_flag = 1;
                                        local_head_block_ptr->obsolete_flag.store(1);
                                        // add block_traverse_ptr and edge array ptr to delete queue

                                    }
                                }
                            }
                            temp_weight = local_head_block_ptr -> adj_vertex[i].weight;
                            goto bottom_check;
                        }


                    }
                }

                //check edge array
                temp_ptr = local_head_block_ptr->next_block.load();
                uint32_t* edge_array;
                if (!temp_ptr) {
                    ptr_block->buffer_locks[index_in_64]--;
                    return false;
                }
                else {
                    edge_array = (uint32_t*)temp_ptr;
                }
                auto edge_array_size = edge_array[0];
                auto edges_ptr = (SpruceTransVer::WeightedOutEdge4B*)(edge_array + 2);

                int64_t low = 0;
                int64_t high = edge_array_size - 1;
                int64_t mid;
                while (high >= low) {
                    mid = (low + high) / 2;
                    if ( edges_ptr[mid].des == UINT32_MAX) {
                        // find actual value
                        uint32_t  new_mid = mid;
                        while (new_mid < high) {
                            new_mid++;
                            if (edges_ptr[new_mid].des != UINT32_MAX) {
                                mid = new_mid;
                                break;
                            }
                        }
                    }
                    if(to_node_id < edges_ptr[mid].des){
                        high = mid - 1;
                    }
                    else if (to_node_id == edges_ptr[mid].des) {
                        edges_ptr[mid].des = UINT32_MAX;    // invalid
                        edge_array[1] += 1; //delete_num += 1
//                        block_traverse_ptr -> next_block_size -= 0;
                        if (edge_array[0] == edge_array[1]) {
                            local_head_block_ptr->next_block.store(0);
                            //add edge array pointer to the delete queue
//                        free(edge_array);
                        }
                        temp_weight = edges_ptr[mid].weight;
                        goto bottom_check;
                    }
                    else {
                        low = mid + 1;
                    }
                }
                ptr_block->buffer_locks[index_in_64]--;
                return false;
            } else  {
                //type = 4;
                for (int i = 0; i < (2 << type); i++) {
                    if (!get_bit(&pre_bv, i)) {
                        if (local_head_block_ptr->adj_vertex[i].des == to_node_id) {
//                        local_head_block_ptr->adj_vertex[i].des = UINT32_MAX;    // invalid
                            set_bit(&local_head_block_ptr->bitvector_64, i);
                            if (local_head_block_ptr->bitvector_64.load() == UINT64_MAX) {
                                //bottom_head_block empty
                                bottom_empty_flag = 1;
                                local_head_block_ptr->obsolete_flag.store(1);
                                //add to delete queue
//                            free(block_traverse_ptr);
                            }
                            temp_weight = local_head_block_ptr -> adj_vertex[i].weight;
                            goto bottom_check;
                        }
                    }
                }
                ptr_block->buffer_locks[index_in_64]--;
                return false;
            }
        } else {
            if (local_head_block_ptr->next_block.load()) {
                //only check edge array (subsequent block is empty)
                temp_ptr = local_head_block_ptr->next_block.load();
                uint32_t* edge_array;
                edge_array = (uint32_t*)temp_ptr;

                auto edge_array_size = edge_array[0];
                auto edges_ptr = (SpruceTransVer::WeightedOutEdge4B*)(edge_array + 2);

                int64_t low = 0;
                int64_t high = edge_array_size - 1;
                int64_t mid;
                while (high >= low) {
                    mid = (low + high) / 2;
                    if ( edges_ptr[mid].des == UINT32_MAX) {
                        // find actual value
                        uint32_t  new_mid = mid;
                        while (new_mid < high) {
                            new_mid++;
                            if (edges_ptr[new_mid].des != UINT32_MAX) {
                                mid = new_mid;
                                break;
                            }
                        }
                    }
                    if(to_node_id < edges_ptr[mid].des){
                        high = mid - 1;
                    }
                    else if (to_node_id == edges_ptr[mid].des) {
                        edges_ptr[mid].des = UINT32_MAX;    // invalid
                        edge_array[1] += 1; //delete_num += 1
//                        block_traverse_ptr -> next_block_size -= 0;
                        if (edge_array[0] == edge_array[1]) {
                            local_head_block_ptr->next_block.store(0);
                            //add edge array pointer to the delete queue
//                        free(edge_array);
                        }
                        temp_weight = edges_ptr[mid].weight;
                        goto bottom_check;
                    }
                    else {
                        low = mid + 1;
                    }
                }
                ptr_block->buffer_locks[index_in_64]--;
                return false;
            }
            else {
                ptr_block->buffer_locks[index_in_64]--;
                return false;
            }
        }
    }
    else {
        //8B case
        auto local_head_block_ptr_8b = (SpruceTransVer::AdjSubsequentBlockFive*)local_head_block_ptr;
        if (local_head_block_ptr_8b->bitvector_64 != UINT64_MAX) {
            //value existed
            uint64_t temp_bv_rev = __builtin_bswap64(local_head_block_ptr_8b->bitvector_64.load());
            auto type = local_head_block_ptr_8b->type.load();
            if (type == 5) {
                //type = 5;
                //firstly check the subsequent block
                for (int i = 0; i < 64; i++) {
                    if (!get_bit(&pre_bv, i)) {
                        if (local_head_block_ptr_8b->adj_vertex[i].des == to_node_id) {
//                        local_head_block_ptr_8b->adj_vertex[i].des = UINT32_MAX;    // invalid
                            set_bit(&local_head_block_ptr_8b->bitvector_64, i);
                            if (local_head_block_ptr_8b->bitvector_64.load() == UINT64_MAX) {
                                temp_ptr = local_head_block_ptr_8b->next_block.load();
                                if (!temp_ptr) {
                                    //bottom_head_block empty
                                    bottom_empty_flag = 1;
                                    local_head_block_ptr_8b->obsolete_flag.store(1);
                                }
                                else {
                                    auto next_block_ptr = (uint64_t*)temp_ptr;
                                    if (next_block_ptr[0] == next_block_ptr[1]) {
                                        //delete_num == edge array size
                                        local_head_block_ptr_8b->next_block.store(0);
                                        bottom_empty_flag = 1;
                                        local_head_block_ptr_8b->obsolete_flag.store(1);
                                        // add block_traverse_ptr and edge array ptr to delete queue

                                    }
                                }
                            }
                            temp_weight = local_head_block_ptr_8b -> adj_vertex[i].weight;
                            goto bottom_check;
                        }


                    }
                }

                //check edge array
                temp_ptr = local_head_block_ptr_8b->next_block.load();
                uint64_t* edge_array;
                if (!temp_ptr) {
                    ptr_block->buffer_locks[index_in_64]--;
                    return false;
                }
                else {
                    edge_array = (uint64_t*)temp_ptr;
                }
                auto edge_array_size = edge_array[0];
                auto edges_ptr = (SpruceTransVer::WeightedOutEdge*)(edge_array + 2);

                int64_t low = 0;
                int64_t high = edge_array_size - 1;
                int64_t mid;
                while (high >= low) {
                    mid = (low + high) / 2;
                    if ( edges_ptr[mid].des == UINT64_MAX) {
                        // find actual value
                        uint64_t  new_mid = mid;
                        while (new_mid < high) {
                            new_mid++;
                            if (edges_ptr[new_mid].des != UINT64_MAX) {
                                mid = new_mid;
                                break;
                            }
                        }
                    }
                    if(to_node_id < edges_ptr[mid].des){
                        high = mid - 1;
                    }
                    else if (to_node_id == edges_ptr[mid].des) {
                        edges_ptr[mid].des = UINT64_MAX;    // invalid
                        edge_array[1] += 1; //delete_num += 1
//                        block_traverse_ptr -> next_block_size -= 0;
                        if (edge_array[0] == edge_array[1]) {
                            local_head_block_ptr_8b->next_block.store(0);
                            //add edge array pointer to the delete queue
//                        free(edge_array);
                        }
                        temp_weight = edges_ptr[mid].weight;
                        goto bottom_check;
                    }
                    else {
                        low = mid + 1;
                    }
                }
                ptr_block->buffer_locks[index_in_64]--;
                return false;
            } else  {
                //type = 4;
                for (int i = 0; i < (2 << type); i++) {
                    if (!get_bit(&pre_bv, i)) {
                        if (local_head_block_ptr_8b->adj_vertex[i].des == to_node_id) {
                            set_bit(&local_head_block_ptr_8b->bitvector_64, i);
                            if (local_head_block_ptr_8b->bitvector_64.load() == UINT64_MAX) {
                                //bottom_head_block empty
                                bottom_empty_flag = 1;
                                local_head_block_ptr_8b->obsolete_flag.store(1);
                                //add to delete queue
//                            free(block_traverse_ptr);
                            }
                            temp_weight = local_head_block_ptr_8b -> adj_vertex[i].weight;
                            goto bottom_check;
                        }
                    }
                }
                ptr_block->buffer_locks[index_in_64]--;
                return false;
            }
        } else {
            if (local_head_block_ptr_8b->next_block.load()) {
                //only check edge array (subsequent block is empty)
                temp_ptr = local_head_block_ptr_8b->next_block.load();
                uint64_t* edge_array;
                edge_array = (uint64_t*)temp_ptr;

                auto edge_array_size = edge_array[0];
                auto edges_ptr = (SpruceTransVer::WeightedOutEdge*)(edge_array + 2);

                int64_t low = 0;
                int64_t high = edge_array_size - 1;
                int64_t mid;
                while (high >= low) {
                    mid = (low + high) / 2;
                    if ( edges_ptr[mid].des == UINT64_MAX) {
                        // find actual value
                        uint64_t  new_mid = mid;
                        while (new_mid < high) {
                            new_mid++;
                            if (edges_ptr[new_mid].des != UINT64_MAX) {
                                mid = new_mid;
                                break;
                            }
                        }
                    }
                    if(to_node_id < edges_ptr[mid].des){
                        high = mid - 1;
                    }
                    else if (to_node_id == edges_ptr[mid].des) {
                        edges_ptr[mid].des = UINT64_MAX;    // invalid
                        edge_array[1] += 1; //delete_num += 1
//                        block_traverse_ptr -> next_block_size -= 0;
                        if (edge_array[0] == edge_array[1]) {
                            local_head_block_ptr_8b->next_block.store(0);
                            //add edge array pointer to the delete queue
//                        free(edge_array);
                        }
                        temp_weight = edges_ptr[mid].weight;
                        goto bottom_check;
                    }
                    else {
                        low = mid + 1;
                    }
                }
                ptr_block->buffer_locks[index_in_64]--;
                return false;
            }
            else {
                ptr_block->buffer_locks[index_in_64]--;
                return false;
            }
        }
    }


    bottom_check:
    SpruceTransVer::AddVersion(DELETE_TYPE,
                               reinterpret_cast<SpruceTransVer::AdjSubsequentBlockFive*>(local_head_block_ptr), {to_node_id, temp_weight, local_head_block_ptr->timestamp});
    //edit blocks and bits from bottom to top
    ptr_block->buffer_locks[index_in_64]--;
//    bottom_empty_flag = 0;
    if (bottom_empty_flag == 1 && ((local_head_block_ptr->fb_flag_log_size & 0x7FFF) == 0 )) {
        //get lock and edit bitmap
        while(!middle_block_ptr->mtx[ptr_block_index].compare_exchange_strong(unlocked_m, write_locked_m)) {
            unlocked_m = UNLOCKED;
        }
        clear_bit(&middle_block_ptr->bitmap_8kb, (uint32_t) from_node_id_low);
        middle_block_ptr->mtx[ptr_block_index].store(UNLOCKED);
        if (ptr_block) {
            //corresponding ptr block exists
            ptr_block->ptr_to_buffer[index_in_64].store(0);
        }
        //Check if bitmap == 0;
        bool bv_empty_flag = 1;

        for (unsigned char &i: middle_block_ptr->bitmap_8kb) {
            if (i != 0) {
                bv_empty_flag = 0;
                break;
            }
        }
        if (bv_empty_flag) {
            // lock all middle block locks and recheck
            for (int l = 0; l < 1024; l++) {
                while(!middle_block_ptr->mtx[l].compare_exchange_strong(unlocked_m, write_locked_m)) {
                    unlocked_m = UNLOCKED;
                }
            }
            // recheck
            bv_empty_flag = 1;
            for (unsigned char &i: middle_block_ptr->bitmap_8kb) {
                if (i != 0) {
                    bv_empty_flag = 0;
                    break;
                }
            }
            if (bv_empty_flag) {
                // add middle block to delete queue
                middle_block_ptr->obsolete_flag.store(1);
            }
            for (int l = 0; l < 1024; l++) {
                middle_block_ptr->mtx[l].store(UNLOCKED);
            }
        }
        else {
            return true;
        }

        if (bv_empty_flag) {
            // lock top block and edit bitmap
            while(!root->mtx[from_node_id_high/64].compare_exchange_strong(unlocked_m, write_locked_m)) {
                unlocked_m = UNLOCKED;
            }
            //add m-b to delete queue
//            free(middle_block_ptr);
            clear_bit(&root->bitmap_8kb, from_node_id_high);
            root->ptr_to_children[from_node_id_high].store(0);
            root->mtx[from_node_id_high/64].store(UNLOCKED);
        }
        else {
            return true;
        }

//        free(bottom_head_block);

        // then check topblock similarly
        bv_empty_flag = 1;
        for (unsigned char &i: root->bitmap_8kb) {
            if (i != 0) {
                bv_empty_flag = 0;
                break;
            }
        }
        if (bv_empty_flag) {
            // lock all middle block locks and recheck
            for (int l = 0; l < 1024; l++) {
                while(!root->mtx[l].compare_exchange_strong(unlocked_m, write_locked_m)) {
                    unlocked_m = UNLOCKED;
                }
            }
            // recheck
            bv_empty_flag = 1;
            for (unsigned char &i: root->bitmap_8kb) {
                if (i != 0) {
                    bv_empty_flag = 0;
                    break;
                }
            }
            if (bv_empty_flag) {
                // add middle block to delete queue
                root->obsolete_flag.store(1);
            }
            for (int l = 0; l < 1024; l++) {
                root->mtx[l].store(UNLOCKED);
            }
        }
        else {
            return true;
        }

        if (bv_empty_flag) {
            if (hash_index!=0) {
                spruce.spruce_hash.erase(hash_index);
                free(root);
            }
        }
    } else {
        return true;
    }
    return true;
}

bool SpruceTransVer::AddVersion(uint8_t type, SpruceTransVer::AdjSubsequentBlockFive* local_head_block_ptr,
                                SpruceTransVer::WeightedOutEdge edge) {
//        return true;
    auto log_size = local_head_block_ptr->fb_flag_log_size.load() & 0x7FFF;
    auto inner_offset = log_size % BLOCK_SIZE;
    auto outer_offset = log_size / BLOCK_SIZE;
    auto log_block_ptr = (SpruceTransVer::LogBlock*)local_head_block_ptr->ptr_to_log_block.load();
    SpruceTransVer::LogBlock* temp_ptr;
    if (inner_offset == 0) {
        // need to malloc new blocks
        if (outer_offset == BLOCK_THRESHOLD) {
            // block full, just reuse the last block
            for (int i = 0; i < BLOCK_THRESHOLD - 1; i++) {
                log_block_ptr = log_block_ptr->next_block;
                if (i == BLOCK_THRESHOLD - 3) {
//                        //debug
//                        for (int j = 0; j < BLOCK_SIZE; j++) {
//                            std::cout << log_block_ptr->versioned_edges[j].timestamp << std::endl;
//                        }
                    temp_ptr = log_block_ptr;
                }
            }
            temp_ptr->next_block = NULL;
            memset(log_block_ptr, 0, sizeof(SpruceTransVer::LogBlock));
            log_block_ptr->next_block = (SpruceTransVer::LogBlock*) local_head_block_ptr->ptr_to_log_block.load();
            local_head_block_ptr->fb_flag_log_size -= (BLOCK_SIZE-1); // -64 + 1
            log_block_ptr->versioned_edges[inner_offset] = {(uint64_t)type << 56 | local_head_block_ptr->timestamp,
                                                            edge}; //delete, type 0
            local_head_block_ptr->ptr_to_log_block = (uint64_t)log_block_ptr;
        }
        else {
            // malloc new blocks
            auto log_block_ptr_new = (SpruceTransVer::LogBlock*)malloc(sizeof(SpruceTransVer::LogBlock));
            if (log_block_ptr_new == NULL) return false;
            memset(log_block_ptr_new, 0, sizeof(SpruceTransVer::LogBlock));
            log_block_ptr_new->next_block = log_block_ptr;
            log_block_ptr_new->versioned_edges[inner_offset] = {((uint64_t)type) << 56 | local_head_block_ptr->timestamp, edge};
            local_head_block_ptr->fb_flag_log_size++;
            //debug
//                std::cout << local_head_block_ptr->log_size << std::endl;
            local_head_block_ptr->ptr_to_log_block = (uint64_t)log_block_ptr_new;
        }
    }
    else {
        log_block_ptr->versioned_edges[inner_offset] = {((uint64_t)type) << 56 | local_head_block_ptr->timestamp, edge};
        local_head_block_ptr->fb_flag_log_size++;
    }
    return true;
}

bool SpruceTransVer::UpdateEdge(SpruceTransVer &spruce, SpruceTransVer::WeightedEdge edge) {



    if (edge.weight < 0) {
        SpruceTransVer::DeleteEdge(spruce, edge.src, edge.des);
    }
    SpruceTransVer::TopBlock * root;
    auto from_node_id = edge.src;
    auto to_node_id = edge.des;
    auto weight = edge.weight;
    auto hash_index = (uint32_t)(from_node_id >> 32);
    auto from_node_id_low = (uint16_t)from_node_id;
    auto from_node_id_high = (uint16_t)(from_node_id >> 16);


    //set up lock status
    int unlocked = UNLOCKED, read_locked = READ_LOCKED, write_locked = WRITE_LOCKED;
    uint8_t unlocked_m = UNLOCKED;
    uint8_t write_locked_m = WRITE_LOCKED;

    //need lock: in case for deletion
    SpruceTransVer::MiddleBlock* middle_block_ptr;
    restart:
    if (hash_index || spruce.fb_flag) {
        // if existed in hash table
        root = spruce.spruce_hash.get(hash_index);
        // if not existed
        if (!root) {
            root = SpruceTransVer::CreateTopBlock();
            spruce.spruce_hash.assign(hash_index,root);
            spruce.fb_flag.store(0x01);
        }
    }
    else {
        root = spruce.top_block;
    }

    restart_middle:
    if (root->obsolete_flag) {
        goto restart;
    }
    if (!get_bit(&(root->bitmap_8kb), (uint32_t) from_node_id_high)) {
        //Noted that when thread get lock, corresponding middle block may be created, thereby recheck is needed;
        //try to get write lock (spin lock)
        while(!root->mtx[from_node_id_high/64].compare_exchange_strong(unlocked_m, write_locked_m)) {
            unlocked_m = UNLOCKED;
        }
        if (!get_bit(&(root->bitmap_8kb), (uint32_t) from_node_id_high)) {
            // Need to malloc new middle block
            middle_block_ptr = (SpruceTransVer::MiddleBlock*) malloc(sizeof(SpruceTransVer::MiddleBlock));
            memset(middle_block_ptr, 0, sizeof(SpruceTransVer::MiddleBlock));
            root->ptr_to_children[from_node_id_high].store((uint64_t)middle_block_ptr);
            set_bit(&(root->bitmap_8kb), (uint32_t) from_node_id_high);
            //Analysis
            middle_block_num_++;
        }
        else {
            //atomically load
            uint64_t temp = root->ptr_to_children[from_node_id_high].load();
            middle_block_ptr = (SpruceTransVer::MiddleBlock*)temp;
        }
        root->mtx[from_node_id_high/64]--/*.store(UNLOCKED)*/;
    }
    else {
        //read atomically
        uint64_t temp = root->ptr_to_children[from_node_id_high].load();
        middle_block_ptr = (SpruceTransVer::MiddleBlock*)temp;
    }

    // Get to middle block
    if (!middle_block_ptr) {
        goto restart_middle;
    }
    restart_bottom:
    if (middle_block_ptr->obsolete_flag) {
        goto restart_middle;
    }
    //now we need to get corresponding uint64_t and check the number of 1
    auto ptr_block_index = from_node_id_low / 64;
    auto auxiliary_ptr = reinterpret_cast<uint64_t*>(&middle_block_ptr->bitmap_8kb);
    uint64_t auxiliary_64 = auxiliary_ptr[ptr_block_index];
    uint64_t auxiliary_64_rev = __builtin_bswap64(auxiliary_64);
//    auto ptr_num = __builtin_popcountl(auxiliary_64);
    auto index_in_64 = (uint32_t)(from_node_id_low % 64);

    //Decide the ptr type due to the number of 1

    SpruceTransVer::AdjSubsequentBlockOne* bottom_head_block;

    int lock_flag = 0;
    unlocked = UNLOCKED;

    if (!get_bit(&auxiliary_64, index_in_64)) {
        // does not exist, lock;
        while(!middle_block_ptr->mtx[ptr_block_index].compare_exchange_strong(unlocked_m, write_locked_m)) {
            unlocked_m = UNLOCKED;
        }
        lock_flag = 1;
    }




    //Get bottom head block
    // change set bit sequence for parallel: firstly edit pointer, then set bit
    // reget values!!!!
    auxiliary_64 = auxiliary_ptr[ptr_block_index];
    auxiliary_64_rev = __builtin_bswap64(auxiliary_64);
    SpruceTransVer::PtrBlock* ptr_block;
    // recheck
    if (!get_bit(&auxiliary_64, index_in_64)) {
        //bottom block does not exist, malloc a new block
        bottom_head_block = (SpruceTransVer::AdjSubsequentBlockOne*) malloc(sizeof(SpruceTransVer::AdjSubsequentBlockOne));
        memset(bottom_head_block, 0, sizeof(SpruceTransVer::AdjSubsequentBlockOne));
        bottom_head_block->bitvector_64 = UINT64_MAX;
        bottom_head_block->type.store(1);

        //Edit middle block bitmap and ptr block
        uint64_t temp = middle_block_ptr->ptr_to_children[ptr_block_index].load();
        ptr_block = (SpruceTransVer::PtrBlock*)temp;

        if(!ptr_block) {
            // + 1 for obsolete flag
            auto new_ptr_block = (SpruceTransVer::PtrBlock*) malloc(sizeof(SpruceTransVer::PtrBlock)) ;
            memset(new_ptr_block, 0, sizeof(SpruceTransVer::PtrBlock));
            new_ptr_block->ptr_to_buffer[index_in_64].store(reinterpret_cast<unsigned long>(bottom_head_block));
            middle_block_ptr->ptr_to_children[ptr_block_index].store((uint64_t)new_ptr_block);
            ptr_block = new_ptr_block;
            sparse_ptr_block_num++;
        }
        else {
            if (ptr_block->obsolete_flag) {
                //obsoleted
                if (lock_flag) {
                    middle_block_ptr->mtx[ptr_block_index]--;
                }
                goto restart_bottom;
            }
            ptr_block->ptr_to_buffer[index_in_64].store(reinterpret_cast<unsigned long>(bottom_head_block));
        }


        set_bit(&(middle_block_ptr->bitmap_8kb), (uint32_t) from_node_id_low);
        //Analysis
        type1++;
    } else {
        //corresponding block exists
        uint64_t temp = middle_block_ptr->ptr_to_children[ptr_block_index].load();
        ptr_block = (SpruceTransVer::PtrBlock*)temp;
    }



    unlocked = UNLOCKED;
    if (lock_flag) {
        middle_block_ptr->mtx[ptr_block_index]--;
    }


    uint32_t insert_index;
    uint8_t type;
    uint64_t temp_ptr;

    while(!ptr_block->buffer_locks[index_in_64].compare_exchange_strong(unlocked_m, write_locked_m)) {
        unlocked_m = UNLOCKED;
    }
    auto local_head_block_ptr = (SpruceTransVer::AdjSubsequentBlockFive4B*)ptr_block->ptr_to_buffer[index_in_64].load();
    if (!local_head_block_ptr) {
        ptr_block->buffer_locks[index_in_64]--;
        goto restart_bottom;
    }
    if (local_head_block_ptr->obsolete_flag) {
        ptr_block->buffer_locks[index_in_64]--;
        goto restart_bottom;
    }

    auto pre_bv = local_head_block_ptr->bitvector_64.load();

    //Edit timestamp
    local_head_block_ptr->timestamp += 1;

    //firstly check whether edge exists
    //value existed
    uint64_t temp_bv_rev = __builtin_bswap64(local_head_block_ptr->bitvector_64.load());
    type = local_head_block_ptr->type.load();

    if (!(local_head_block_ptr->fb_flag_log_size >> 15)) {
        //4B case
        if (type == 5) {
            //type = 5;
            //firstly check the subsequent block
            for (int i = 0; i < 64; i++) {
                if (!get_bit(&pre_bv, i)) {
                    if (local_head_block_ptr->adj_vertex[i].des == (uint32_t)to_node_id) {
                        AddVersion(UPDATE_TYPE,
                                   reinterpret_cast<SpruceTransVer::AdjSubsequentBlockFive*>(local_head_block_ptr), {to_node_id, local_head_block_ptr->adj_vertex[i].weight, local_head_block_ptr->timestamp});
                        local_head_block_ptr->adj_vertex[i].weight = weight;
                        ptr_block->buffer_locks[index_in_64]--;
                        return true;
                    }
                }
            }

            //check edge array
            temp_ptr = local_head_block_ptr->next_block.load();
            uint32_t* edge_array;
            if (!temp_ptr) {
                goto insert_update;
            }
            else {
                edge_array = (uint32_t*)temp_ptr;
            }
            auto edge_array_size = edge_array[0];
            auto edges_ptr = (SpruceTransVer::WeightedOutEdge4B*)(edge_array + 2);

            int64_t low = 0;
            int64_t high = edge_array_size - 1;
            int64_t mid;
            while (high >= low) {
                mid = (low + high) / 2;
                if ( edges_ptr[mid].des == UINT32_MAX) {
                    // find actual value
                    uint32_t  new_mid = mid;
                    while (new_mid < high) {
                        new_mid++;
                        if (edges_ptr[new_mid].des != UINT32_MAX) {
                            mid = new_mid;
                            break;
                        }
                    }
                }
                if(to_node_id < edges_ptr[mid].des){
                    high = mid - 1;
                }
                else if (to_node_id == edges_ptr[mid].des) {
                    edges_ptr[mid].weight = weight;
                    AddVersion(UPDATE_TYPE,
                               reinterpret_cast<SpruceTransVer::AdjSubsequentBlockFive*>(local_head_block_ptr), {to_node_id, edges_ptr[mid].weight, local_head_block_ptr->timestamp});
                    ptr_block->buffer_locks[index_in_64]--;
                    return true;
                }
                else {
                    low = mid + 1;
                }
            }
            goto insert_update;
        } else {
            for (int i = 0; i < (2 << type); i++) {
                if (!get_bit(&pre_bv, i)) {
                    if (local_head_block_ptr->adj_vertex[i].des == (uint32_t)to_node_id) {
                        local_head_block_ptr->adj_vertex[i].weight = weight;
                        ptr_block->buffer_locks[index_in_64]--;
                        return true;                    }
                }
            };
            goto insert_update;
        }
    }
    else {
        //8B case
        auto local_head_block_ptr_8b = (AdjSubsequentBlockFive*)local_head_block_ptr;
        if (type == 5) {
            //type = 5;
            //firstly check the subsequent block
            for (int i = 0; i < 64; i++) {
                if (!get_bit(&pre_bv, i)) {
                    if (local_head_block_ptr_8b->adj_vertex[i].des == to_node_id) {
                        AddVersion(UPDATE_TYPE,
                                   reinterpret_cast<SpruceTransVer::AdjSubsequentBlockFive*>(local_head_block_ptr), {to_node_id,  local_head_block_ptr->adj_vertex[i].weight, local_head_block_ptr->timestamp});
                        local_head_block_ptr_8b->adj_vertex[i].weight = weight;
                        ptr_block->buffer_locks[index_in_64]--;
                        return true;
                    }
                }
            }

            //check edge array
            temp_ptr = local_head_block_ptr_8b->next_block.load();
            uint64_t* edge_array;
            if (!temp_ptr) {
                goto insert_update;
            }
            else {
                edge_array = (uint64_t*)temp_ptr;
            }
            auto edge_array_size = edge_array[0];
            auto edges_ptr = (SpruceTransVer::WeightedOutEdge*)(edge_array + 2);

            int64_t low = 0;
            int64_t high = edge_array_size - 1;
            int64_t mid;
            while (high >= low) {
                mid = (low + high) / 2;
                if ( edges_ptr[mid].des == UINT32_MAX) {
                    // find actual value
                    uint32_t  new_mid = mid;
                    while (new_mid < high) {
                        new_mid++;
                        if (edges_ptr[new_mid].des != UINT32_MAX) {
                            mid = new_mid;
                            break;
                        }
                    }
                }
                if(to_node_id < edges_ptr[mid].des){
                    high = mid - 1;
                }
                else if (to_node_id == edges_ptr[mid].des) {
                    edges_ptr[mid].weight = weight;
                    AddVersion(UPDATE_TYPE,
                               reinterpret_cast<SpruceTransVer::AdjSubsequentBlockFive*>(local_head_block_ptr), {to_node_id, edges_ptr[mid].weight, local_head_block_ptr->timestamp});
                    ptr_block->buffer_locks[index_in_64]--;
                    return true;
                }
                else {
                    low = mid + 1;
                }
            }
            goto insert_update;
        } else {
            for (int i = 0; i < (2 << type); i++) {
                if (!get_bit(&pre_bv, i)) {
                    if (local_head_block_ptr_8b->adj_vertex[i].des == to_node_id) {
                        local_head_block_ptr_8b->adj_vertex[i].weight = weight;
                        ptr_block->buffer_locks[index_in_64]--;
                        return true;                    }
                }
            };
            goto insert_update;
        }
    }




    //head insertion
    insert_update:
    if ((local_head_block_ptr->fb_flag_log_size>>15)^spruce.fb_flag) {
        // 4B->8B case
        WeightedOutEdge out_edge = {edge.des, edge.weight, local_head_block_ptr->timestamp};
        SpruceTransVer::AdjSubsequentBlockFive* local_head_block_ptr_8b;
        if (local_head_block_ptr->bitvector_64 != 0) {
            //not full yet
            uint64_t temp_bv_rev = __builtin_bswap64(local_head_block_ptr->bitvector_64.load());
            insert_index = __builtin_clzl(temp_bv_rev);
            type = local_head_block_ptr->type;
            clear_bit(&local_head_block_ptr->bitvector_64, insert_index);
            uint64_t temp;
            if (insert_index > ((1 << (type - 1)) * 4 - 1)) {
                //need to malloc new space;
                switch (type) {
                    case 1: {
                        auto new_block = (SpruceTransVer::AdjSubsequentBlockTwo*) malloc(
                                sizeof(AdjSubsequentBlockTwo));
                        memset(new_block, 0, sizeof(AdjSubsequentBlockTwo));
                        temp = (uint64_t)local_head_block_ptr;

                        //copy values
//                        memcpy(new_block,
//                               (SpruceTransVer::AdjSubsequentBlockOne*)(temp),
//                               sizeof(SpruceTransVer::AdjSubsequentBlockOne));
                        new_block->type = local_head_block_ptr->type + 1;
                        new_block->fb_flag_log_size = local_head_block_ptr->fb_flag_log_size.load() & 0x7FFF;
                        new_block->bitvector_64 = local_head_block_ptr->bitvector_64.load();
                        new_block->timestamp = local_head_block_ptr->timestamp.load();
                        for (int c = 0; c < ((1 << (type - 1)) * 4); c++) {
                            new_block->adj_vertex[c] = {local_head_block_ptr->adj_vertex->des, local_head_block_ptr->adj_vertex->weight};
                        }

                        new_block->adj_vertex[insert_index] = out_edge;
                        local_head_block_ptr->obsolete_flag = 1;
                        // add previous block to delete queue
                        ptr_block->ptr_to_buffer[index_in_64].store(reinterpret_cast<uint64_t>(new_block));
                        free((SpruceTransVer::AdjSubsequentBlockOne4B*)(temp));
                        local_head_block_ptr_8b = (SpruceTransVer::AdjSubsequentBlockFive*)new_block;
                        //analysis
                        type2++;
                        type1--;

                        break;
                    }
                    case 2: {
                        auto new_block = (SpruceTransVer::AdjSubsequentBlockThree*) malloc(
                                sizeof(AdjSubsequentBlockThree));
                        memset(new_block, 0, sizeof(AdjSubsequentBlockThree));
                        temp = (uint64_t)local_head_block_ptr;
                        //copy values
                        new_block->type = local_head_block_ptr->type + 1;
                        new_block->fb_flag_log_size = local_head_block_ptr->fb_flag_log_size.load() & 0x7FFF;
                        new_block->bitvector_64 = local_head_block_ptr->bitvector_64.load();
                        new_block->timestamp = local_head_block_ptr->timestamp.load();
                        for (int c = 0; c < ((1 << (type - 1)) * 4); c++) {
                            new_block->adj_vertex[c] = {local_head_block_ptr->adj_vertex->des, local_head_block_ptr->adj_vertex->weight};
                        }
                        new_block->adj_vertex[insert_index] = out_edge;
                        new_block->type.store(3);
                        local_head_block_ptr->obsolete_flag = 1;
                        // add previous block to delete queue
                        ptr_block->ptr_to_buffer[index_in_64].store(reinterpret_cast<uint64_t>(new_block));
                        free((SpruceTransVer::AdjSubsequentBlockTwo*)temp);
                        //analysis
                        type3++;
                        type2--;

                        break;
                    }
                    case 3: {
                        auto new_block = (SpruceTransVer::AdjSubsequentBlockFour*) malloc(
                                sizeof(AdjSubsequentBlockFour));
                        memset(new_block, 0, sizeof(AdjSubsequentBlockFour));
                        temp = (uint64_t)local_head_block_ptr;
                        //copy values
                        new_block->type = local_head_block_ptr->type + 1;
                        new_block->fb_flag_log_size = local_head_block_ptr->fb_flag_log_size.load() & 0x7FFF;
                        new_block->bitvector_64 = local_head_block_ptr->bitvector_64.load();
                        new_block->timestamp = local_head_block_ptr->timestamp.load();
                        for (int c = 0; c < ((1 << (type - 1)) * 4); c++) {
                            new_block->adj_vertex[c] = {local_head_block_ptr->adj_vertex->des, local_head_block_ptr->adj_vertex->weight};
                        }
                        new_block->adj_vertex[insert_index] = out_edge;
                        local_head_block_ptr->obsolete_flag = 1;
                        // add previous block to delete queue
                        ptr_block->ptr_to_buffer[index_in_64].store(reinterpret_cast<uint64_t>(new_block));
                        free((SpruceTransVer::AdjSubsequentBlockThree*)temp);
                        //analysis
                        type4++;
                        type3--;

                        break;

                    }
                    case 4: {
                        auto new_block = (SpruceTransVer::AdjSubsequentBlockFive*) malloc(
                                sizeof(AdjSubsequentBlockFive));
                        temp = (uint64_t)local_head_block_ptr;
                        memset(new_block, 0, sizeof(AdjSubsequentBlockFive));
                        // new_block->bitvector_64 = UINT64_MAX; //no subsequent block
                        //copy values
                        new_block->type = local_head_block_ptr->type + 1;
                        new_block->fb_flag_log_size = local_head_block_ptr->fb_flag_log_size.load() & 0x7FFF;
                        new_block->bitvector_64 = local_head_block_ptr->bitvector_64.load();
                        new_block->timestamp = local_head_block_ptr->timestamp.load();
                        for (int c = 0; c < ((1 << (type - 1)) * 4); c++) {
                            new_block->adj_vertex[c] = {local_head_block_ptr->adj_vertex->des, local_head_block_ptr->adj_vertex->weight};
                        }
                        new_block->adj_vertex[insert_index] = out_edge;
                        local_head_block_ptr->obsolete_flag = 1;
                        // add previous block to delete queue
                        ptr_block->ptr_to_buffer[index_in_64].store(reinterpret_cast<uint64_t>(new_block));
                        free((SpruceTransVer::AdjSubsequentBlockFour*)temp);
                        //analysis
                        type5++;
                        type4--;

                        break;
                    }
                    default: {
                        std::cout << "Impossible type!" << std::endl;
                        break;
                    }

                }
            } else {
                // just change type
                switch (type) {
                    case 1: {
                        auto new_block = (SpruceTransVer::AdjSubsequentBlockOne*) malloc(
                                sizeof(AdjSubsequentBlockOne));
                        memset(new_block, 0, sizeof(AdjSubsequentBlockOne));
                        temp = (uint64_t)local_head_block_ptr;

                        //copy values
//                        memcpy(new_block,
//                               (SpruceTransVer::AdjSubsequentBlockOne*)(temp),
//                               sizeof(SpruceTransVer::AdjSubsequentBlockOne));
                        new_block->type = local_head_block_ptr->type.load();
                        new_block->fb_flag_log_size = local_head_block_ptr->fb_flag_log_size.load() & 0x7FFF;
                        new_block->bitvector_64 = local_head_block_ptr->bitvector_64.load();
                        new_block->timestamp = local_head_block_ptr->timestamp.load();
                        for (int c = 0; c < ((1 << (type - 1)) * 4); c++) {
                            new_block->adj_vertex[c] = {local_head_block_ptr->adj_vertex->des, local_head_block_ptr->adj_vertex->weight};
                        }

                        new_block->adj_vertex[insert_index] = out_edge;
                        local_head_block_ptr->obsolete_flag = 1;
                        // add previous block to delete queue
                        ptr_block->ptr_to_buffer[index_in_64].store(reinterpret_cast<uint64_t>(new_block));
                        free((SpruceTransVer::AdjSubsequentBlockOne4B*)(temp));
                        local_head_block_ptr_8b = (SpruceTransVer::AdjSubsequentBlockFive*)new_block;
                        break;
                    }
                    case 2: {
                        auto new_block = (SpruceTransVer::AdjSubsequentBlockTwo*) malloc(
                                sizeof(AdjSubsequentBlockTwo));
                        memset(new_block, 0, sizeof(AdjSubsequentBlockTwo));
                        temp = (uint64_t)local_head_block_ptr;
                        //copy values
                        new_block->type = local_head_block_ptr->type.load();
                        new_block->fb_flag_log_size = local_head_block_ptr->fb_flag_log_size.load() & 0x7FFF;
                        new_block->bitvector_64 = local_head_block_ptr->bitvector_64.load();
                        new_block->timestamp = local_head_block_ptr->timestamp.load();
                        for (int c = 0; c < ((1 << (type - 1)) * 4); c++) {
                            new_block->adj_vertex[c] = {local_head_block_ptr->adj_vertex->des, local_head_block_ptr->adj_vertex->weight};
                        }
                        new_block->adj_vertex[insert_index] = out_edge;
                        new_block->type.store(3);
                        local_head_block_ptr->obsolete_flag = 1;
                        // add previous block to delete queue
                        ptr_block->ptr_to_buffer[index_in_64].store(reinterpret_cast<uint64_t>(new_block));
                        free((SpruceTransVer::AdjSubsequentBlockTwo*)temp);
                        break;
                    }
                    case 3: {
                        auto new_block = (SpruceTransVer::AdjSubsequentBlockThree*) malloc(
                                sizeof(AdjSubsequentBlockThree));
                        memset(new_block, 0, sizeof(AdjSubsequentBlockThree));
                        temp = (uint64_t)local_head_block_ptr;
                        //copy values
                        new_block->type = local_head_block_ptr->type.load();
                        new_block->fb_flag_log_size = local_head_block_ptr->fb_flag_log_size.load() & 0x7FFF;
                        new_block->bitvector_64 = local_head_block_ptr->bitvector_64.load();
                        new_block->timestamp = local_head_block_ptr->timestamp.load();
                        for (int c = 0; c < ((1 << (type - 1)) * 4); c++) {
                            new_block->adj_vertex[c] = {local_head_block_ptr->adj_vertex->des, local_head_block_ptr->adj_vertex->weight};
                        }
                        new_block->adj_vertex[insert_index] = out_edge;
                        local_head_block_ptr->obsolete_flag = 1;
                        // add previous block to delete queue
                        ptr_block->ptr_to_buffer[index_in_64].store(reinterpret_cast<uint64_t>(new_block));
                        free((SpruceTransVer::AdjSubsequentBlockThree*)temp);
                        //analysis
                        break;

                    }
                    case 4: {
                        auto new_block = (SpruceTransVer::AdjSubsequentBlockFour*) malloc(
                                sizeof(AdjSubsequentBlockFour));
                        temp = (uint64_t)local_head_block_ptr;
                        memset(new_block, 0, sizeof(AdjSubsequentBlockFour));
                        //copy values
                        new_block->type = local_head_block_ptr->type.load();
                        new_block->fb_flag_log_size = local_head_block_ptr->fb_flag_log_size.load() & 0x7FFF;
                        new_block->bitvector_64 = local_head_block_ptr->bitvector_64.load();
                        new_block->timestamp = local_head_block_ptr->timestamp.load();
                        for (int c = 0; c < ((1 << (type - 1)) * 4); c++) {
                            new_block->adj_vertex[c] = {local_head_block_ptr->adj_vertex->des, local_head_block_ptr->adj_vertex->weight};
                        }
                        new_block->adj_vertex[insert_index] = out_edge;
                        local_head_block_ptr->obsolete_flag = 1;
                        // add previous block to delete queue
                        ptr_block->ptr_to_buffer[index_in_64].store(reinterpret_cast<uint64_t>(new_block));
                        free((SpruceTransVer::AdjSubsequentBlockFour*)temp);
                        //analysis
                        break;
                    }
                    case 5: {
                        auto new_block = (SpruceTransVer::AdjSubsequentBlockFive*) malloc(
                                sizeof(AdjSubsequentBlockFive));
                        temp = (uint64_t)local_head_block_ptr;
                        memset(new_block, 0, sizeof(AdjSubsequentBlockFive));
                        //copy values
                        new_block->type = local_head_block_ptr->type.load();
                        new_block->fb_flag_log_size = local_head_block_ptr->fb_flag_log_size.load() & 0x7FFF;
                        new_block->bitvector_64 = local_head_block_ptr->bitvector_64.load();
                        new_block->timestamp = local_head_block_ptr->timestamp.load();
                        for (int c = 0; c < ((1 << (type - 1)) * 4); c++) {
                            new_block->adj_vertex[c] = {local_head_block_ptr->adj_vertex->des, local_head_block_ptr->adj_vertex->weight};
                        }
                        new_block->adj_vertex[insert_index] = out_edge;
                        local_head_block_ptr->obsolete_flag = 1;

                        //also change the edge array from 4B->8B
                        uint32_t old_block_size, old_delete_num;
                        uint64_t temp1 = (uint64_t)local_head_block_ptr;
                        uint64_t temp2 = local_head_block_ptr->next_block.load();
                        auto old_edge_array = (uint32_t*)temp2;
                        if (old_edge_array) {
                            // does not exist, set initial size
                            old_block_size = old_edge_array[0];
                            old_delete_num = old_edge_array[1];
                            auto old_edges = (SpruceTransVer::WeightedOutEdge4B*) (old_edge_array + 2);
                            int64_t old_delete_64 = 0;
                            if (old_delete_num > 64) {
                                old_delete_64 = old_delete_num / 64;
                            }

                            // + 1 for block_size and delete_num (equals to weighted out edge)
                            // resize block according to delete_num
                            auto new_block_size = old_block_size + 1;
                            auto new_edge_array = (uint64_t*) malloc(
                                    sizeof(SpruceTransVer::WeightedOutEdge) * (new_block_size));
                            memset(new_edge_array, 0, sizeof(SpruceTransVer::WeightedOutEdge) * (new_block_size));

                            //then use merge sort to place spaces in new block
                            uint64_t k = 0;
                            uint64_t start1 = 0, start2 = 0;
                            uint64_t end1 = 64, end2 = old_block_size;

                            auto new_edges = (SpruceTransVer::WeightedOutEdge*) (new_edge_array + 2);

                            while (start2 < end2) {
                                if (old_edges[start2].des == UINT32_MAX) {
                                    // skip invalid data
                                    new_edges[k++] = convert_to_8b(old_edges[start2++]);
                                    new_edges[k++].des = UINT64_MAX;
                                    continue;
                                }
                                new_edges[k++] = convert_to_8b(old_edges[start2++]);
                            }

                            //reset subsequent block status
                            new_edge_array[0] = new_block_size -
                                                1; // remember: block_size does not include first 2 uint32_t for information
                            new_edge_array[1] = old_delete_num;  // copy delete_num
                            new_block->next_block.store(reinterpret_cast<uint64_t>(new_edge_array));
                            //do not execute free function in new;
                            // Noted that when execute on large datasets, use free function to avoid memory exceeded
                            // when execute parallel test, comment it to avoid errors
                            free(old_edge_array);
                        }
                        // add previous block to delete queue
                        ptr_block->ptr_to_buffer[index_in_64].store(reinterpret_cast<uint64_t>(new_block));
                        free((SpruceTransVer::AdjSubsequentBlockFour*)temp);
                        break;
                    }
                    default: {
                        std::cout << "Impossible type!" << std::endl;
                        break;
                    }
                }
            }
        } else {
            type = 5;
            // while ends, no spaces
            // Firstly check the block size,we assume that we only check invalid label in deletion
            uint32_t old_block_size, old_delete_num;
            uint64_t temp1 = (uint64_t)local_head_block_ptr;
            uint64_t temp2 = local_head_block_ptr->next_block.load();
            auto old_edge_array = (uint32_t*)temp2;
            if (!old_edge_array) {
                // does not exist, set initial size

                old_block_size = 0;
                old_delete_num = 0;
                // analysis
                edge_array_num++;
            }
            else {
                old_block_size = old_edge_array[0];
                old_delete_num = old_edge_array[1];
            }

            auto old_edges = (SpruceTransVer::WeightedOutEdge4B*)(old_edge_array + 2);

            int64_t old_delete_64 = 0;
            if (old_delete_num > 64) {
                old_delete_64 = old_delete_num/64;
            }

            // + 1 for block_size and delete_num (equals to weighted out edge)
            // resize block according to delete_num
            auto new_block_size = (64 + 1 + old_block_size - old_delete_64 * 64);
            auto new_edge_array = (uint64_t*) malloc(
                    sizeof(SpruceTransVer::WeightedOutEdge) * (new_block_size));
            memset(new_edge_array, 0, sizeof(SpruceTransVer::WeightedOutEdge) * (new_block_size));
//

            //sort the vertex using bubble sort
            for (int i = 0; i < 64; i++) {
                bool flag = false;
                for (int j = 0; j < 64 - i - 1; j++) {
                    if (local_head_block_ptr->adj_vertex[j].des > local_head_block_ptr->adj_vertex[j + 1].des) {
                        auto temp = local_head_block_ptr->adj_vertex[j];
                        local_head_block_ptr->adj_vertex[j] = local_head_block_ptr->adj_vertex[j + 1];
                        local_head_block_ptr->adj_vertex[j + 1] = temp;
                        flag = true;
                    }
                }
                if (!flag) {
                    break;
                }
            }

            //then use merge sort to place spaces in new block
            uint64_t k = 0;
            uint64_t start1 = 0, start2 = 0;
            uint64_t end1 = 64, end2 = old_block_size;

            auto new_edges = (SpruceTransVer::WeightedOutEdge*)(new_edge_array + 2);

            while (start1 < end1 && start2 < end2) {
                if (old_edges[start2].des == UINT32_MAX) {
                    // skip invalid data
                    start2++;
                    continue;
                }
                new_edges[k++] = local_head_block_ptr->adj_vertex[start1].des < old_edges[start2].des ?
                                 convert_to_8b(local_head_block_ptr->adj_vertex[start1++]) : convert_to_8b(old_edges[start2++]);
            }
            while (start1 < end1) {
                new_edges[k++] = convert_to_8b(local_head_block_ptr->adj_vertex[start1++]);
            }
            while (start2 < end2) {
                if (old_edges[start2].des == UINT32_MAX) {
                    // skip invalid data
                    start2++;
                    continue;
                }
                new_edges[k++] = convert_to_8b(old_edges[start2++]);
            }
            while (k < new_block_size - 1) {
                // shift invalid values with size < 64
                new_edges[k++].des = UINT64_MAX;
            }

            auto new_block = (SpruceTransVer::AdjSubsequentBlockFive*) malloc(
                    sizeof(AdjSubsequentBlockFive));
            auto temp = (uint64_t)local_head_block_ptr;
            memset(new_block, 0, sizeof(AdjSubsequentBlockFive));
            //copy values
            new_block->type = local_head_block_ptr->type.load();
            new_block->fb_flag_log_size = local_head_block_ptr->fb_flag_log_size.load() & 0x7FFF;
            new_block->bitvector_64.store(UINT64_MAX);
            new_block->timestamp = local_head_block_ptr->timestamp.load();
            local_head_block_ptr->obsolete_flag = 1;
            // add previous block to delete queue
            ptr_block->ptr_to_buffer[index_in_64].store(reinterpret_cast<uint64_t>(new_block));
            free((SpruceTransVer::AdjSubsequentBlockFour*)temp);

            //reset subsequent block status
            new_edge_array[0] = new_block_size - 1; // remember: block_size does not include first 2 uint32_t for information
            new_edge_array[1] = old_delete_num % 64;  // copy delete_num
            local_head_block_ptr->next_block.store(reinterpret_cast<uint64_t>(new_edge_array));

            // Generally you need to use your own garbage collecting function in your system instead of simply free it.
            if(old_edge_array) free(old_edge_array);

            //Then insert new edge to subsequent block
            new_block->adj_vertex[0] = out_edge;
            clear_bit(&new_block->bitvector_64, 0);

            //analysis
            type64++;

        }
    }
    else if (!spruce.fb_flag) {
        //4B case
        WeightedOutEdge4B out_edge = {(uint32_t)edge.des, (uint32_t)local_head_block_ptr->timestamp, weight};
        if (local_head_block_ptr->bitvector_64 != 0) {
            //not full yet
            uint64_t temp_bv_rev = __builtin_bswap64(local_head_block_ptr->bitvector_64.load());
            insert_index = __builtin_clzl(temp_bv_rev);
            type = local_head_block_ptr->type;
            clear_bit(&local_head_block_ptr->bitvector_64, insert_index);
            uint64_t temp;
            if (insert_index > ((1 << (type - 1)) * 4 - 1)) {
                //need to malloc new space;
                switch (type) {
                    case 1: {
                        auto new_block = (SpruceTransVer::AdjSubsequentBlockTwo4B*) malloc(
                                sizeof(AdjSubsequentBlockTwo4B));
                        memset(new_block, 0, sizeof(AdjSubsequentBlockTwo4B));
                        temp = (uint64_t)local_head_block_ptr;
                        memcpy(new_block,
                               (SpruceTransVer::AdjSubsequentBlockOne4B*)(temp),
                               sizeof(SpruceTransVer::AdjSubsequentBlockOne4B));
                        new_block->adj_vertex[insert_index] = out_edge;
                        new_block->type.store(2);
                        local_head_block_ptr->obsolete_flag = 1;
                        // add previous block to delete queue
                        ptr_block->ptr_to_buffer[index_in_64].store(reinterpret_cast<uint64_t>(new_block));
//                        free((SpruceTransVer::AdjSubsequentBlockOne4B*)(temp));
                        local_head_block_ptr = (SpruceTransVer::AdjSubsequentBlockFive4B*)new_block;
                        //analysis
                        type2++;
                        type1--;

                        break;
                    }
                    case 2: {
                        auto new_block = (SpruceTransVer::AdjSubsequentBlockThree4B*) malloc(
                                sizeof(AdjSubsequentBlockThree4B));
                        memset(new_block, 0, sizeof(AdjSubsequentBlockThree4B));
                        temp = (uint64_t)local_head_block_ptr;
                        memcpy(new_block,
                               (SpruceTransVer::AdjSubsequentBlockTwo4B*)temp,
                               sizeof(SpruceTransVer::AdjSubsequentBlockTwo4B));
                        new_block->adj_vertex[insert_index] = out_edge;
                        new_block->type.store(3);
                        local_head_block_ptr->obsolete_flag = 1;
                        // add previous block to delete queue
                        ptr_block->ptr_to_buffer[index_in_64].store(reinterpret_cast<uint64_t>(new_block));
//                        free((SpruceTransVer::AdjSubsequentBlockTwo4B*)temp);
                        //analysis
                        type3++;
                        type2--;

                        break;
                    }
                    case 3: {
                        auto new_block = (SpruceTransVer::AdjSubsequentBlockFour4B*) malloc(
                                sizeof(AdjSubsequentBlockFour4B));
                        memset(new_block, 0, sizeof(AdjSubsequentBlockFour4B));
                        temp = (uint64_t)local_head_block_ptr;
                        memcpy(new_block,
                               (SpruceTransVer::AdjSubsequentBlockThree4B*)temp,
                               sizeof(SpruceTransVer::AdjSubsequentBlockThree4B));
                        new_block->adj_vertex[insert_index] = out_edge;
                        new_block->type.store(4);
                        local_head_block_ptr->obsolete_flag = 1;
                        // add previous block to delete queue
                        ptr_block->ptr_to_buffer[index_in_64].store(reinterpret_cast<uint64_t>(new_block));
//                        free((SpruceTransVer::AdjSubsequentBlockThree4B*)temp);
                        //analysis
                        type4++;
                        type3--;

                        break;

                    }
                    case 4: {
                        auto new_block = (SpruceTransVer::AdjSubsequentBlockFive4B*) malloc(
                                sizeof(AdjSubsequentBlockFive4B));
                        temp = (uint64_t)local_head_block_ptr;
                        memset(new_block, 0, sizeof(AdjSubsequentBlockFive4B));
                        // new_block->bitvector_64 = UINT64_MAX; //no subsequent block
                        memcpy(new_block,
                               (SpruceTransVer::AdjSubsequentBlockFour4B*)(temp),
                               sizeof(SpruceTransVer::AdjSubsequentBlockFour4B));
                        new_block->adj_vertex[insert_index] = out_edge;
                        new_block->type.store(5);
                        local_head_block_ptr->obsolete_flag = 1;
                        // add previous block to delete queue
                        ptr_block->ptr_to_buffer[index_in_64].store(reinterpret_cast<uint64_t>(new_block));
//                        free((SpruceTransVer::AdjSubsequentBlockFour4B*)temp);
                        //analysis
                        type5++;
                        type4--;

                        break;
                    }
                    default: {
                        std::cout << "Impossible type!" << std::endl;
                        break;
                    }

                }
            } else {
                local_head_block_ptr->adj_vertex[insert_index] = out_edge;
            }
        } else {
            type = 5;
            // while ends, no spaces
            // Firstly check the block size,we assume that we only check invalid label in deletion
            uint32_t old_block_size, old_delete_num;
            uint64_t temp1 = (uint64_t)local_head_block_ptr;
            uint64_t temp2 = local_head_block_ptr->next_block.load();
            auto old_edge_array = (uint32_t*)temp2;
            if (!old_edge_array) {
                // does not exist, set initial size

                old_block_size = 0;
                old_delete_num = 0;
                // analysis
                edge_array_num++;
            }
            else {
                old_block_size = old_edge_array[0];
                old_delete_num = old_edge_array[1];
            }

            auto old_edges = (SpruceTransVer::WeightedOutEdge4B*)(old_edge_array + 2);

            uint32_t old_delete_64 = 0;
            if (old_delete_num > 64) {
                old_delete_64 = old_delete_num/64;
            }

            // + 1 for block_size and delete_num (equals to weighted out edge)
            // resize block according to delete_num
            auto new_block_size = (64 + 1 + old_block_size - old_delete_64 * 64);
            auto new_edge_array = (uint32_t*) malloc(
                    sizeof(SpruceTransVer::WeightedOutEdge4B) * (new_block_size));
            memset(new_edge_array, 0, sizeof(SpruceTransVer::WeightedOutEdge4B) * (new_block_size));
//

            //sort the vertex using bubble sort
            for (int i = 0; i < 64; i++) {
                bool flag = false;
                for (int j = 0; j < 64 - i - 1; j++) {
                    if (local_head_block_ptr->adj_vertex[j].des > local_head_block_ptr->adj_vertex[j + 1].des) {
                        auto temp = local_head_block_ptr->adj_vertex[j];
                        local_head_block_ptr->adj_vertex[j] = local_head_block_ptr->adj_vertex[j + 1];
                        local_head_block_ptr->adj_vertex[j + 1] = temp;
                        flag = true;
                    }
                }
                if (!flag) {
                    break;
                }
            }

            //then use merge sort to place spaces in new block
            uint32_t k = 0;
            uint32_t start1 = 0, start2 = 0;
            uint32_t end1 = 64, end2 = old_block_size;

            auto new_edges = (SpruceTransVer::WeightedOutEdge4B*)(new_edge_array + 2);

            while (start1 < end1 && start2 < end2) {
                if (old_edges[start2].des == UINT32_MAX) {
                    // skip invalid data
                    start2++;
                    continue;
                }
                new_edges[k++] = local_head_block_ptr->adj_vertex[start1].des < old_edges[start2].des ?
                                 local_head_block_ptr->adj_vertex[start1++] : old_edges[start2++];
            }
            while (start1 < end1) {
                new_edges[k++] = local_head_block_ptr->adj_vertex[start1++];
            }
            while (start2 < end2) {
                if (old_edges[start2].des == UINT32_MAX) {
                    // skip invalid data
                    start2++;
                    continue;
                }
                new_edges[k++] = old_edges[start2++];
            }
            while (k < new_block_size - 1) {
                // shift invalid values with size < 64
                new_edges[k++].des = UINT32_MAX;
            }


            //reset head block status
            local_head_block_ptr->bitvector_64.store(UINT64_MAX);

            //reset subsequent block status
            memset(local_head_block_ptr->adj_vertex, 0, sizeof(SpruceTransVer::WeightedOutEdge4B)*64);
            new_edge_array[0] = new_block_size - 1; // remember: block_size does not include first 2 uint32_t for information
            new_edge_array[1] = old_delete_num % 64;  // copy delete_num
            local_head_block_ptr->next_block.store(reinterpret_cast<uint64_t>(new_edge_array));
            // Generally you need to use your own garbage collecting function in your system instead of simply free itrrors
//            if(old_edge_array) free(old_edge_array);

            //Then insert new edge to subsequent block
            local_head_block_ptr->adj_vertex[0] = out_edge;
            clear_bit(&local_head_block_ptr->bitvector_64, 0);

            //analysis
            type64++;

        }
    }
    else {
        //8B case
        WeightedOutEdge out_edge = {edge.des, edge.weight, local_head_block_ptr->timestamp};
        auto local_head_block_ptr_8b = (AdjSubsequentBlockFive*)local_head_block_ptr;
        if (local_head_block_ptr_8b->bitvector_64 != 0) {
            //not full yet
            uint64_t temp_bv_rev = __builtin_bswap64(local_head_block_ptr_8b->bitvector_64.load());
            insert_index = __builtin_clzl(temp_bv_rev);
            type = local_head_block_ptr_8b->type;
            clear_bit(&local_head_block_ptr_8b->bitvector_64, insert_index);
            uint64_t temp;
            if (insert_index > ((1 << (type - 1)) * 4 - 1)) {
                //need to malloc new space;
                switch (type) {
                    case 1: {
                        auto new_block = (SpruceTransVer::AdjSubsequentBlockTwo*) malloc(
                                sizeof(AdjSubsequentBlockTwo));
                        memset(new_block, 0, sizeof(AdjSubsequentBlockTwo));
                        temp = (uint64_t)local_head_block_ptr_8b;
                        memcpy(new_block,
                               (SpruceTransVer::AdjSubsequentBlockOne*)(temp),
                               sizeof(SpruceTransVer::AdjSubsequentBlockOne));
                        new_block->adj_vertex[insert_index] = out_edge;
                        new_block->type.store(2);
                        local_head_block_ptr_8b->obsolete_flag = 1;
                        // add previous block to delete queue
                        ptr_block->ptr_to_buffer[index_in_64].store(reinterpret_cast<uint64_t>(new_block));
                        free((SpruceTransVer::AdjSubsequentBlockOne*)(temp));
                        local_head_block_ptr_8b = (SpruceTransVer::AdjSubsequentBlockFive*)new_block;
                        //analysis
                        type2++;
                        type1--;

                        break;
                    }
                    case 2: {
                        auto new_block = (SpruceTransVer::AdjSubsequentBlockThree*) malloc(
                                sizeof(AdjSubsequentBlockThree));
                        memset(new_block, 0, sizeof(AdjSubsequentBlockThree));
                        temp = (uint64_t)local_head_block_ptr_8b;
                        memcpy(new_block,
                               (SpruceTransVer::AdjSubsequentBlockTwo*)temp,
                               sizeof(SpruceTransVer::AdjSubsequentBlockTwo));
                        new_block->adj_vertex[insert_index] = out_edge;
                        new_block->type.store(3);
                        local_head_block_ptr_8b->obsolete_flag = 1;
                        // add previous block to delete queue
                        ptr_block->ptr_to_buffer[index_in_64].store(reinterpret_cast<uint64_t>(new_block));
                        free((SpruceTransVer::AdjSubsequentBlockTwo*)temp);
                        //analysis
                        type3++;
                        type2--;

                        break;
                    }
                    case 3: {
                        auto new_block = (SpruceTransVer::AdjSubsequentBlockFour*) malloc(
                                sizeof(AdjSubsequentBlockFour));
                        memset(new_block, 0, sizeof(AdjSubsequentBlockFour));
                        temp = (uint64_t)local_head_block_ptr_8b;
                        memcpy(new_block,
                               (SpruceTransVer::AdjSubsequentBlockThree*)temp,
                               sizeof(SpruceTransVer::AdjSubsequentBlockThree));
                        new_block->adj_vertex[insert_index] = out_edge;
                        new_block->type.store(4);
                        local_head_block_ptr_8b->obsolete_flag = 1;
                        // add previous block to delete queue
                        ptr_block->ptr_to_buffer[index_in_64].store(reinterpret_cast<uint64_t>(new_block));
                        free((SpruceTransVer::AdjSubsequentBlockThree*)temp);
                        //analysis
                        type4++;
                        type3--;

                        break;

                    }
                    case 4: {
                        auto new_block = (SpruceTransVer::AdjSubsequentBlockFive*) malloc(
                                sizeof(AdjSubsequentBlockFive));
                        temp = (uint64_t)local_head_block_ptr_8b;
                        memset(new_block, 0, sizeof(AdjSubsequentBlockFive));
                        // new_block->bitvector_64 = UINT64_MAX; //no subsequent block
                        memcpy(new_block,
                               (SpruceTransVer::AdjSubsequentBlockFour*)(temp),
                               sizeof(SpruceTransVer::AdjSubsequentBlockFour));
                        new_block->adj_vertex[insert_index] = out_edge;
                        new_block->type.store(5);
                        local_head_block_ptr_8b->obsolete_flag = 1;
                        // add previous block to delete queue
                        ptr_block->ptr_to_buffer[index_in_64].store(reinterpret_cast<uint64_t>(new_block));
                        free((SpruceTransVer::AdjSubsequentBlockFour*)temp);
                        //analysis
                        type5++;
                        type4--;

                        break;
                    }
                    default: {
                        std::cout << "Impossible type!" << std::endl;
                        break;
                    }

                }
            } else {
                local_head_block_ptr_8b->adj_vertex[insert_index] = out_edge;
            }
        } else {
            type = 5;
            // while ends, no spaces
            // Firstly check the block size,we assume that we only check invalid label in deletion
            uint64_t old_block_size, old_delete_num;
            uint64_t temp1 = (uint64_t)local_head_block_ptr_8b;
            uint64_t temp2 = local_head_block_ptr_8b->next_block.load();
            auto old_edge_array = (uint64_t*)temp2;
            if (!old_edge_array) {
                // does not exist, set initial size

                old_block_size = 0;
                old_delete_num = 0;
                // analysis
                edge_array_num++;
            }
            else {
                old_block_size = old_edge_array[0];
                old_delete_num = old_edge_array[1];
            }

            auto old_edges = (SpruceTransVer::WeightedOutEdge*)(old_edge_array + 2);

            int64_t old_delete_64 = 0;
            if (old_delete_num > 64) {
                old_delete_64 = old_delete_num/64;
            }

            // + 1 for block_size and delete_num (equals to weighted out edge)
            // resize block according to delete_num
            auto new_block_size = (64 + 1 + old_block_size - old_delete_64 * 64);
            auto new_edge_array = (uint64_t*) malloc(
                    sizeof(SpruceTransVer::WeightedOutEdge) * (new_block_size));
            memset(new_edge_array, 0, sizeof(SpruceTransVer::WeightedOutEdge) * (new_block_size));
//

            //sort the vertex using bubble sort
            for (int i = 0; i < 64; i++) {
                bool flag = false;
                for (int j = 0; j < 64 - i - 1; j++) {
                    if (local_head_block_ptr_8b->adj_vertex[j].des > local_head_block_ptr_8b->adj_vertex[j + 1].des) {
                        auto temp = local_head_block_ptr_8b->adj_vertex[j];
                        local_head_block_ptr_8b->adj_vertex[j] = local_head_block_ptr_8b->adj_vertex[j + 1];
                        local_head_block_ptr_8b->adj_vertex[j + 1] = temp;
                        flag = true;
                    }
                }
                if (!flag) {
                    break;
                }
            }

            //then use merge sort to place spaces in new block
            uint64_t k = 0;
            uint64_t start1 = 0, start2 = 0;
            uint64_t end1 = 64, end2 = old_block_size;

            auto new_edges = (SpruceTransVer::WeightedOutEdge*)(new_edge_array + 2);

            while (start1 < end1 && start2 < end2) {
                if (old_edges[start2].des == UINT64_MAX) {
                    // skip invalid data
                    start2++;
                    continue;
                }
                new_edges[k++] = local_head_block_ptr_8b->adj_vertex[start1].des < old_edges[start2].des ?
                                 local_head_block_ptr_8b->adj_vertex[start1++] : old_edges[start2++];
            }
            while (start1 < end1) {
                new_edges[k++] = local_head_block_ptr_8b->adj_vertex[start1++];
            }
            while (start2 < end2) {
                if (old_edges[start2].des == UINT64_MAX) {
                    // skip invalid data
                    start2++;
                    continue;
                }
                new_edges[k++] = old_edges[start2++];
            }
            while (k < new_block_size - 1) {
                // shift invalid values with size < 64
                new_edges[k++].des = UINT64_MAX;
            }


            //reset head block status
            local_head_block_ptr_8b->bitvector_64.store(UINT64_MAX);

            //reset subsequent block status
            memset(local_head_block_ptr_8b->adj_vertex, 0, sizeof(SpruceTransVer::WeightedOutEdge)*64);
            new_edge_array[0] = new_block_size - 1; // remember: block_size does not include first 2 uint32_t for information
            new_edge_array[1] = old_delete_num % 64;  // copy delete_num
            local_head_block_ptr_8b->next_block.store(reinterpret_cast<uint64_t>(new_edge_array));
            // Generally you need to use your own garbage collecting function in your system instead of simply free it
            if(old_edge_array) free(old_edge_array);

            //Then insert new edge to subsequent block
            local_head_block_ptr_8b->adj_vertex[0] = out_edge;
            clear_bit(&local_head_block_ptr_8b->bitvector_64, 0);

            //analysis
            type64++;

        }
    }
    ptr_block->buffer_locks[index_in_64]--;
    return true;
}

uint64_t SpruceTransVer::get_global_timestamp() {
    return global_timestamp++;
}

bool SpruceTransVer::get_neighbours_snapshot(SpruceTransVer &spruce, const uint64_t from_node_id,
                                    std::vector<WeightedOutEdgeSimple> &neighbours) {
    auto read_timestamp = SpruceTransVer::get_global_timestamp();

    //do not use locks
    SpruceTransVer::TopBlock* root;
    uint32_t restart_num = 0;
    auto hash_index = (uint32_t)(from_node_id >> 32);
    auto from_node_id_low = (uint16_t)from_node_id;
    auto from_node_id_high = (uint16_t)(from_node_id >> 16);
    restart:
    if (hash_index || spruce.fb_flag) {
        // if existed in hash table
        root = spruce.spruce_hash.get(hash_index);
        // if not existed
        if (!root) {
            return false;
        }
    }
    else {
        root = spruce.top_block;
    }
    neighbours.clear();
    uint64_t temp_ptr;
    if (!get_bit(&(root->bitmap_8kb), (uint32_t) from_node_id_high)) {
        return false;
    }
    temp_ptr = root->ptr_to_children[from_node_id_high].load();
    auto* middle_block_ptr = (SpruceTransVer::MiddleBlock*)temp_ptr;
    if (!middle_block_ptr) {
        return false;
    }
    if (!get_bit(&(middle_block_ptr->bitmap_8kb), (uint32_t)from_node_id_low)) {
        return false;
    }
    auto ptr_block_index = from_node_id_low / 64;
    auto auxiliary_ptr = reinterpret_cast<uint64_t*>(&middle_block_ptr->bitmap_8kb);
    uint64_t auxiliary_64 = auxiliary_ptr[ptr_block_index];
    uint64_t auxiliary_64_rev = __builtin_bswap64(auxiliary_64);
    auto ptr_num = __builtin_popcountl(auxiliary_64);
    auto index_in_64 = (uint32_t) (from_node_id_low % 64);

    SpruceTransVer::AdjSubsequentBlockFive4B* bottom_head_block;

    SpruceTransVer::PtrBlock* ptr_block = (SpruceTransVer::PtrBlock*)middle_block_ptr->ptr_to_children[ptr_block_index].load();
    bottom_head_block = (SpruceTransVer::AdjSubsequentBlockFive4B*)(ptr_block->ptr_to_buffer[index_in_64].load());

    if (!bottom_head_block) {
        return false;
    }

    SpruceTransVer::AdjSubsequentBlockFive4B* local_head_block_ptr = bottom_head_block;
    uint32_t get_index = 0;
    uint64_t temp_bitvector;

    recheck_bottom_lock:
    while (ptr_block->buffer_locks[index_in_64].load() != UNLOCKED);
    uint32_t timestamp_before = bottom_head_block->timestamp.load();
    // also recheck
    if (ptr_block->buffer_locks[index_in_64].load() == WRITE_LOCKED) {
        goto recheck_bottom_lock;
    }
    // allow reading from obsoleted block
    if (!(local_head_block_ptr->fb_flag_log_size >> 15)){
        //4B case
        if (local_head_block_ptr->bitvector_64 != UINT64_MAX) {
            //value existed
            uint64_t temp_bv_rev = __builtin_bswap64(local_head_block_ptr->bitvector_64.load());
            auto pre_bv = bottom_head_block->bitvector_64.load();
            auto type = local_head_block_ptr->type.load();
            if (type == 5) {
                //type = 5;
                //first check subsequent block
                for (int i = 0; i < 64; i++) {
                    if (!get_bit(&pre_bv, i)) {
                        if (local_head_block_ptr->adj_vertex[i].timestamp <= read_timestamp) {
                            neighbours.push_back(convert_to_8b_simple(local_head_block_ptr->adj_vertex[i]));
                        }
                    }
                }
                //then check edge array
                temp_ptr = local_head_block_ptr->next_block.load();
                auto edge_array = (uint32_t*)temp_ptr;
                uint32_t block_size, delete_space;
                if (!edge_array) {
                    block_size = 0;
                    delete_space = 0;
                }
                else {
                    block_size = edge_array[0];
                    delete_space = edge_array[1];
                }
                uint32_t new_index = neighbours.size(); // use new index as array index to improve efficiency (instead of push_back)
                neighbours.resize(neighbours.size() + block_size - delete_space);
                auto edges_ptr = (SpruceTransVer::WeightedOutEdge4B*)(edge_array + 2);

                for (int i = 0; i < block_size; i++) {
                    if (edges_ptr[i].des != UINT32_MAX) {
                        neighbours[new_index++] = convert_to_8b_simple(edges_ptr[i]);
                    }
                }


            } else {
                //type = 4;
                for (int i = 0; i < (2 << type) ; i++) {
                    if (!get_bit(&pre_bv, i)) {
                        neighbours.push_back(convert_to_8b_simple(local_head_block_ptr->adj_vertex[i]));
                    }
                }
            }
        } else {
            // check edge array
            if (!local_head_block_ptr->next_block.load()) {
                return false;
            }
            else {
                temp_ptr = local_head_block_ptr->next_block.load();
                auto edge_array = (uint32_t*)temp_ptr;
                uint32_t block_size, delete_space;
                block_size = edge_array[0];
                delete_space = edge_array[1];
                uint32_t new_index = neighbours.size(); // use new index as array index to improve efficiency (instead of push_back)
                neighbours.resize(neighbours.size() + block_size - delete_space);
                auto edges_ptr = (SpruceTransVer::WeightedOutEdge4B*)(edge_array + 2);

                for (int i = 0; i < block_size; i++) {
                    if (edges_ptr[i].des != UINT32_MAX) {
                        neighbours[new_index++] = convert_to_8b_simple(edges_ptr[i]);
                    }
                }
            }
        }
    }
    else {
        // 8B case
        auto local_head_block_ptr_8b = (SpruceTransVer::AdjSubsequentBlockFive*)local_head_block_ptr;
        if (local_head_block_ptr_8b->bitvector_64 != UINT64_MAX) {
            //value existed
            uint64_t temp_bv_rev = __builtin_bswap64(local_head_block_ptr_8b->bitvector_64.load());
            auto pre_bv = bottom_head_block->bitvector_64.load();
            auto type = local_head_block_ptr_8b->type.load();
            if (type == 5) {
                //type = 5;
                //first check subsequent block
                for (int i = 0; i < 64; i++) {
                    if (!get_bit(&pre_bv, i)) {
                        if (local_head_block_ptr->adj_vertex[i].timestamp <= read_timestamp) {
                            neighbours.push_back({local_head_block_ptr_8b->adj_vertex[i].des, local_head_block_ptr_8b->adj_vertex[i].weight});
                        }
                    }
                }
                //then check edge array
                temp_ptr = local_head_block_ptr_8b->next_block.load();
                auto edge_array = (uint64_t*)temp_ptr;
                uint64_t block_size, delete_space;
                if (!edge_array) {
                    block_size = 0;
                    delete_space = 0;
                }
                else {
                    block_size = edge_array[0];
                    delete_space = edge_array[1];
                }
                uint64_t new_index = neighbours.size(); // use new index as array index to improve efficiency (instead of push_back)
                neighbours.resize(neighbours.size() + block_size - delete_space);
                auto edges_ptr = (SpruceTransVer::WeightedOutEdge*)(edge_array + 2);

                for (int i = 0; i < block_size; i++) {
                    if (edges_ptr[i].des != UINT64_MAX) {
                        neighbours[new_index++] = {edges_ptr[i].des, edges_ptr[i].weight};
                    }
                }


            } else {
                //type = 4;
                for (int i = 0; i < (2 << type) ; i++) {
                    if (!get_bit(&pre_bv, i)) {
                        neighbours.push_back({local_head_block_ptr_8b->adj_vertex[i].des, local_head_block_ptr_8b->adj_vertex[i].weight});
                    }
                }
            }
        } else {
            // check edge array
            if (!local_head_block_ptr_8b->next_block.load()) {
                return false;
            }
            else {
                temp_ptr = local_head_block_ptr_8b->next_block.load();
                auto edge_array = (uint64_t*)temp_ptr;
                uint64_t block_size, delete_space;
                block_size = edge_array[0];
                delete_space = edge_array[1];
                uint64_t new_index = neighbours.size(); // use new index as array index to improve efficiency (instead of push_back)
                neighbours.resize(neighbours.size() + block_size - delete_space);
                auto edges_ptr = (SpruceTransVer::WeightedOutEdge*)(edge_array + 2);

                for (int i = 0; i < block_size; i++) {
                    if (edges_ptr[i].des != UINT64_MAX) {
                        neighbours[new_index++] = {edges_ptr[i].des, edges_ptr[i].weight};
                    }
                }
            }
        }
    }

    uint32_t timestamp_after = bottom_head_block->timestamp.load();
    if (timestamp_before != timestamp_after) {
        if (restart_num < RESTART_THRESHOLD){
            restart_num++;
            goto restart;
        }
        else {
            // read bottom head block exclusively
            return get_neighbours_exclusively(spruce, from_node_id, neighbours);
        }
    }

    restart_snapshot:
    timestamp_before = bottom_head_block->timestamp.load();
    // check version
    if (local_head_block_ptr->timestamp > read_timestamp) {
        // need to read previous version



        // then restore updated edges/deleted edges
        auto log_block_ptr = (SpruceTransVer::LogBlock*)local_head_block_ptr->ptr_to_log_block.load();
        auto outer_offset = (local_head_block_ptr->fb_flag_log_size & 0x7FFF) / 64;
        for (int i = 0; i < outer_offset - 1; i++) {
            if (log_block_ptr->versioned_edges[BLOCK_SIZE-1].timestamp && TIMESTAMP_MASK <= read_timestamp) {
                continue;
            }
            for (int j = BLOCK_SIZE - 1; j >=0; j--) {
                if (log_block_ptr->versioned_edges[j].timestamp && TIMESTAMP_MASK > read_timestamp) {
                    // cases
                    auto current_type = (uint32_t) (log_block_ptr->versioned_edges[j].timestamp >> 56);
                    if (current_type == DELETE_TYPE) {
                        neighbours.push_back({log_block_ptr->versioned_edges[j].ver.des,log_block_ptr->versioned_edges[j].ver.weight});
                    }
                    else if (current_type == UPDATE_TYPE) {
                        if (neighbours.size() <= 64) {
                            for (int k = 0; k < neighbours.size(); k++) {
                                if (neighbours[k].des == log_block_ptr->versioned_edges[j].ver.des) {
                                    neighbours[k].weight = log_block_ptr->versioned_edges[j].ver.weight;
                                }
                            }
                        }
                        else {
                            for (int k = 0; k < 64; k++) {
                                if (neighbours[k].des == log_block_ptr->versioned_edges[j].ver.des) {
                                    neighbours[k].weight = log_block_ptr->versioned_edges[j].ver.weight;
                                }
                            }
                            // binary search
                            int64_t low = 64;
                            int64_t high = neighbours.size() - 1;
                            int64_t mid;
                            while (high >= low) {
                                mid = (low + high) / 2;
                                if(log_block_ptr->versioned_edges[j].ver.des < neighbours[mid].des){
                                    high = mid - 1;
                                }
                                else if (log_block_ptr->versioned_edges[j].ver.des == neighbours[mid].des) {
                                    neighbours[mid].weight = log_block_ptr->versioned_edges[j].ver.weight;
                                    break;
                                }
                                else {
                                    low = mid + 1;
                                }
                            }
                        }
                    }
                }
                else {
                    break;
                }
            }
            log_block_ptr = log_block_ptr->next_block;
        }
    }
    timestamp_after = bottom_head_block->timestamp.load();
    if (timestamp_before != timestamp_after) {
        goto restart_snapshot;
    }

    return true;
}
