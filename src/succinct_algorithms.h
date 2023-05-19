//
// Created by sjf on 3/26/2022.
//

#ifndef GRAPHINDEX_SUCCINCT_ALGORITHMS_H
#define GRAPHINDEX_SUCCINCT_ALGORITHMS_H

#include "header.h"

inline bool get_bit(void* bitmap, uint32_t index) {
    return *((uint8_t *)bitmap + index / 8) & (1 << (7 - index % 8));
}

inline void set_bit(void* bitmap, uint32_t index) {
    *((uint8_t *)bitmap + index / 8) |= (1 << (7 - index % 8));
}

inline void clear_bit(void* bitmap, uint32_t index) {
     *((uint8_t *)bitmap + index / 8) &= (~(1 << (7 - index % 8)));
}

inline uint32_t Rank_64bit (bool flag, uint32_t pos, uint64_t bitvector){
    if (pos == 0) {
        return 0;
    }
    //clear bits
    bitvector >>= 64 - pos;
    bitvector <<= 64 - pos;
    uint64_t x = bitvector;
    //rank by Hemming Weight
    const uint64_t m1 = 0x5555555555555555; //binary: 0101...
    const uint64_t m2 = 0x3333333333333333; //binary: 00110011..
    const uint64_t m4 = 0x0f0f0f0f0f0f0f0f; //binary:  4 zeros,  4 ones ...
    const uint64_t h01 = 0x0101010101010101; //the sum of 256 to the power of 0,1,2,3...
    x -= (x >> 1) & m1;             //put count of each 2 bits into those 2 bits
    x = (x & m2) + ((x >> 2) & m2); //put count of each 4 bits into those 4 bits
    x = (x + (x >> 4)) & m4;        //put count of each 8 bits into those 8 bits
    uint64_t num_1 = (x * h01) >> 56;  //returns left 8 bits of x + (x<<8) + (x<<16) + (x<<24) + ...
    if (flag == 1) {
        return (uint32_t) num_1;
    } else if (flag == 0) {
        return (uint32_t) (pos + 1 - num_1);
    } else {
        printf("Error in rank!\n");
        return 0;
    }
}

#endif //GRAPHINDEX_SUCCINCT_ALGORITHMS_H
