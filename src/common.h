#ifndef COMMON_H
#define COMMON_H

#define FUSIONSCAN_VER "0.8.0"

#define _DEBUG true

typedef long int64;
typedef unsigned long uint64;

typedef int int32;
typedef unsigned int uint32;

typedef short int16;
typedef unsigned short uint16;

typedef char int8;
typedef unsigned char uint8;


#pragma pack(2) 
// if contig is -1, means this is a dupe entry, and position will be the position in the dupList
struct GenePos{
    short contig;
    int position;
};
#pragma pack() 

// the limit of the queue to store the packs
// error may happen if it generates more packs than this number
static const int PACK_NUM_LIMIT  = 5000000;

// how many reads one pack has
static const int PACK_SIZE = 1000;

// if one pack is produced, but not consumed, it will be kept in the memory
// this number limit the number of in memory packs
// if the number of in memory packs is full, the producer thread should sleep
static const int PACK_IN_MEM_LIMIT = 100;

// the key dup in normal level will be kept, in high level will be skipped
static const int DUPE_NORMAL_LEVEL = -1;
static const int DUPE_HIGH_LEVEL = -2;


#endif /* COMMON_H */
