#ifndef MEMORY_USAGE_H
#define MEMORY_USAGE_H

#include "sys/types.h"
#include "sys/resource.h"


#ifdef __APPLE__
    #define UNIT 1
#else
    #define UNIT 0
#endif

#include "stdlib.h"
#include "stdio.h"
#include "math.h"
#include "string.h"
#include <iostream>

class MemoryUsage{

private:
    int who = RUSAGE_SELF;
    struct rusage usage;
    int ret;


public:

    inline float getValue_in_KB(bool output){

         ret=getrusage(who,&usage);
         float memSize = usage.ru_maxrss/pow(1024.0,UNIT);
         if(output)  std::cout << memSize << " KB\n";

         return memSize;
    }

    inline float getValue_in_MB(bool output){

         ret=getrusage(who,&usage);
         float memSize = usage.ru_maxrss/pow(1024.0,UNIT+1);
         if(output) std::cout << "Memory Usage: " << memSize << " MB\n";

         return memSize;
    }

    inline float getValue_in_GB(bool output){

         ret=getrusage(who,&usage);
         float memSize = usage.ru_maxrss/pow(1024.0,UNIT+1);
         if(output) std::cout << "Memory Usage: " << memSize << " GB\n";

         return memSize;
    }
};



#endif // MEMORY_USAGE_H