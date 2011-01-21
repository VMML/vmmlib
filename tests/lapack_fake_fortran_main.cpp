
#include <iostream>
#include <cstdlib>

extern void missing_MAIN__(void) 
{ 
    std::cout << "missing fortrain MAIN__ function." << std::endl;
    abort(); 
} 

#define MAIN__ missing_MAIN__
