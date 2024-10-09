// function to terminate code execution
// these functions are defined in main.cpp

#ifndef COH_TOPLEVEL
#include <sstream>
extern std::ostringstream  message;
#endif

/**************************************/
/*      main.cpp                      */
/**************************************/
int     cohTerminateCode      (std::string);
void    cohNotice             (std::string);
