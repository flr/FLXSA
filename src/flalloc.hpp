#ifndef _INC_flalloc
#define _INC_flalloc
                    
#include <malloc.h>  
#include "fl__types.hpp"

#define MALLOC   malloc
#define FREE     free

bool flallocArray(     LPCHAR  *, long, long); 
bool flallocArray(     LP2CHAR *, long, long);
bool flallocArray(     LP2CHAR *, long, long, long, long);

void free_flallocArray(LPCHAR   , long, long); 
void free_flallocArray(LP2CHAR  , long, long); 
void free_flallocArray(LP2CHAR  , long, long, long, long); 


bool flallocArray(     LP3CHAR *, long, long);
void free_flallocArray(LP3CHAR  , long, long); 


bool flallocArray(     LPINT  *, long, long); 
void free_flallocArray(LPINT   , long, long); 


bool flallocArray(     LPSHORT  *, long, long); 
bool flallocArray(     LP3SHORT *, long, long);
bool flallocArray(     LP2SHORT *, long, long, long, long);

void free_flallocArray(LPSHORT  , long, long); 
void free_flallocArray(LP3SHORT , long, long); 
void free_flallocArray(LP2SHORT , long, long, long, long); 


bool flallocArray(LPDOUBLE  *, long, long); 
bool flallocArray(LP2DOUBLE *, long, long);
bool flallocArray(LP3DOUBLE *, long, long);
bool flallocArray(LP2DOUBLE *, long, long, long, long);
bool flallocArray(LP3DOUBLE *, long, long, long, long);

void free_flallocArray(LPDOUBLE , long, long); 
void free_flallocArray(LP2DOUBLE, long, long);
void free_flallocArray(LP3DOUBLE, long, long);
void free_flallocArray(LP2DOUBLE, long, long, long, long);
void free_flallocArray(LP3DOUBLE, long, long, long, long);

#endif /* _INC_flalloc */

