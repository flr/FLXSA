#ifndef _INC_fl__types
#define _INC_fl__types

#define FiFiErrNull      0 
#define FiFiErrInput     2
#define FiFiErrInput2    3 

#define FALSE 0
#define TRUE  1

#ifndef max
#define max(a,b)  (((a) > (b)) ? (a) : (b))
#endif

#ifndef min
#define min(a,b)  (((a) < (b)) ? (a) : (b))
#endif

#define SQR(a)   ((a)*(a))

typedef enum tagConstQuantity 
	{
   FLcBiomass = 1,
   FLcNumber  = 2
	} FLConstQuantity;

typedef enum tagConstSelectivity 
	{
   FLcKnown        = 1,
   FLcPartialCatch = 2,
   FLcCPUEbyAge    = 3
	} FLConstSelectivity;

typedef char      * LPCHAR;
typedef LPCHAR    * LP2CHAR;
typedef LP2CHAR   * LP3CHAR;
typedef LP3CHAR   * LP4CHAR;

typedef short     * LPSHORT;
typedef LPSHORT   * LP2SHORT;
typedef LP2SHORT  * LP3SHORT;
typedef LP3SHORT  * LP4SHORT;

typedef int       * LPINT;
typedef LPINT     * LP2INT;
typedef LP2INT    * LP3INT;
typedef LP3INT    * LP4INT;

typedef long      * LPLONG;
typedef LPLONG    * LP2LONG;
typedef LP2LONG   * LP3LONG;
typedef LP3LONG   * LP4LONG;

typedef float     * LPFLOAT;  
typedef LPFLOAT   * LP2FLOAT;  
typedef LP2FLOAT  * LP3FLOAT;  
typedef LP3FLOAT  * LP4FLOAT;  

typedef double    * LPDOUBLE;
typedef LPDOUBLE  * LP2DOUBLE;
typedef LP2DOUBLE * LP3DOUBLE;
typedef LP3DOUBLE * LP4DOUBLE;
typedef LP4DOUBLE * LP5DOUBLE;
typedef LP5DOUBLE * LP6DOUBLE;
typedef LP6DOUBLE * LP7DOUBLE;

typedef char      * LPSTR;
typedef LPSTR     * LP2STR;
typedef LP2STR    * LP3STR;
typedef LP3STR    * LP4STR; 


//FLTypeXSAControl
struct XSAControl
   {
   double MSSTol,
          MinNSE,
          FShkSE,
          FPreSpwn,
          MPreSpwn;      

   short TSRange,
         TSPower,
         LRcrtAge,
         ConQAge,
         Shk2N,
         Shk2F,
         Shk2FYr,         
         Shk2FAge,
         MaxIters,
         PlusGroup,          
         VPA,
         SmoothQ,
         Shk2FTarget,
         NSEWts,
         ReturnType,
         TuningWindow;
};
#endif /* _INC_fl__types */



