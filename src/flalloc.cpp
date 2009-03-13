#include "flalloc.hpp"

/*********************************************************************************/
/*                                                                               */
/* Dynamic Array Allocation Routines                                             */
/*                                                                               */
/*********************************************************************************/

/**********************************************/
/* char                          */
/**********************************************/

/**************/
/* allocation */
/**************/

bool flallocArray(LPCHAR * pv, long nl, long nh)
   {
   if (nl > nh) return FALSE;

   *pv=(LPCHAR)MALLOC((nh-nl+1)*sizeof(char));
   
   if (!*pv) 
      return FALSE;
   
   *pv-=nl;

   return TRUE;
   }

bool flallocArray(LP2CHAR * pm, long nrl, long nrh, long ncl, long nch)
   {  
   if (nrl > nrh || ncl > nch) return FALSE;

   *pm=(char**) MALLOC((nrh-nrl+1)*sizeof(LPCHAR));
   
   if (!*pm) 
     return FALSE;
   
   *pm -= nrl;

   for(long i=nrl;i<=nrh;i++) 
      {
      (*pm)[i]=(LPCHAR) MALLOC((nch-ncl+1)*sizeof(char));
      
      if (!(*pm)[i]) 
         return FALSE;
      
      (*pm)[i] -= ncl;
      }                                                           

   return TRUE;
   }                                                              

bool flallocArray(LP2CHAR * pv, long nl, long nh)
   {
   if (nl > nh) return FALSE;
   
   *pv=(LP2CHAR)MALLOC((nh-nl+1)*sizeof(LPCHAR));
   
   if (!*pv) 
      return FALSE;
   
   *pv-=nl;

   return TRUE;
   }

/**************/
/* free       */
/**************/

void free_flallocArray(LPCHAR v, long nl, long nh)
   {  
   if (nh >= nl)
      FREE((char * ) (v+nl));
   }

void free_flallocArray(LP2CHAR m, long nrl, long nrh, long ncl, long nch)
   {
   for(long i=nrh;i>=nrl;i--)
       FREE((char * ) (m[i]+ncl));
         
   if (nrh >= nrl)
      FREE((char * ) (m+nrl));
   }                                                               

void free_flallocArray(LP2CHAR v, long nl, long nh)
   {  
   if (nh >= nl)
      FREE((char * ) (v+nl));
   }


/**********************************************/
/* shorts                        */
/**********************************************/

/**************/
/* allocation */
/**************/

bool flallocArray(LP2SHORT * pm, long nrl, long nrh, long ncl, long nch)
   {  
   if (nrl > nrh || ncl > nch) return FALSE;

   *pm=(LP2SHORT) MALLOC((nrh-nrl+1)*sizeof(LPSHORT));
   
   if (!*pm) 
      return FALSE;
   
   *pm -= nrl;

   for(long i=nrl;i<=nrh;i++) 
      {
      (*pm)[i]=(LPSHORT) MALLOC((nch-ncl+1)*sizeof(long));
 
      if (!(*pm)[i]) 
         return FALSE;
 
      (*pm)[i] -= ncl;
      }
   
   return TRUE;
   }       

bool flallocArray(LP3SHORT * pv, long nl, long nh)
   {
   if (nl > nh) return FALSE;
   
   *pv=(LP3SHORT)MALLOC((nh-nl+1)*sizeof(LP2SHORT));
   
   if (!*pv) 
      return FALSE;
   
   *pv-=nl;

   return TRUE;
   }       


/**************/
/* free       */
/**************/

void free_flallocArray(LP3SHORT v, long nl, long nh)
   {  
   if (nh >= nl)
      FREE((char * ) (v+nl));
   }

void free_flallocArray(LP2SHORT m, long nrl, long nrh, long ncl, long nch)
   {
   for(long i=nrh;i>=nrl;i--)
         FREE((char * ) (m[i]+ncl));
      
   if (nrh >= nrl)
      FREE((char * ) (m+nrl));
   }
 

/**********************************************/
/* ints                          */
/**********************************************/

/**************/
/* allocation */
/**************/

bool flallocArray(LPINT * pv, long nl, long nh)
   {
   if (nl > nh) return FALSE;

   *pv=(LPINT)MALLOC((nh-nl+1)*sizeof(int));
   
   if (!*pv) 
      return FALSE;
   
   *pv-=nl;

   return TRUE;
   }


/**************/
/* free       */
/**************/

void free_flallocArray(LPINT v, long nl, long nh)
   {  
   if (nh >= nl)
      FREE((char * ) (v+nl));
   }




/**********************************************/
/* ints                          */
/**********************************************/

/**************/
/* allocation */
/**************/
bool flallocArray(LPSHORT * pv, long nl, long nh)
   {
   if (nl > nh) return FALSE;

   *pv=(LPSHORT)MALLOC((nh-nl+1)*sizeof(short));
   
   if (!*pv) 
      return FALSE;
   
   *pv-=nl;

   return TRUE;
   }

/**************/
/* free       */
/**************/

void free_flallocArray(LPSHORT v, long nl, long nh)
   {  
   if (nh >= nl)
      FREE((char * ) (v+nl));
   }
          

 
/**********************************************/
/* doubles                                    */
/**********************************************/

/**************/
/* allocation */
/**************/

bool flallocArray(LPDOUBLE * pv, long nl, long nh)
   {
   if (nl > nh) return FALSE;

   *pv=(LPDOUBLE)MALLOC((nh-nl+1)*sizeof(double));
   
   if (!*pv) 
      return FALSE;
   
   *pv-=nl;

   return TRUE;
   }

bool flallocArray(LP2DOUBLE * pv, long nl, long nh)
   {
   if (nl > nh) return FALSE;
   
	*pv=(LP2DOUBLE)MALLOC((nh-nl+1)*sizeof(LPDOUBLE));
   
   if (!*pv) 
      return FALSE;
   
   *pv-=nl;

   return TRUE;
   }

bool flallocArray(LP3DOUBLE * pv, long nl, long nh)
   {
   if (nl > nh) return FALSE;
   
   *pv=(LP3DOUBLE)MALLOC((nh-nl+1)*sizeof(LP2DOUBLE));
   
   if (!*pv) 
      return FALSE;
   
   *pv-=nl;

   return TRUE;
   }       

bool flallocArray(LP2DOUBLE * pm, long nrl, long nrh, long ncl, long nch)
   {  
   if (nrl > nrh || ncl > nch) return FALSE;

   *pm=(double**) MALLOC((nrh-nrl+1)*sizeof(LPDOUBLE));
   //*pm=(LP2DOUBLE) MALLOC((nrh-nrl+1)*sizeof(LPDOUBLE));
  
   if (!*pm) 
     return FALSE;
   
   *pm -= nrl;

   for(long i=nrl;i<=nrh;i++) 
      {
      (*pm)[i]=(LPDOUBLE) MALLOC((nch-ncl+1)*sizeof(double));
      
      if (!(*pm)[i]) 
         return FALSE;
      
      (*pm)[i] -= ncl;
      }                                                           

   return TRUE;
   }                                                              

bool flallocArray(LP3DOUBLE * pm, long nrl, long nrh, long ncl, long nch)
   {  
   if (nrl > nrh || ncl > nch) return FALSE;

   *pm=(LP3DOUBLE) MALLOC((nrh-nrl+1)*sizeof(LPDOUBLE));
   
   if (!*pm) 
      return FALSE;
   
   *pm -= nrl;

   for(long i=nrl;i<=nrh;i++) 
      {
      (*pm)[i]=(LP2DOUBLE) MALLOC((nch-ncl+1)*sizeof(double));
      
      if (!(*pm)[i]) 
         return FALSE;
      
      (*pm)[i] -= ncl;
      }                                                           

   return TRUE;
   }                                                              


/**************/
/* free       */
/**************/

void free_flallocArray(LPDOUBLE v, long nl, long nh)
   {  
   if (nh >= nl)
      FREE((char * ) (v+nl));
   }
           
void free_flallocArray(LP2DOUBLE v, long nl, long nh)
   {  
   if (nh >= nl)
      FREE((char * ) (v+nl));
   }

void free_flallocArray(LP3DOUBLE v, long nl, long nh)
   {  
   if (nh >= nl)
      FREE((char * ) (v+nl));
   }

void free_flallocArray(LP2DOUBLE m, long nrl, long nrh, long ncl, long nch)
   {
   for(long i=nrh;i>=nrl;i--)
       FREE((char * ) (m[i]+ncl));
         
   if (nrh >= nrl)
      FREE((char * ) (m+nrl));
   }                                                               
                                                               
void free_flallocArray(LP3DOUBLE m, long nrl, long nrh, long ncl, long nch)
   {
   for(long i=nrh;i>=nrl;i--)
       FREE((char * ) (m[i]+ncl));
         
   if (nrh >= nrl)
      FREE((char * ) (m+nrl));
   }
   
   
