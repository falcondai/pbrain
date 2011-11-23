/*     Version  7-NOV-1990 18:27  */

/*                                                                           *
 *                          Copyright (c) 1990:                              *
 *                                                                           *
 *                       The University of Chicago,                          *
 *                     Department of Radiation Oncology;                     *
 *                      Michael Reese / U of Chicago                         *
 *                      Center for Radiation Therapy.                        *
 *                                                                           *
 *                         All rights reserved.                              *
 *                                                                           *
 *      For information contact:                                             *
 *      C. Pelizzari                                                         *
 *      5841 S. Maryland Avenue, Box 442                                     *
 *      Chicago, IL 60637                                                    *
 *      Tel. (312) 702-6841                                                  *
 *      Email chuck@rover.UChicago.EDU                                       *
 *                                                                           *
 *                                                                           */

#include <stdio.h>
#include <curses.h>

#ifdef THINKC
#include <unix.h>
#include <strings.h>
#include <storage.h>
#endif

#ifndef MAX
#define MAX(a,b) ( (a) > (b) ? (a) : (b))
#endif

#include <string.h>
#include <malloc.h>

int GetRealList(numlist, NullValue)

  float **numlist, NullValue;

{

  char inbuf[80];

  int ParseRealList();
  char format[127], *ptr;
  int i;

  ptr = format;
  strcpy(format, "%[^");
  ptr += strlen(format);

  for (i = 1; i <= 31; i++) {

    *ptr = (char) i;
    ptr++;  
  } /*   for   */
  strcpy(ptr, "]");
  
  inbuf[0] = '\0';
  
   *numlist = NULL;  

   scanf(format, inbuf);
   i = getchar();

    if (!inbuf[0]) return(0);

    return (ParseRealList(numlist, NullValue, inbuf));

} /*  */


int ParseRealList(numlist, NullValue, inbuf)

  float **numlist, NullValue;
  char *inbuf;

{
  char token[80], *remainder;
  int i, buflen, toklen;
  int  howmany;

  int ParseRealList();

    buflen = strlen(inbuf);
    *numlist = NULL;


    if (!buflen) return(0);

/*    printf("string length = %d\n", buflen);*/
    howmany = 0;
    
    *numlist=(float *) malloc(sizeof(float));
    remainder = inbuf;
    if (*remainder == ','){
         *numlist = (float *) realloc(*numlist,(howmany + 1) 
          * sizeof(float));
           (*numlist)[howmany++] = NullValue; }/*flag skipped fields */

    while (buflen > 0) {
    
      while (*remainder == ',' || *remainder == ' ') 
        { if (*remainder == ',' && *(remainder - 1) == ',') {
         *numlist = (float *) realloc(*numlist,(howmany + 1) 
          * sizeof(float));
           (*numlist)[howmany++] = NullValue; }/*flag skipped fields */
            remainder++; buflen--;}
          
      toklen = strcspn(remainder, ", ");
/*      printf("%s %d\n", remainder, toklen);*/

      strncpy(token, remainder, toklen);
      token[toklen] = 0;

/*      printf("this token %s \n",token);*/

      remainder += (toklen + 1);
      buflen -= toklen + 1;
/*      printf("length now %d\n", buflen);*/
 

      if (strspn(token, "-+.0123456789") == toklen)  {

        *numlist = (float *) realloc(*numlist,(howmany + 1) * sizeof(float));
	    (*numlist)[howmany] = NullValue;
        sscanf(token, "%f", &((*numlist)[howmany++]));  /* single number  */

      } /* if (strspn(token, "-+. 0123456789") == toklen)   */


             

    }  /* while buflen...*/
    
    while ((*numlist)[howmany - 1] == NullValue  && howmany > 0) howmany--; 

/*    printf("%3d numbers:\n",howmany);
    for (i = 0; i < howmany; printf("%f ",(*numlist)[i++])); printf("\n");*/

    if (!howmany) {

        free(*numlist);
        remainder= (char*) malloc (strlen(inbuf) *sizeof(char));
	strcpy(remainder, inbuf);
        *numlist = (float *) remainder;

      } /*        */
    return(howmany);
  
}




int GetRealListWithPrompt(numlist, bogus, prompt)

float **numlist, bogus;
char *prompt;

{
  printf("%s",prompt);
  return( GetRealList(numlist, bogus) );
}



int WGetRealListWithPrompt(win, numlist, bogus, prompt)

WINDOW *win;
float **numlist;
float bogus;
char *prompt;

{

  char tbuf[80];
  char format[127], *ptr;
  int i;

  ptr = format;
  strcpy(format, "%[^");
  ptr += strlen(format);

  for (i = 1; i < 32; i++) {
    if (i == 3 || i == 26) continue; 
    *ptr = (char) i;
    ptr++;  
  } /*      */
  strcpy(ptr, "]s");

  tbuf[0] = 0;
  wclrtoeol(win);
  if (prompt) wprintw(win, "%s", prompt);
  wrefresh(win);

  wscanw(win, format, tbuf);

  return (ParseRealList(numlist, bogus, tbuf));


} /*  */

#ifdef DTEST

#define BOGUS 999.9

main ()

{
  int i,n;
  float *RealList;


  while (1) {



    n = GetRealListWithPrompt(&RealList, BOGUS, "Enter list:");
    for (i = 0; i < n; i++) printf( "%f ",RealList[i]); 
    printf("\n");

    if (n > 0 && RealList[0] == -1.0) break;
    if (RealList) free(RealList);

  } /* while (1)  */


} /*  */



#endif




