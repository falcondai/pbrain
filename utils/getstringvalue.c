
/*     Version 10-SEP-1990 15:45  */

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

#ifndef THINK_C
#ifdef ULTRIX
#include <cursesX.h>
#else
#include <curses.h>
#endif
#endif

#include <string.h>
#include <malloc.h>


char *GetaString(String, Prompt)

char **String;
char *Prompt;

{

  char tbuf[80];
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

  tbuf[0] = 0;

  if (Prompt) printf("%s", Prompt);
  scanf(format, tbuf);
  i=getchar();

  i = strlen(tbuf);
  if (!i) return (0);

  *String = (char*) malloc ((i + 1) * sizeof(char));
  strcpy(*String, tbuf);

  return (*String);   

} /*  */


char *WGetaString(win, String, Prompt)

WINDOW *win;
char **String;
char *Prompt;

{

  char tbuf[80];
  char format[127], *ptr;
  int i;

  ptr = format;
  strcpy(format, "%[^");
  ptr += strlen(format);

  for (i = 1; i < 32; i++) {
    if (i == 27) continue; 
    *ptr = (char) i;
    ptr++;  
  } /*      */
  strcpy(ptr, "]s");

  tbuf[0] = 0;
  if (Prompt) wprintw(win, "%s", Prompt);
  wclrtoeol(win);
  noecho();
  wrefresh(win);
  wscanw(win, "%s", tbuf);

  echo();


  i = strlen(tbuf);
  if (!i) return (0);

  *String = (char*) malloc ((i + 1) * sizeof(char));
  strcpy(*String, tbuf);

  return (*String);   

} /*  */

int GetValue(Value, Format, Prompt)

int *Value;
char *Format, *Prompt;

{

  char tbuf[80];
  char format[127], *ptr;
  int i;

  ptr = format;
  strcpy(format, "%[^");
  ptr += strlen(format);

  for (i = 1; i < 32; i++) {

    *ptr = ( char) i;
    ptr++;  
  } /*   for   */
  strcpy(ptr, "]s");

  tbuf[0] = 0;
  if (Prompt) printf("%s", Prompt);
  scanf(format, tbuf);
  i=getchar();

  i = strlen(tbuf);
  if (!i) return (0);

  sscanf(tbuf, Format, Value);

  return (1);   

} /*  */

#ifdef DTEST

main ()

{
int i;

  char *string;


  while ( (string = GetaString( &string, "enter string") )){

    printf( "%s\n", string);
    free(string);
  } /* */

} /*  */
#endif

