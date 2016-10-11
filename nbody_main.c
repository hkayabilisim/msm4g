/** @mainpage
 *
 * This is the introduction.
 *
 */

#include <stdio.h>
#include "nbody_lib.h"


/** 
 * Main entrance to the code.
 */
int main()
{
    /*
    BodyList *list = NULL;
    int n;
    
    nbody_read_bodies("data/eight.ini", &n, &list);
    
  
    nbody_printbodylist(list);
    
    
    nbody_freebodylist(&list);
    */
    
    nbody_linkedlist_test1();

    
    return 0;
}
