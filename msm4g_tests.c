/** @file msm4g_tests.c
 * @brief The test functions to check the functionality of the MSM4G package.
 *
 * This is probably the best place one can learn how to use the package.
 * These test functions are used to check the outcomes of the functions
 * by comparing with expected results.
 */
#include "msm4g_lib.h"
#include "msm4g_tests.h"

void msm4g_unit_test_all()
{
    Boolean status;
    int index, numberoftests,failed;
    typedef Boolean (testFunctionType)();
    testFunctionType *testFunction;
    LinkedList *list = msm4g_linkedlist_new();
    
    msm4g_linkedlist_add(list, (void *)&msm4g_unit_test_1);
    msm4g_linkedlist_add(list, (void *)&msm4g_unit_test_2);
    msm4g_linkedlist_add(list, (void *)&msm4g_unit_test_3);
    msm4g_linkedlist_add(list, (void *)&msm4g_unit_test_4);
    msm4g_linkedlist_add(list, (void *)&msm4g_unit_test_5);
    msm4g_linkedlist_add(list, (void *)&msm4g_unit_test_6);
    
    numberoftests = msm4g_linkedlist_size(list);

    failed=0;
    for (index = 0; index < numberoftests ; index++)
    {
        testFunction = (testFunctionType *)msm4g_linkedlist_get(list, index);
        status = testFunction();
        if (status == false)
        {
            failed++;
            fprintf(stderr,"Unit test %3d: %s\n",index+1,"failed");
        }
    }
    printf("%d of %d passed, %d failed\n",numberoftests-failed,numberoftests,failed);
    
    msm4g_linkedlist_destroy(list);

}

Boolean msm4g_unit_test_1()
{
    LinkedList *list;
    Body x[10];
    int i;
    
    list = msm4g_linkedlist_new();
    
    for (i=0;i<10;i++)
    {
        msm4g_body_reset(x+i);
        msm4g_linkedlist_add(list,&x[i]);
    }
    if (msm4g_linkedlist_size(list) != 10)
    {
        msm4g_linkedlist_destroy(list);
        return false;
    }
    msm4g_linkedlist_destroy(list);
    return true;
}

Boolean msm4g_unit_test_2()
{
    LinkedList *list;
    LinkedListElement *curr;
    Body x[10];
    Body *y;
    Boolean status ;
    int i;
    
    status = true;
    list = msm4g_linkedlist_new();
    
    for (i=0;i<10;i++)
    {
        x[i].m = i;
        msm4g_linkedlist_add(list,&x[i]);
    }
    
    /**< Iterating from tail to head */
    i=9;
    curr=list->tail;
    while (curr != NULL)
    {
        y =  (Body *)curr->data;
        if (fabs(i - (y->m)) > DBL_EPSILON)
        {
            status = false;
            break;
        }
        curr = curr->prev;
        i--;
    }
    msm4g_linkedlist_destroy(list);
    
    return status;
}

Boolean msm4g_unit_test_3()
{
    LinkedList *list;
    LinkedListElement *curr;
    Body x[10];
    Body *y;
    Boolean status ;
    int i;
    
    status = true;
    list = msm4g_linkedlist_new();
    
    for (i=0;i<10;i++)
    {
        x[i].m = i;
        msm4g_linkedlist_add(list,&x[i]);
    }
    
    i=0;
    curr=list->head;
    while (curr != NULL)
    {
        y =  (Body *)curr->data;
        if (fabs(i - (y->m)) > DBL_EPSILON)
        {
            status = false;
            break;
        }
        curr = curr->next;
        i++;
    }
    msm4g_linkedlist_destroy(list);
    
    return status;
}


Boolean msm4g_unit_test_4()
{
    const int DIM = 3;
    const int N = 3; /**< Number of bodies in eight.ini */
    Body bodies[3];
    LinkedList *bodiesInFile;
    Body *bodyInFile;
    int i,j;
    
    for (i=0;i<N;i++)
    {
        msm4g_body_reset(&bodies[i]);
        bodies[i].m = 1.0;
    }
    bodies[0].r[0] =   0.97000436 ;
    bodies[0].r[1] =  -0.24308753 ;
    bodies[0].v[0] =   0.466203685 ;
    bodies[0].v[1] =   0.43236573 ;
    
    bodies[1].r[0] =  -0.97000436;
    bodies[1].r[1] =   0.24308753;
    bodies[1].v[0] =   0.466203685;
    bodies[1].v[1] =   0.43236573;
    
    bodies[2].v[0] =  -0.93240737;
    bodies[2].v[1] =  -0.86473146;
    
    
    bodiesInFile=msm4g_body_read("data/eight.ini");
   
    
    for (i=0;i<N;i++) /* for each body */
    {
        bodyInFile = msm4g_linkedlist_get(bodiesInFile,i);
        if (fabs(bodies[i].m - bodyInFile->m)>DBL_EPSILON )
            return false;
        for (j=0;j<DIM;j++) /* for each dimension */
        {
            if (fabs(bodies[i].r[j]-bodyInFile->r[j]) > DBL_EPSILON)
                return false;
            if (fabs(bodies[i].v[j]-bodyInFile->v[j]) > DBL_EPSILON)
                return false;
        }
    }
    
    /* Deallocate the bodies */
    for (i=0;i<N;i++)
    {
        bodyInFile = msm4g_linkedlist_get(bodiesInFile, i);
        free(bodyInFile);
    }
    msm4g_linkedlist_destroy(bodiesInFile);
    
    return true;
}

Boolean msm4g_unit_test_5()
{
    int i,n;
    Boolean status ;
    LinkedList *list ;
    Body **bodies;
    SimulationBox *box;
    D3Vector locationError;
    D3Vector locationExpected ;
    D3Vector widthExpected;
    D3Vector widthError;
    
    msm4g_d3vector_set(&locationExpected, -1.5, -3.0, -4.5);
    msm4g_d3vector_set(&widthExpected,     3.0,  6.0,  9.0);
    
    n=10;
    status = true;
    bodies = msm4g_body_rand(n);
    
    list = msm4g_linkedlist_new();
    for (i=0;i<n;i++)
    {
        msm4g_linkedlist_add(list, bodies[i]);
    }
    
    bodies[5]->r[0] = -1.0;
    bodies[5]->r[1] = -2.0;
    bodies[5]->r[2] = -3.0;
    
    bodies[7]->r[0] =  1.0;
    bodies[7]->r[1] =  2.0;
    bodies[7]->r[2] =  3.0;
    
    box = msm4g_box_new();
    msm4g_box_update(box, list, 0.5);

    msm4g_d3vector_daxpy(&locationError, &locationExpected,-1.0,&(box->location));
    msm4g_d3vector_daxpy(&widthError,    &widthExpected,   -1.0,&(box->width));
    
    if (msm4g_d3vector_norm(&locationError) > DBL_EPSILON ||
        msm4g_d3vector_norm(&widthError)    > DBL_EPSILON )
    {
        status = false;
    }
    
    msm4g_linkedlist_destroy(list);
    msm4g_body_destroyarray(bodies, n);
    msm4g_box_destroy(box);
    
    return status;
}

Boolean msm4g_unit_test_6()
{
    Boolean status = true;
    LinkedList *bodylist;
    LinkedList *binlist;
    SimulationBox box;
    double binwidth = 10;
    int i;
    
    bodylist = msm4g_body_read("data/bintest.ini");
    for (i=0;i<3;i++)
    {
        box.location.value[i] = 0.0;
        box.width.value[i]    = 30.0;
    }
    
    binlist=msm4g_bin_generate(&box,bodylist,binwidth);
    msm4g_force_short(binlist, 10.0);
    
    msm4g_bin_destroy(binlist);
    msm4g_linkedlist_destroyWithData(bodylist);
    return status;
}
