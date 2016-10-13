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
    FILE *fp;
    const int DIM = 3;
    const int N = 3; /**< Number of bodies in eight.ini */
    double mass,r[3],v[3];
    Body bodies[3], bodiesInFile[3];
    int i,j;
    int ibody;
    
    for (i=0;i<N;i++)
    {
        msm4g_body_reset(&bodies[i]);
        msm4g_body_reset(&bodiesInFile[i]);
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
    
    fp = fopen("data/eight.ini","r");
    if (fp == NULL) return false;
    
    ibody = 0;
    while (true)
    {
        int ret = fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf",&mass,&r[0],&r[1],&r[2],&v[0],&v[1],&v[2]);
        if (ret == 7)
        {
            bodiesInFile[ibody].m = mass;
            for (i=0;i<DIM;i++)
            {
                bodiesInFile[ibody].r[i] = r[i];
                bodiesInFile[ibody].v[i] = v[i];
            }
            ibody++;
        } else if (ret == EOF)
            break;
        else
        {
            fclose(fp);
            return false;
        }
    }
    
    for (i=0;i<N;i++) /* for each body */
    {
        if (fabs(bodies[i].m-bodiesInFile[i].m)>DBL_EPSILON )
            return false;
        for (j=0;j<DIM;j++) /* for each dimension */
        {
            if (fabs(bodies[i].r[j]-bodiesInFile[i].r[j]) > DBL_EPSILON)
                return false;
            if (fabs(bodies[i].v[j]-bodiesInFile[i].v[j]) > DBL_EPSILON)
                return false;
        }
    }
    fclose(fp);
    return true;
}

Boolean msm4g_unit_test_5()
{
    int i,n;
    Boolean status ;
    LinkedList *list ;
    Body *bodies;
    SimulationBox *box;
    D3Vector z;
    D3Vector expectedLocation ;
    D3Vector expectedWidth    ;
    
    msm4g_d3vector_set(&expectedLocation, -1.5, -3.0, -4.5);
    msm4g_d3vector_set(&expectedWidth,     3.0,  6.0,  9.0);
    
    n=10;
    status = true;
    bodies = msm4g_body_rand(n);
    
    list = msm4g_linkedlist_new();
    for (i=0;i<n;i++)
    {
        msm4g_linkedlist_add(list, bodies+i);
    }
    
    (bodies+5)->r[0] = -1.0;
    (bodies+5)->r[1] = -2.0;
    (bodies+5)->r[2] = -3.0;
    
    (bodies+7)->r[0] =  1.0;
    (bodies+7)->r[1] =  2.0;
    (bodies+7)->r[2] =  3.0;
    
    box = msm4g_box_new();
    msm4g_box_update(box, list, 0.0);
    msm4g_box_print(box);
    msm4g_box_update(box, list, 0.5);
    msm4g_box_print(box);

    msm4g_d3vector_daxpy(&z, -1.0, &(box->location), &expectedLocation);
    
    
    free(bodies);
    msm4g_box_destroy(box);
    msm4g_linkedlist_destroy(list);
    
    return status;
}
