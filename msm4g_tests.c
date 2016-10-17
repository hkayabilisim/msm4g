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
    msm4g_linkedlist_add(list, (void *)&msm4g_unit_test_7);
    
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
    Particle x[10];
    int i;
    
    list = msm4g_linkedlist_new();
    
    for (i=0;i<10;i++)
    {
        msm4g_particle_reset(x+i);
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
    Particle x[10];
    Particle *y;
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
        y =  (Particle *)curr->data;
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
    Particle x[10];
    Particle *y;
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
        y =  (Particle *)curr->data;
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
    const int N = 3; /**< Number of particles in eight.ini */
    Particle particles[3];
    LinkedList *particlesInFile;
    Particle *particleInFile;
    int i,j;
    
    for (i=0;i<N;i++)
    {
        msm4g_particle_reset(&particles[i]);
        particles[i].m = 1.0;
    }
    particles[0].r[0] =   0.97000436 ;
    particles[0].r[1] =  -0.24308753 ;
    particles[0].v[0] =   0.466203685 ;
    particles[0].v[1] =   0.43236573 ;
    
    particles[1].r[0] =  -0.97000436;
    particles[1].r[1] =   0.24308753;
    particles[1].v[0] =   0.466203685;
    particles[1].v[1] =   0.43236573;
    
    particles[2].v[0] =  -0.93240737;
    particles[2].v[1] =  -0.86473146;
    
    
    particlesInFile=msm4g_particle_read("data/eight.ini");
   
    
    for (i=0;i<N;i++) /* for each particle */
    {
        particleInFile = msm4g_linkedlist_get(particlesInFile,i);
        if (fabs(particles[i].m - particleInFile->m)>DBL_EPSILON )
            return false;
        for (j=0;j<DIM;j++) /* for each dimension */
        {
            if (fabs(particles[i].r[j]-particleInFile->r[j]) > DBL_EPSILON)
                return false;
            if (fabs(particles[i].v[j]-particleInFile->v[j]) > DBL_EPSILON)
                return false;
        }
    }
    
    /* Deallocate the particles */
    for (i=0;i<N;i++)
    {
        particleInFile = msm4g_linkedlist_get(particlesInFile, i);
        free(particleInFile);
    }
    msm4g_linkedlist_destroy(particlesInFile);
    
    return true;
}

Boolean msm4g_unit_test_5()
{
    int i,n;
    Boolean status ;
    LinkedList *list ;
    Particle **particles;
    SimulationBox *box;
    D3Vector locationError;
    D3Vector locationExpected ;
    D3Vector widthExpected;
    D3Vector widthError;
    
    msm4g_d3vector_set(&locationExpected, -1.5, -3.0, -4.5);
    msm4g_d3vector_set(&widthExpected,     3.0,  6.0,  9.0);
    
    n=10;
    status = true;
    particles = msm4g_particle_rand(n);
    
    list = msm4g_linkedlist_new();
    for (i=0;i<n;i++)
    {
        msm4g_linkedlist_add(list, particles[i]);
    }
    
    particles[5]->r[0] = -1.0;
    particles[5]->r[1] = -2.0;
    particles[5]->r[2] = -3.0;
    
    particles[7]->r[0] =  1.0;
    particles[7]->r[1] =  2.0;
    particles[7]->r[2] =  3.0;
    
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
    msm4g_particle_destroyarray(particles, n);
    msm4g_box_destroy(box);
    
    return status;
}

Boolean msm4g_unit_test_6()
{
    Boolean status = true;
    LinkedList *particlelist;
    LinkedList *binlist;
    SimulationBox box;
    double binwidth = 10;
    int i;
    
    particlelist = msm4g_particle_read("data/bintest.ini");
    for (i=0;i<3;i++)
    {
        box.location.value[i] = 0.0;
        box.width.value[i]    = 30.0;
    }
    
    binlist=msm4g_bin_generate(&box,particlelist,binwidth);
    msm4g_force_short(binlist, 10.0);
    
    msm4g_bin_destroy(binlist);
    msm4g_linkedlist_destroyWithData(particlelist);
    return status;
}

Boolean msm4g_unit_test_7()
{
    Boolean status = true;
    int i,j;
    int inp[6] = {0, 1, 2, 3,  4, 5  };
    int out[6] = {1, 1, 2, 6, 24, 120};
    int x[4][4];
    int vector2d[2];
    
    
    
    /* Checking factorial function */
    for (i=0; i<6; i++)
        if (msm4g_math_factorial(inp[i]) != out[i])
            return false;
    for (i=0; i<4; i++)
    {
        for (j=0; j<4; j++)
        {
            vector2d[0] = i; vector2d[1]=j;
            x[i][j] = msm4g_math_cantor(vector2d,2);
            printf("%2d ",x[i][j]);
            
        }
        printf("\n");
        
    }
    return status;
}
