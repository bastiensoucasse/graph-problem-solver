/**
 * @file EdgeList.c
 * @author Vincent Penelle (vincent.penelle@u-bordeaux.fr)
 * @brief  Structure to store a list of graph edges that can be dynamically modified. Used as a temporary structure during parsing before translating into a more static structure.
 * @version 1
 * @date 2018-11-10
 * 
 * @copyright Creative Commons.
 * 
 */


#include "EdgeList.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/**
 * @brief Allocates space for edge list
 * @return The expression or NULL if not enough memory
 */
static SEdgeList *allocateEdgeList()
{
    SEdgeList *b = (SEdgeList *)malloc(sizeof(SEdgeList));

    if (b == NULL)
        return NULL;

    b->node1 = NULL;
    b->node2 = NULL;

    b->next = NULL;

    return b;
}


SEdgeList *addEdge(char *n1, char *n2, SEdgeList *list)
{
    SEdgeList *b = allocateEdgeList();

    if (b == NULL)
        return NULL;

    b->node1 = (char *)malloc((strlen(n1)+1)*sizeof(char));
    strcpy(b->node1,n1);
    b->node2 = (char *)malloc((strlen(n2)+1)*sizeof(char));
    strcpy(b->node2,n2);

    b->next = list;
	
    return b;
}

void printEdgeList(SEdgeList *e)
{
	if(e == NULL) { printf("\n"); return;}
	printf("(%s,%s) -- ",e->node1,e->node2);
	printEdgeList(e->next);
}

void deleteExpression(SEdgeList *b)
{
    if (b == NULL)
        return;

    deleteExpression(b->next);

    free(b->node1);
    free(b->node2);

    free(b);
}

