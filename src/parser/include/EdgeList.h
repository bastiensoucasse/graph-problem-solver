/**
 * @file EdgeList.h
 * @author Vincent Penelle (vincent.penelle@u-bordeaux.fr)
 * @brief  Structure to store a list of graph edges that can be dynamically modified. Used as a temporary structure during parsing before translating into a more static structure.
 * @version 1
 * @date 2018-11-10
 * 
 * @copyright Creative Commons.
 * 
 */

#ifndef COCA_EDGELIST_H_
#define COCA_EDGELIST_H_


/**
 * @brief The EdgeList structure
 */
typedef struct tagSEdgeList
{
	char* node1;
	char* node2;
    struct tagSEdgeList *next;
} SEdgeList;

/**
 * @brief Adds an edge in front of a list (works if list is null).
 * @param n1 the left node
 * @param n2 the right node
 * @param list the list to append to
 * @return the new list or NULl in case of no memory.
 */

SEdgeList *addEdge(char* n1, char* n2, SEdgeList *list);

/**
 * @brief Prints an EdgeList.
 * 
 * @param e an EdgeList to print.
 */
void printEdgeList(SEdgeList *e);

/**
 * @brief Deletes an edge list
 * @param b The edge list
 */
void deleteExpression(SEdgeList *b);


#endif /* DOT_PARSER_EDGELIST_H_ */
