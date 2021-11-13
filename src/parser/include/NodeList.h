/**
 * @file NodeList.h
 * @author Vincent Penelle (vincent.penelle@u-bordeaux.fr)
 * @brief  Structure to store a list of graph nodes that can be dynamically modified. Used as a temporary structure during parsing before translating into a more static structure.
 *         Includes automata features (initial and final nodes), and color of nodes.
 * @version 2
 * @date 2019-06-24
 * 
 * @copyright Creative Commons.
 * 
 */

#ifndef COCA_NODELIST_H_
#define COCA_NODELIST_H_

#include <stdbool.h>

/**
 * @brief The NodeList structure
 */
typedef struct tagSNodeList
{
	char* node;
    bool initial;
    bool final;
    char* color;
    struct tagSNodeList *next;
} SNodeList;

/**
 * @brief Adds a node in front of a list (works if list is null).
 * @param n1 the node
 * @param isInit tells if the node is initial
 * @param isFinal tells if the node is final
 * @param list the list to append to
 * @param col the color of the node
 * @return the new list or NULl in case of no memory.
 */

SNodeList *addNode(char* n1,bool isInit, bool isFinal, char* col, SNodeList *list);

/**
 * @brief If n is present in the list, updates its initial and final status with the arguments given. Otherwise, adds the node at the end of the list.
 * 
 * @param n the node to modify or add.
 * @param isInit tells if the node is initial.
 * @param isFinal tells if the node is final.
 * @param col the color of the node.
 * @param list the list to modify.
 */
void addOrUpdateNode(char *n, bool isInit, bool isFinal, char* col, SNodeList *list);

/**
 * @brief Prints a NodeList.
 * 
 * @param e a NodeList to print.
 */
void printNodeList(SNodeList *e);

/**
 * @brief Deletes a node list
 * @param b The node list
 */
void deleteNodeList(SNodeList *b);


#endif /* DOT_PARSER_NODELIST_H_ */
