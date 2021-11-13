/**
 * @file GraphListToGraph.h
 * @author Vincent Penelle (vincent.penelle@u-bordeaux.fr)
 * @brief  Functions to translate a GraphList into a graph.
 * @version 2
 * @date 2020-06-24
 * 
 * @copyright Creative Commons.
 * 
 */

#ifndef COCA_GRAPHLISTTOGRAPH_H_
#define COCA_GRAPHLISTTOGRAPH_H_

#include "GraphList.h"
#include "Graph.h"


/**
 * @brief Creates a Graph object from a GraphList. Does NOT free the source, so it must be destroyed independently.
 * 
 * @param source the GraphList to reinterpret as a graph.
 * @return Graph the graph corresponding to the source.
 */
Graph createGraph(GraphList source);

/**
 * @brief Auxilary function to determine the index of a node name.
 * 
 * @param list the array of node to search in.
 * @param size the size of list.
 * @param target the node to search
 * @return int the index of target in list.
 */
int findNode(char **list,int size,char *target);

#endif /* DOT_PARSER_GRAPHLISTTOGRAPH_H_ */