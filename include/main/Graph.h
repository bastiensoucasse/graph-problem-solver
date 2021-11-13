/**
 * @file Graph.h
 * @author Vincent Penelle (vincent.penelle@u-bordeaux.fr)
 * @brief  Structure to store a graph statically, and to access to its informations easily.
  		   Includes source and target node informations, and colors.
 * @version 2
 * @date 2018-11-18, 2019-07-23, 2019-08-02, 2020-06-24
 * 
 * @copyright Creative Commons.
 * 
 */


#ifndef COCA_GRAPH_H_
#define COCA_GRAPH_H_

#include <stdbool.h>


/** @brief: the graph type. The first four fields are needed to represent a directed graph. The rest depends on needs. Here, the rest represents initial and final states of an automaton.*/
typedef struct {
	int numNodes; ///< The number of nodes of the graph.
	int numEdges; ///< The number of edges of the graph.
	char** nodes; ///< The names of nodes of the graph.
	bool* edges;	  ///< The edges of the graph.

//These are to handle colored graphs.
	int numColor;	///< The number of different colors.
	int* color;		///< The color of each node of the graph.
	char** colorNames;	///< The string associated with the colors (for display).

//This is only for dealing with automata. May be changed according to needs.
	bool *initial;	///< Array of source nodes.
	bool *final;	///< Array of target nodes.
} Graph;

/**
 * @brief Creates a copy of the graph passed in argument.
 * 
 * @param graph A graph.
 * @return graph A copy of graph.
 * @pre @p graph must be a valid graph.
 */
Graph copyGraph(Graph graph);

/**
 * @brief Displays a graph with a list of nodes and a matrix of edges.
 * 
 * @param graph the graph to display.
 * 
 * @pre @p graph must be a valid graph.
 */
void printGraph(Graph graph);

/**
 * @brief Frees all memory occupied by a graph.
 * 
 * @param graph The graph to delete.
 * 
 * @pre @p graph must be a valid graph.
 */
void deleteGraph(Graph graph);

/**
 * @brief Returns the number of nodes of @p graph.
 * 
 * @param graph A graph.
 * @return int Its number of nodes.
 * 
 * @pre @p graph must be a valid graph.
 */
int orderG(Graph graph);

/**
 * @brief Returns the number of edges of @p graph.
 * 
 * @param graph A graph.
 * @return int Its number of edges.
 * 
 * @pre @p graph must be a valid graph.
 */
int sizeG(Graph graph);

/**
 * @brief Tells if (@p source, @p target) is an edge in @p graph.
 * 
 * @param graph A graph.
 * @param source The source of the edge.
 * @param target The target of the edge.
 * @return true If the edge is present in @p graph.
 * @return false Otherwise.
 * @pre @p graph must be a valid graph.
 * @pre 0 <= @p source < @p graph.numNodes
 * @pre 0 <= @p target < @p graph.numNodes
 */
bool isEdge(Graph graph, int source, int target);

/**
 * @brief Get the Color of a node
 * 
 * @param graph A graph.
 * @param node The node.
 * @return int The color number of @p node.
 * @pre @p graph must be a valid graph.
 * @pre 0 <= @p node < @p graph.numNodes
 */
int getColor(Graph graph, int node);

/**
 * @brief Get the string representing the color @p col.
 * 
 * @param graph A graph.
 * @param col	A color number.
 * @return char* The string representing the color @p col.
 * @pre @p graph must be a valid graph.
 * @pre 0 <= @p col < @p graph.numColor
 */
char* getColorString(Graph graph, int col);

/**
 * @brief Get the number of colors of @p graph.
 * 
 * @param graph A graph
 * @return int The number of colors of its nodes.
 * @pre @p graph must be a valid graph.
 */
int getNumColor(Graph graph);

/**
 * @brief Tells if @p node is source in @p graph.
 * 
 * @param graph A graph.
 * @param node A node.
 * @return true If @p node is a source.
 * @return false Otherwise.
 * @pre @p graph must be a valid graph.
 * @pre 0 <= @p node < @p graph.numNodes
 */
bool isSource(Graph graph, int node);

/**
 * @brief Tells if @p node is target in «graph».
 * 
 * @param graph A graph.
 * @param node A node.
 * @return true If @p node is a target.
 * @return false Otherwise.
 * @pre @p graph must be a valid graph.
 * @pre 0 <= @p node < @p graph.numNodes
 */
bool isTarget(Graph graph, int node);

/**
 * @brief Returns the name of a node given its identifier.
 * 
 * @param graph A graph.
 * @param node A node identifier. Must be lower than orderG(@p graph).
 * @return char* The name of @p node.
 * @pre @p graph must be a valid graph.
 * @pre 0 <= @p node < @p graph.numNodes
 */
char* getNodeName(Graph graph, int node);


#endif /* DOT_PARSER_GRAPH_H_ */
