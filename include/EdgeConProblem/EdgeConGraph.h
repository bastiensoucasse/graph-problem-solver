/**
 * @file EdgeConGraph.h
 * @author Vincent Penelle (vincent.penelle@u-bordeaux.fr)
 * @brief  Structure to store a graph with heterogeneous edges informations for the EdgeCon Problem
 * @version 1
 * @date 2021-10-18
 * 
 * @copyright Creative Commons.
 * 
 */

#ifndef COCA_EDGECONGRAPH_H
#define COCA_EDGECONGRAPH_H

#include <stdbool.h>
#include "Graph.h"

/**
 * @brief The struct containing a graph, its heteregeneous edges, homogenous components, and current translator set.
 */

typedef struct EdgeConGraph_s *EdgeConGraph;

/**
 * @brief Initializes the struct from a graph. Does NOT creates a copy of the graph. Starts with empty translators.
 * 
 * @param graph A graph
 * @return EdgeConGraph The corresponding EdgeConGraph with no translator.
 * @pre @p graph must be a valid graph.
 */
EdgeConGraph initializeGraph(Graph graph);

/**
 * @brief Get the Graph to manipulate it.
 * 
 * @param graph An EdgeConGraph
 * @return Graph The graph without translators and components informations.
 */
Graph getGraph(const EdgeConGraph graph);

/**
 * @brief Resets the translator set to empty, and recomputes heterogeneous edges and connected components.
 * 
 * @param graph A EdgeConGraph.
 * @pre @p graph must be valid.
 */
void resetTranslator(EdgeConGraph graph);

/**
 * @brief Frees the additional informations of a EdgeConGraph, but NOT the graph (cf initialize).
 * 
 * @param graph A EdgeConGraph.
 * 
 * @pre @p graph must be valid.
 */
void deleteEdgeConGraph(EdgeConGraph graph);

/**
 * @brief Set the edge (@p node1,@p node2) as heterogeneous (if it exists, otherwise does nothing).
 * 
 * @param graph A EdgeConGraph.
 * @param node1 A node.
 * @param node2 A node.
 * 
 * @pre @p graph must be valid.
 * @pre @p node1 < @p graph->graph.numNodes.
 * @pre @p node2 < @p graph->graph.numNodes.
 */
void setHeterogeneousEdge(EdgeConGraph graph, int node1, int node2);

/**
 * @brief Set the edge (@p node1,@p node2) as heterogeneous (if it exists, otherwise does nothing).
 * 
 * @param graph A EdgeConGraph.
 * @param node1 A node.
 * @param node2 A node.
 * 
 * @pre @p graph must be valid.
 * @pre @p node1 < @p graph->graph.numNodes.
 * @pre @p node2 < @p graph->graph.numNodes.
 */
void SetHomogeneousEdge(EdgeConGraph graph, int node1, int node2);

/**
 * @brief
 * 
 * @param graph A EdgeConGraph.
 * @param node1 A node.
 * @param node2 A node.
 * @return True if (@p node1, @p node2) is an homogeneous edge.
 * 
 * @pre @p graph must be valid.
 * @pre @p node1 < @p graph.graph.numNodes.
 * @pre @p node2 < @p graph.graph.numNodes.
 */
bool isEdgeHomogeneous(const EdgeConGraph graph, int node1, int node2);

/**
 * @brief
 * 
 * @param graph A EdgeConGraph.
 * @param node1 A node.
 * @param node2 A node.
 * @return True if (@p node1, @p node2) is an heterogeneous edge.
 * 
 * @pre @p graph must be valid.
 * @pre @p node1 < @p graph.graph.numNodes.
 * @pre @p node2 < @p graph.graph.numNodes.
 */
bool isEdgeHeterogeneous(const EdgeConGraph graph, int node1, int node2);

/**
 * @brief Get the number of heterogeneous edges.
 * 
 * @param graph An EdgeConGraph
 * @return int the number of heterogeneous edges.
 */
int getNumHeteregeneousEdges(const EdgeConGraph graph);

/**
 * @brief Sets the edge (@p node1,@p node2) to be a translator. Warning: DOES NOT recompute homogeneous components, you must use computesHomogeneousComponents if you wont them to be updated.
 * 
 * @param graph A EdgeConGraph.
 * @param node1 A node.
 * @param node2 A node.
 * 
 * @pre graph must be valid.
 * @pre node < graph.graph.numNodes.
 */
void addTranslator(EdgeConGraph graph, int node1, int node2);

/**
 * @brief Sets the edge (@p node1,@p node2) to not be a translator. Warning: DOES NOT recompute homogeneous components, you must use computesHomogeneousComponents if you wont them to be updated.
 * 
 * @param graph A EdgeConGraph.
 * @param node1 A node.
 * @param node2 A node.
 * 
 * @pre graph must be valid.
 * @pre node < graph.graph.numNodes.
 */
void removeTranslator(EdgeConGraph graph, int node1, int node2);

/**
 * @brief Tells if (@p node1,@p node2) bears a translator.
 * 
 * @param graph An EdgeConGraph
 * @param node1 A node
 * @param node2 A node
 * @return true if (@p node1,@p node2) bears a translator
 * @return false otherwise
 */
bool isTranslator(const EdgeConGraph graph, int node1, int node2);

/**
 * @brief Computes the homogeneous components.
 * 
 * @param graph a EdgeConGraph.
 * 
 * @pre graph must be a valid EdgeConGraph.
 */
void computesHomogeneousComponents(EdgeConGraph graph);

/**
 * @brief Tests if @p node1 and @p node2 are in the same homogeneous component.
 * 
 * @param graph A EdgeConGraph.
 * @param node1 A node.
 * @param node2 A node.
 * @return true If @p node1 and @p node2 are in the same homogeneous component.
 * @return false Otherwise.
 * 
 * @pre @p graph must be a valid EdgeConGraph
 * @pre @p node1 < @p graph.graph.numNodes
 * @pre @p node2 < @p graph.graph.numNodes
 */
bool areInSameComponent(const EdgeConGraph graph, int node1, int node2);

/**
 * @brief Tells if @p node is in the homogeneous component @p component
 * 
 * @param graph An EdgeConGraph
 * @param node A node
 * @param component The number of a component
 * @return true if @p node is in component @p component
 * @return false otherwise
 */
bool isNodeInComponent(const EdgeConGraph graph, int node, int component);

/**
 * @brief Get the number of homogeneous components in @p graph.
 * 
 * @param graph An EdgeConGraph
 * @return int The number of homogeneous components
 */
int getNumComponents(const EdgeConGraph graph);

/**
 * @brief Displays the translator set of the EdgeConGraph.
 * 
 * @param graph A EdgeConGraph.
 * @pre @p graph must be a valid EdgeConGraph.
 */
void printTranslator(const EdgeConGraph graph);

/**
 * @brief Creates the file ("%s.dot",name) representing the solution to the problem described by @p set, or ("result.dot") if name is NULL.
 * 
 * @param graph A graph.
 * @param name The name of the output file.
 * @pre @p graph must be a valid graph.
 */
void createDotOfEdgeConGraph(const EdgeConGraph graph, char *name);

#endif