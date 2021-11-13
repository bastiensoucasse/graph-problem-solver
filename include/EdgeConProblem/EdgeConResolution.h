/**
 * @file EdgeConResolution.h
 * @author Vincent Penelle (vincent.penelle@u-bordeaux.fr)
 * @brief  Algorithms to solve directly the EdgeCon problem
 * @version 1
 * @date 2021-10-18
 * 
 * @copyright Creative Commons.
 * 
 */

#ifndef COCA_EDGECONRESOLUTION_H
#define COCA_EDGECONRESOLUTION_H

#include "EdgeConGraph.h"

/**
 * @brief Brute Force Algorithm. If there is a result, the solution will be stored in @param graph, and its homogeneous components updated. If no solution, graph won't be modified. Returns the maximal cost of communication for any choice of translators.
 * 
 * @param graph An instance of the problem.
 * @return the maximal cost that two nodes communicate with for any possible set of transducers of minimal size. Returns -1 if there is no solution (i.e., the graph is not connex).
 * 
 * @pre graph must be valid.
 */
int BruteForceEdgeCon(EdgeConGraph graph);

#endif