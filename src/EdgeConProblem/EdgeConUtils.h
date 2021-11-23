/**
 * @file EdgeConUtils.h
 * @author Quentin Desire
 * @author Iantsa Provost
 * @author Bastien Soucasse
 * @brief TODO…
 * @version 1
 * @date 2021-11-26
 *
 * @copyright Creative Commons
 *
 */
#ifndef EDGE_CON_UTILS_H
#define EDGE_CON_UTILS_H

#include "EdgeConGraph.h"

/**
 * @brief Gets the number of translators in the graph.
 *
 * @param graph An instance of the problem.
 * @return The number of translators in @p graph.
 * @pre @p graph must be a valid graph.
 */
int getNumTranslators(const EdgeConGraph graph);

#endif  // EDGE_CON_UTILS_H
