/**
 * @file EdgeConUtils.c
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
#include "EdgeConUtils.h"

#include <stdio.h>

#include "Graph.h"

int getNumTranslators(const EdgeConGraph graph) {
    int numTranslators = 0;
    printf("\n*** NUM TRANSLATORS ***\n");

    Graph g = getGraph(graph);
    int numNodes = orderG(g);
    printf("numNodes: %d\n", numNodes);
    for (int node1 = 0; node1 < numNodes; node1++)
        for (int node2 = node1 + 1; node2 < numNodes; node2++) {
            printf("Edge %d,%d", node1, node2);
            if (isTranslator(graph, node1, node2))
                printf(" is a translator.\n"), numTranslators++;
            else
                printf("…\n");
        }
    return numTranslators;
}
