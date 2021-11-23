/**
 * @file EdgeConResolution.c
 * @author Quentin Desire
 * @author Iantsa Provost
 * @author Bastien Soucasse
 * @brief Algorithms to solve directly the EdgeCon problem
 * @version 1
 * @date 2021-11-26
 *
 * @copyright Creative Commons.
 *
 */
#include "EdgeConResolution.h"

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

#include "EdgeConUtils.h"
#include "Graph.h"

/**
 * @brief Compare each depth, stock and then return the highest value.
 *
 * @param numNodes Number of nodes in the EdgeConGraph.
 * @param depth An array filled with depth of each nodes.
 * @return the maximum of depth array.
 */
static int computeMaxDepth(int numNodes, int *depth) {
    int max = -1;

    for (int i = 0; i < numNodes; i++)
        if (depth[i] > max) max = depth[i];

    return max;
}

/**
 * @brief Browse neighbours of node u given and register the depth of each nodes
 * (one depth per homogeneous components). Recursively loop on all the graph.
 *
 * @param graph An instance of the problem.
 * @param u Node u that we will explore each neighbours (if it hasn't been
 * explored yet).
 * @param depth An array filled with depth of each nodes.
 */
static void browse(EdgeConGraph graph, int u, int *depth) {
    Graph g = getGraph(graph);
    int numNodes = orderG(g);

    for (int v = 0; v < numNodes; v++) {
        if (v != u && isEdge(g, u, v) && depth[v] == -1) {
            if (isTranslator(graph, u, v)) {  // Should be heterogeneous
                depth[v] = depth[u] + 1;
                browse(graph, v, depth);
            } else if (isEdgeHomogeneous(graph, u, v)) {
                depth[v] = depth[u];
                browse(graph, v, depth);
            }
        }
    }
}

/**
 * @brief Get the highest cost of communication between two nodes with a set of
 * choosen transducers.
 *
 * @param graph An instance of the problem.
 * @return the maximal cost that cost that two nodes communicate with for the
 * set of transducers translators.
 */
static int getMaxPathLength(EdgeConGraph graph) {
    int max = 0;

    Graph g = getGraph(graph);
    int numNodes = orderG(g);
    int depth[numNodes];

    for (int i = 0; i < numNodes; i++) {
        for (int j = 0; j < numNodes; j++) depth[j] = -1;
        depth[i] = 0;
        browse(graph, i, depth);
        for (int k = 0; k < numNodes; k++) {
            if (depth[k] == -1) return 0;
            int maxDepth = computeMaxDepth(numNodes, depth);
            max = maxDepth > max ? maxDepth : max;
        }
    }

    return max;
}

/**
 * @brief Apply the next set of translators of size N in the graph.
 *
 * @param graph An instance of the problem.
 * @param N The number of translators we need to have.
 * @param n The number of vertices in the graph.
 * @param H_t The number of heterogeneous edges in the graph.
 * @param hE The custom heteregeneous edges matrix.
 * @return true if a next translator set has been found, false otherwise.
 */
static bool nextTranslatorSet(EdgeConGraph graph, int N, int n, int H_t,
                              int **hE) {
    do {
        for (int i = 0; i < H_t; i++) {
            int u = hE[0][i];
            int v = hE[1][i];
            if (!isTranslator(graph, u, v)) {
                addTranslator(graph, u, v);
                break;
            } else {
                removeTranslator(graph, u, v);
                if (i == H_t - 1) return false;
            }
        }
    } while (getNumTranslators(graph) != N);

    return true;
}

/**
 * @brief Creates the custom heterogeneous edges matrix.
 *
 * @param graph An instance of the problem.
 * @param n The number of vertices.
 * @param H_t The number of heterogeneous edges.
 * @return a pointer on the created matrix.
 */
static int **createHE(EdgeConGraph graph, int n, int H_t) {
    int **hE = (int **)malloc(2 * sizeof(int *));
    hE[0] = (int *)malloc(H_t * sizeof(int));
    hE[1] = (int *)malloc(H_t * sizeof(int));

    int i = 0;
    for (int u = 0; u < n; u++)
        for (int v = u + 1; v < n; v++)
            if (isEdgeHeterogeneous(graph, u, v)) {
                hE[0][i] = u;
                hE[1][i] = v;
                i++;
            }

    return hE;
}

/**
 * @brief Frees the heterogeneous edges matrix.
 *
 * @param hE The heterogeneous edges matrix.
 */
static void freeHE(int **hE) {
    for (int i = 0; i < 2; i++) free(hE[i]), hE[i] = NULL;
    free(hE), hE = NULL;
}

/**
 * @brief Brute Force Algorithm. If there is a result, the solution will be
 * stored in @param graph, and its homogeneous components updated. If no
 * solution, graph won't be modified. Returns the maximal cost of communication
 * for any choice of translators.
 *
 * @param graph An instance of the problem.
 * @return the maximal cost that two nodes communicate with for any possible set
 * of transducers of minimal size. Returns -1 if there is no solution (i.e., the
 * graph is not connex).
 *
 * @pre graph must be valid.
 */
int BruteForceEdgeCon(EdgeConGraph graph) {
    int min = INT_MAX;

    int n = orderG(getGraph(graph));            // Number of vertices
    int H_t = getNumHeteregeneousEdges(graph);  // Number of heterogeneous edges
    int N = getNumComponents(graph) - 1;  // Number of homogeneous components
    int **hE = createHE(graph, n, H_t);   // Heterogeneous edges custom

    for (int e = 0; e < N; e++) addTranslator(graph, hE[0][e], hE[1][e]);

    do {
        int k = getMaxPathLength(graph);
        if (k < min && k != 0) min = k;
    } while (nextTranslatorSet(graph, N, n, H_t, hE));

    resetTranslator(graph);
    freeHE(hE);
    return min;
}
