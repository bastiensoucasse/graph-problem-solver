/**
 * @file EdgeConResolution.c
 * @author Quentin Desire
 * @author Iantsa Provost
 * @author Bastien Soucasse
 * @brief  Algorithms to solve directly the EdgeCon problem
 * @version 1
 * @date 2021-11-19
 * 
 * @copyright Creative Commons.
 * 
 */
#include "EdgeConResolution.h"
#include "Graph.h"
#include <stdio.h>
#include <stdlib.h>

/**
 * @brief 
 * 
 * @param numNodes 
 * @param depth 
 * @return int 
 */
static int computeMaxDepth(int numNodes, int depth[])
{
    int max = -1;
    for (int i = 0; i < numNodes; i++)
    {
        if(depth[i] > max)
            max = depth[i];
    }
    return max;
}

/**
 * @brief 
 * 
 * @param graph 
 * @param u 
 * @param depth 
 * 
 * @pre graph must be valid.
 */
static void browse(EdgeConGraph graph, int u, int depth[])
{
    Graph graphInit = getGraph(graph);
    int numNodes = orderG(graphInit);
    for (int v = 0; v < numNodes; v++)
    {
        if (v != u && isEdge(graphInit, u, v) && depth[v] == - 1)
        {
            if (isEdgeHeterogeneous(graph, u, v))
                depth[v] = depth[u] + 1;
            else 
                depth[v] = depth[u];
            browse(graph, v, depth);
        }
    }
}

/**
 * @brief Brute Force Algorithm. If there is a result, the solution will be stored in @param graph, and its homogeneous components updated. If no solution, graph won't be modified. Returns the maximal cost of communication for any choice of translators.
 * 
 * @param graph An instance of the problem.
 * @return the maximal cost that two nodes communicate with for any possible set of transducers of minimal size. Returns -1 if there is no solution (i.e., the graph is not connex).
 * 
 * @pre graph must be valid.
 */
int BruteForceEdgeCon(EdgeConGraph graph)
{
    int max = 0;
    Graph graphInit = getGraph(graph);
    int numNodes = orderG(graphInit);
    int depth[numNodes];

    for (int i = 0; i < numNodes; i++)
    {
        for (int j = 0; j < numNodes; j++)
            depth[j] = -1;
        depth[i] = 0;
        browse(graph, i, depth);
        for (int k = 0; k < numNodes; k++)
        {
            if (depth[k] == -1)
                return 0;
            int maxDepth = computeMaxDepth(numNodes, depth);
            max = maxDepth > max ? maxDepth : max;
        }
    }
    printf("Result is : %d !\n", max);
    return max;
}
