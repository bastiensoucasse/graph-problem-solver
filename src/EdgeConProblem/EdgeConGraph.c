#include "EdgeConGraph.h"
#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <string.h>

struct EdgeConGraph_s
{
    Graph graph;                  ///< The graph.
    bool *heterogeneousEdges;     ///< The edges that are heterogeneous (not taking into account the translators).
    bool **homogeneousComponents; ///< The homogeneous connected components (taking into account the translators). Each member is an array of size @p graph.numNodes that sets all nodes in a connected component to true and other to false. There are @p numComponents such arrays, and a node is true in exactly one array.
    bool *translators;            ///< The translator edges.
    int numComponents;            ///< The number of homogeneous components.
};

EdgeConGraph initializeGraph(Graph graph)
{
    EdgeConGraph result = (EdgeConGraph)malloc(sizeof(*result));
    result->graph = graph;
    result->heterogeneousEdges = (bool *)malloc(orderG(graph) * orderG(graph) * sizeof(bool));
    result->homogeneousComponents = (bool **)malloc(orderG(graph) * sizeof(bool *));
    for (int i = 0; i < orderG(graph); i++)
        result->homogeneousComponents[i] = (bool *)malloc(orderG(graph) * sizeof(bool));
    result->translators = (bool *)malloc(orderG(graph) * orderG(graph) * sizeof(bool));
    result->numComponents = 0;

    for (int i = 0; i < orderG(graph); i++)
    {
        for (int j = 0; j < orderG(graph); j++)
        {
            result->heterogeneousEdges[i * orderG(graph) + j] = isEdge(graph, i, j) && (getColor(graph, i) != getColor(graph, j));
        }
    }

    for (int i = 0; i < orderG(graph) * orderG(graph); i++)
        result->translators[i] = false;

    computesHomogeneousComponents(result);

    return result;
}

Graph getGraph(const EdgeConGraph graph) { return graph->graph; }

void resetTranslator(EdgeConGraph graph)
{
    graph->numComponents = 0;

    for (int i = 0; i < orderG(graph->graph) * orderG(graph->graph); i++)
        graph->translators[i] = false;

    computesHomogeneousComponents(graph);
}

void deleteEdgeConGraph(EdgeConGraph graph)
{
    free(graph->heterogeneousEdges);
    free(graph->translators);
    for (int i = 0; i < orderG(graph->graph); i++)
        free(graph->homogeneousComponents[i]);
    free(graph->homogeneousComponents);
    free(graph);
}

void setHeterogeneousEdge(EdgeConGraph graph, int node1, int node2)
{
    if (isEdge(graph->graph, node1, node2))
    {
        graph->heterogeneousEdges[node1 * orderG(graph->graph) + node2] = true;
        graph->heterogeneousEdges[node2 * orderG(graph->graph) + node1] = true;
    }
}

void SetHomogeneousEdge(EdgeConGraph graph, int node1, int node2)
{
    if (isEdge(graph->graph, node1, node2))
    {
        graph->heterogeneousEdges[node1 * orderG(graph->graph) + node2] = false;
        graph->heterogeneousEdges[node2 * orderG(graph->graph) + node1] = false;
    }
}

bool isEdgeHomogeneous(const EdgeConGraph graph, int node1, int node2)
{
    return (isEdge(graph->graph, node1, node2) && !graph->heterogeneousEdges[node1 * orderG(graph->graph) + node2]);
}

bool isEdgeHeterogeneous(const EdgeConGraph graph, int node1, int node2)
{
    return (isEdge(graph->graph, node1, node2) && graph->heterogeneousEdges[node1 * orderG(graph->graph) + node2]);
}

int getNumHeteregeneousEdges(const EdgeConGraph graph)
{
    int numNodes = orderG(graph->graph);
    int result = 0;
    for (int u = 0; u < numNodes; u++)
    {
        for (int v = u + 1; v < numNodes; v++)
        {
            if (isEdgeHeterogeneous(graph, u, v))
                result++;
        }
    }
    return result;
}

void addTranslator(EdgeConGraph graph, int node1, int node2)
{
    graph->translators[node1 * orderG(graph->graph) + node2] = true;
    graph->translators[node2 * orderG(graph->graph) + node1] = true;
}

void removeTranslator(EdgeConGraph graph, int node1, int node2)
{
    graph->translators[node1 * orderG(graph->graph) + node2] = false;
    graph->translators[node2 * orderG(graph->graph) + node1] = false;
}

bool isTranslator(const EdgeConGraph graph, int node1, int node2)
{
    return graph->translators[node1 * orderG(graph->graph) + node2];
}

void computesComponent(EdgeConGraph graph, int node, int component)
{
    if (graph->homogeneousComponents[component][node])
        return;
    graph->homogeneousComponents[component][node] = true;
    for (int i = 0; i < orderG(graph->graph); i++)
    {
        if (isEdgeHomogeneous(graph, node, i) || isTranslator(graph, node, i))
            computesComponent(graph, i, component);
    }
    return;
}

void computesHomogeneousComponents(EdgeConGraph graph)
{
    graph->numComponents = 0;
    int currentNode = 0;
    while (currentNode < orderG(graph->graph))
    {
        for (int i = 0; i < orderG(graph->graph); i++)
            graph->homogeneousComponents[graph->numComponents][i] = false;
        computesComponent(graph, currentNode, graph->numComponents);
        bool acc = true;
        while (currentNode < orderG(graph->graph) && acc)
        {
            acc = false;
            for (int i = 0; i <= graph->numComponents; i++)
            {
                acc = acc || graph->homogeneousComponents[i][currentNode];
            }
            if (acc)
                currentNode++;
        }
        (graph->numComponents)++;
    }
}

bool areInSameComponent(const EdgeConGraph graph, int node1, int node2)
{
    bool result = false;
    for (int i = 0; i < graph->numComponents; i++)
    {
        result = result || (graph->homogeneousComponents[i][node1] && graph->homogeneousComponents[i][node2]);
    }
    return result;
}

bool isNodeInComponent(const EdgeConGraph graph, int node, int component)
{
    return graph->homogeneousComponents[component][node];
}

int getNumComponents(const EdgeConGraph graph)
{
    return graph->numComponents;
}

void printTranslator(const EdgeConGraph graph)
{
    if (graph->numComponents != 1)
        printf("Warning, printing a set of translator that yields %d homogeneous components\n", graph->numComponents);
    printf("Translator edges: ");
    for (int i = 0; i < orderG(graph->graph); i++)
    {
        for (int j = i + 1; j < orderG(graph->graph); j++)
        {
            if (graph->translators[i * orderG(graph->graph) + j])
                printf("%s(%d)-%s(%d), ", getNodeName(graph->graph, i), i, getNodeName(graph->graph, j), j);
        }
    }
    printf("\n");
}

void createDotOfEdgeConGraph(const EdgeConGraph graph, char *name)
{
    if (graph->numComponents != 1)
        printf("Warning, printing a set of translator that yields %d homogeneous components\n", graph->numComponents);

    FILE *file;

    struct stat st = {0};
    if (stat("./sol", &st) == -1)
        mkdir("./sol", 0777);

    if (name == NULL)
    {
        char nameFile[30];
        snprintf(nameFile, 30, "sol/result.dot");
        file = fopen(nameFile, "w");
        fprintf(file, "graph Sol{\n");
    }
    else
    {
        int length = strlen(name) + 12;
        char nameFile[length];
        snprintf(nameFile, length, "sol/%s.dot", name);
        file = fopen(nameFile, "w");
        fprintf(file, "graph %s{\n", name);
    }

    for (int node = 0; node < orderG(graph->graph); node++)
    {
        fprintf(file, "%s", getNodeName(graph->graph, node));
        if (strcmp(getColorString(graph->graph, getColor(graph->graph, node)), "") != 0)
        {
            fprintf(file, "[color=%s]", getColorString(graph->graph, getColor(graph->graph, node)));
        }
    }

    for (int node = 0; node < orderG(graph->graph); node++)
    {
        for (int node2 = 0; node2 < node; node2++)
        {
            if (isEdge(graph->graph, node, node2))
            {
                fprintf(file, "%s -- %s", getNodeName(graph->graph, node), getNodeName(graph->graph, node2));
                if (isTranslator(graph, node, node2))
                    fprintf(file, "[color=blue]");
                fprintf(file, ";\n");
            }
        }
    }

    fprintf(file, "}\n");

    fclose(file);
}