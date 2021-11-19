/**
 * @file EdgeConReduction.c
 * @author Quentin Desire
 * @author Iantsa Provost
 * @author Bastien Soucasse
 * @brief An implementation of the Coca project for 2021 year. Converts a 2-coloured graph g and an integer k to a formula true if and only if it exists a translator set of minimal size which forces two nodes to communicate with cost strictly bigger than k.
 * Provides functions to generate the formula and the necessary variables, alongside function to decode a solution from a model of the formula.
 * @version 1
 * @date 2021-11-26
 * 
 * @copyright Creative Commons
 * 
 */
#include "EdgeConReduction.h"
#include "Z3Tools.h"
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

/**
 * @brief Creates a variable representing that the edge (@p node1, @p node2) has translator number @p number. The variable name will always start with the smaller node.
 * 
 * @param ctx The solver context
 * @param node1 A node
 * @param node2 A node
 * @param number A number
 * @return Z3_ast The lowest variable in lexicographic order among [(@p node1, @p node2),number] and [(@p node2, @p node1),number]
 */
Z3_ast getVariableIsIthTranslator(Z3_context ctx, int node1, int node2, int number)
{
    char name[40];
    if (node1 < node2)
        snprintf(name, 40, "[(%d,%d),%d]", node1, node2, number);
    else
        snprintf(name, 40, "[(%d,%d),%d]", node2, node1, number);
    return mk_bool_var(ctx, name);
}

/**
 * @brief Get the Variable p_{@p child,@p parent} representing that the connected component @p parent is the parent of the connected component @p child in the reduction.
 * 
 * @param ctx The solver context.
 * @param child The number of a connected component.
 * @param parent The number of a connected component.
 * @return Z3_ast A formula representing p_{@p child,@p parent}.
 */
Z3_ast getVariableParent(Z3_context ctx, int child, int parent)
{
    char name[40];
    snprintf(name, 40, "p_[%d,%d]", child, parent);
    return mk_bool_var(ctx, name);
}

/**
 * @brief Get the Variable l_{@p component,@p level} representing that the connected component @p component is at level @p level in the spanning tree.
 * 
 * @param ctx The solver context.
 * @param level The level.
 * @param component The number of the component.
 * @return Z3_ast A formula representing l_{@p component,@p level}
 */
Z3_ast getVariableLevelInSpanningTree(Z3_context ctx, int level, int component)
{
    char name[40];
    snprintf(name, 40, "l_[%d,%d]", component, level);
    return mk_bool_var(ctx, name);
}

/**
 * @brief Makes the set of x_{e,i} from the graph.
 * 
 * @param ctx The solver context.
 * @param g The graph as an EdgeConGraph.
 * @param n The number of vertices.
 * @param H_t The number of heterogenous edges.
 * @param N The number of translators.
 * @return The set of x_{e,i} created.
 */
Z3_ast **makeX(Z3_context ctx, EdgeConGraph g, int n, int H_t, int N)
{
    Z3_ast **X = (Z3_ast **)malloc(H_t * sizeof(Z3_ast *));
    for (int i = 0; i < H_t; i++)
        X[i] = (Z3_ast *)malloc(N * sizeof(Z3_ast));

    int e = 0;
    for (int u = 0; u < n; u++)
        for (int v = n - 1; v > u; v--)
            if (isEdgeHeterogeneous(g, u, v))
            {
                for (int i = 0; i < N; i++)
                {
                    X[e][i] = getVariableIsIthTranslator(ctx, u, v, i);
                    printf("Creating x_%d,%d (e = %d,%d)…\n", e, i, u, v);
                }
                e++;
            }
    return X;
}

/**
 * @brief Frees the set of x_{e,i}.
 * 
 * @param X The set of x_{e,i}.
 * @param H_t The number of heterogenous edges.
 */
void freeX(Z3_ast **X, int H_t)
{
    for (int i = 0; i < H_t; i++)
        free(X[i]), X[i] = NULL;
    free(X), X = NULL;
}

/**
 * @brief Creates a formula for at least one literal set to TRUE.
 * 
 * @param ctx The solver context.
 * @param X The set of x_{e,i}.
 * @param size The size of the set of x_{e,i}.
 * @return The final formula.
 */
Z3_ast atLeast(Z3_context ctx, Z3_ast *X, int size)
{
    Z3_ast formula = Z3_mk_false(ctx);
    for (int i = 0; i < size; i++)
    {
        Z3_ast args[2] = {formula, X[i]};
        formula = Z3_mk_or(ctx, 2, args);
    }
    return formula;
}

/**
 * @brief Creates a formula for at most one literal set to TRUE.
 * 
 * @param ctx The solver context.
 * @param X The set of x_{e,i}.
 * @param size The size of the set of x_{e,i}.
 * @return The final formula.
 */
Z3_ast atMost(Z3_context ctx, Z3_ast *X, int size)
{
    Z3_ast formula = Z3_mk_true(ctx);
    for (int i = 0; i < size; i++)
        for (int j = i + 1; j < size; j++)
        {
            Z3_ast args[2] = {Z3_mk_not(ctx, X[i]), Z3_mk_not(ctx, X[j])};
            Z3_ast clause = Z3_mk_or(ctx, 2, args);

            args[0] = formula;
            args[1] = clause;
            formula = Z3_mk_and(ctx, 2, args);
        }
    return formula;
}

/**
 * @brief Creates a formula that ensures that 
 * each translator can only be associated to exactly one edge,
 * and each edge can only receive at most one translator.
 
 * @param ctx The solver context.
 * @param X The set of x_{e,i}.
 * @param H_t The number of heterogeneous edges.
 * @param N The max number of translators.
 * @return Z3_ast The final formula.
 */
Z3_ast phi1(Z3_context ctx, Z3_ast **X, int H_t, int N)
{
    Z3_ast left = Z3_mk_true(ctx);
    Z3_ast right = Z3_mk_true(ctx);

    Z3_ast Xi[H_t];
    for (int i = 0; i < N; i++)
    {
        // Retrieve every x_e,i with i settled (the column i)
        for (int e = 0; e < H_t; e++)
            Xi[e] = X[e][i];

        Z3_ast argsLeft[2] = {left, atMost(ctx, Xi, H_t)};
        left = Z3_mk_and(ctx, 2, argsLeft);
    }

    for (int e = 0; e < H_t; e++)
    {
        Z3_ast argsRight[2] = {right, atMost(ctx, X[e], N)};
        right = Z3_mk_and(ctx, 2, argsRight);
    }

    Z3_ast args[2] = {left, right};
    return Z3_mk_and(ctx, 2, args);
}

/**
 * @brief Generates a SAT formula satisfiable if and only if there is a set of translators of cost @p cost such that the graph admits a valid path between any two nodes.
 * 
 * @param ctx The solver context.
 * @param graph A EdgeConGraph.
 * @param cost The cost of the translator set.
 * @return Z3_ast The formula.
 * @pre graph must be an initialized EdgeConGraph with computed connected components.
 */
Z3_ast EdgeConReduction(Z3_context ctx, EdgeConGraph edgeGraph, int cost)
{
    // Z3_ast X[5];
    // for (int i = 0; i < 5; i++)
    //     X[i] = Z3_mk_true(ctx);
    // return atMost(ctx, X, 5);

    int n = orderG(getGraph(edgeGraph));
    int H_t = getNumHeteregeneousEdges(edgeGraph);
    int N = getNumComponents(edgeGraph) - 1;
    printf("n = %d, H_t = %d, N = %d\n", n, H_t, N);

    Z3_ast **X = makeX(ctx, edgeGraph, n, H_t, N);

    // TODO…
    Z3_ast test = phi1(ctx, X, H_t, N);
    printf("%s\n", Z3_ast_to_string(ctx, test));

    freeX(X, H_t);
    return Z3_mk_true(ctx);
}

/**
 * @brief Gets the translator set from a model and adds it to the EdgeConGraph. It also computes the homogeneous components taking into account the translators (this number should be one).
 * 
 * @param ctx The solver context.
 * @param model A variable assignment.
 * @param graph A EdgeConGraph.
 * 
 * @pre @p model must be a valid model.
 * @pre @p graph must be a valid EdgeConGraph with no translators (or at least the original number of homogeneous components).
 */
void getTranslatorSetFromModel(Z3_context ctx, Z3_model model, EdgeConGraph graph)
{
    return;
}
