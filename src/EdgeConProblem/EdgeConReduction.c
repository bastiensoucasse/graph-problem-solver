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


Z3_ast atMost(Z3_context ctx, Z3_ast *X, int size) 
{
    Z3_ast formula = Z3_mk_true(ctx);
    for (int i = 0; i < size; i++)
        for (int j = i + 1; j < size; j++)
        {
            Z3_ast clause = Z3_mk_or(ctx, Z3_mk_not(ctx, X[i]), Z3_mk_not(ctx, X[j]));
            formula = Z3_mk_and(formula, clause);
        }
    return formula;
}

Z3_ast atLeast(Z3_context ctx, Z3_ast *X, int size) 
{
    Z3_ast formula = Z3_mk_false(ctx);
    for (int i = 0; i < size; i++)
        formula = Z3_mk_or(formula, X[i]);
    return formula;
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
    return Z3_mk_false(ctx);
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
