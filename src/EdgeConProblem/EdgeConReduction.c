/**
 * @file EdgeConReduction.c
 * @author Quentin Desire
 * @author Iantsa Provost
 * @author Bastien Soucasse
 * 
 * @brief TODO
 */

#include "EdgeConReduction.h"
#include "Z3Tools.h"
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

Z3_ast getVariableIsIthTranslator(Z3_context ctx, int node1, int node2, int number)
{
    char name[40];
    if (node1 < node2)
        snprintf(name, 40, "[(%d,%d),%d]", node1, node2, number);
    else
        snprintf(name, 40, "[(%d,%d),%d]", node2, node1, number);
    return mk_bool_var(ctx, name);
}

Z3_ast getVariableParent(Z3_context ctx, int child, int parent)
{
    char name[40];
    snprintf(name, 40, "p_[%d,%d]", child, parent);
    return mk_bool_var(ctx, name);
}

Z3_ast getVariableLevelInSpanningTree(Z3_context ctx, int level, int component)
{
    char name[40];
    snprintf(name, 40, "l_[%d,%d]", component, level);
    return mk_bool_var(ctx, name);
}

Z3_ast EdgeConReduction(Z3_context ctx, EdgeConGraph edgeGraph, int cost)
{
    return Z3_mk_false(ctx);
}

void getTranslatorSetFromModel(Z3_context ctx, Z3_model model, EdgeConGraph graph)
{
    return;
}
