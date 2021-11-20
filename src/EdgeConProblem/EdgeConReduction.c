/**
 * @file EdgeConReduction.c
 * @author Quentin Desire
 * @author Iantsa Provost
 * @author Bastien Soucasse
 * @brief An implementation of the Coca project for 2021 year. Converts a
 * 2-coloured graph g and an integer k to a formula true if and only if it
 * exists a translator set of minimal size which forces two nodes to communicate
 * with cost strictly bigger than k. Provides functions to generate the formula
 * and the necessary variables, alongside function to decode a solution from a
 * model of the formula.
 * @version 1
 * @date 2021-11-26
 *
 * @copyright Creative Commons
 *
 */
#include "EdgeConReduction.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "Z3Tools.h"

/**
 * @brief Creates a variable representing that the edge (@p node1, @p node2) has
 * translator number @p number. The variable name will always start with the
 * smaller node.
 *
 * @param ctx The solver context
 * @param node1 A node
 * @param node2 A node
 * @param number A number
 * @return The lowest variable in lexicographic order among [(@p node1,
 * @p node2),number] and [(@p node2, @p node1),number]
 */
Z3_ast getVariableIsIthTranslator(Z3_context ctx, int node1, int node2,
                                  int number) {
    char name[40];
    if (node1 < node2)
        snprintf(name, 40, "[(%d,%d),%d]", node1, node2, number);
    else
        snprintf(name, 40, "[(%d,%d),%d]", node2, node1, number);
    return mk_bool_var(ctx, name);
}

/**
 * @brief Get the Variable p_{@p child,@p parent} representing that the
 * connected component @p parent is the parent of the connected component @p
 * child in the reduction.
 *
 * @param ctx The solver context.
 * @param child The number of a connected component.
 * @param parent The number of a connected component.
 * @return A formula representing p_{@p child,@p parent}.
 */
Z3_ast getVariableParent(Z3_context ctx, int child, int parent) {
    char name[40];
    snprintf(name, 40, "p_[%d,%d]", child, parent);
    return mk_bool_var(ctx, name);
}

/**
 * @brief Get the Variable l_{@p component,@p level} representing that the
 * connected component @p component is at level @p level in the spanning tree.
 *
 * @param ctx The solver context.
 * @param level The level.
 * @param component The number of the component.
 * @return A formula representing l_{@p component,@p level}
 */
Z3_ast getVariableLevelInSpanningTree(Z3_context ctx, int level,
                                      int component) {
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
Z3_ast **makeX(Z3_context ctx, EdgeConGraph g, int n, int H_t, int N) {
    Z3_ast **X = (Z3_ast **)malloc(H_t * sizeof(Z3_ast *));
    for (int i = 0; i < H_t; i++) X[i] = (Z3_ast *)malloc(N * sizeof(Z3_ast));

    int e = 0;
    for (int u = 0; u < n; u++)
        for (int v = n - 1; v > u; v--)
            if (isEdgeHeterogeneous(g, u, v)) {
                for (int i = 0; i < N; i++) {
                    X[e][i] = getVariableIsIthTranslator(ctx, u, v, i);
                    // printf("Creating x_(%d,%d),%d…\n", u, v, i);
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
void freeX(Z3_ast **X, int H_t) {
    for (int i = 0; i < H_t; i++) free(X[i]), X[i] = NULL;
    free(X), X = NULL;
}

/**
 * @brief Makes the set of p_{i,j} from the graph.
 *
 * @param ctx The solver context.
 * @param g The graph as an EdgeConGraph.
 * @param C_H The number of homogeneous components.
 * @return The set of p_{i,j} created.
 */
Z3_ast **makeP(Z3_context ctx, EdgeConGraph g, int C_H) {
    Z3_ast **P = (Z3_ast **)malloc(C_H * sizeof(Z3_ast *));
    for (int i = 0; i < C_H; i++) P[i] = (Z3_ast *)malloc(C_H * sizeof(Z3_ast));

    for (int j1 = 0; j1 < C_H; j1++)
        for (int j2 = 0; j2 < C_H; j2++)
        {
            P[j1][j2] = getVariableParent(ctx, j1, j2);
            // printf("Creating p_%d,%d…\n", j1, j2);
        }

    return P;
}

/**
 * @brief Frees the set of p_{i,j}.
 *
 * @param X The set of p_{i,j}.
 * @param H_t The number of heterogenous edges.
 */
void freeP(Z3_ast **P, int C_H) {
    for (int i = 0; i < C_H; i++) free(P[i]), P[i] = NULL;
    free(P), P = NULL;
}

/**
 * @brief Makes the set of l_{j,h} from the graph.
 *
 * @param ctx The solver context.
 * @param g The graph as an EdgeConGraph.
 * @param C_H The number of homogeneous components.
 * @return The set of l_{j,h} created.
 */
Z3_ast **makeL(Z3_context ctx, EdgeConGraph g, int C_H) {
    Z3_ast **L = (Z3_ast **)malloc(C_H * sizeof(Z3_ast *));
    for (int i = 0; i < C_H; i++) L[i] = (Z3_ast *)malloc(C_H * sizeof(Z3_ast));

    for (int j = 0; j < C_H; j++)
        for (int h = 0; h < C_H; h++) {
            L[j][h] = getVariableLevelInSpanningTree(ctx, h, j);
            // printf("Creating l_%d,%d…\n", j, h);
        }

    return L;
}

/**
 * @brief Frees the set of l_{e,i}.
 *
 * @param L The set of l_{e,i}.
 * @param C_H The number of homogeneous components.
 */
void freeL(Z3_ast **L, int C_H) {
    for (int i = 0; i < C_H; i++) free(L[i]), L[i] = NULL;
    free(L), L = NULL;
}

/**
 * @brief Creates a formula for at least one literal set to TRUE.
 *
 * @param ctx The solver context.
 * @param X The set of x_{e,i}.
 * @param size The size of the set of x_{e,i}.
 * @return The final formula.
 */
Z3_ast atLeast(Z3_context ctx, Z3_ast *X, int size) {
    Z3_ast formula = Z3_mk_false(ctx);
    for (int i = 0; i < size; i++) {
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
Z3_ast atMost(Z3_context ctx, Z3_ast *X, int size) {
    Z3_ast formula = Z3_mk_true(ctx);
    for (int i = 0; i < size; i++)
        for (int j = i + 1; j < size; j++) {
            Z3_ast args[2] = {Z3_mk_not(ctx, X[i]), Z3_mk_not(ctx, X[j])};
            Z3_ast clause = Z3_mk_or(ctx, 2, args);

            args[0] = formula;
            args[1] = clause;
            formula = Z3_mk_and(ctx, 2, args);
        }
    return formula;
}

/**
 * @brief Creates a formula that ensures that each translator can only be
 associated to at most one edge, and each edge can only receive at most one
 translator.

 * @param ctx The solver context.
 * @param X The set of x_{e,i}.
 * @param H_t The number of heterogeneous edges.
 * @param N The max number of translators.
 * @return The final formula phi1.
 */
Z3_ast phi1(Z3_context ctx, Z3_ast **X, int H_t, int N) {
    Z3_ast left = Z3_mk_true(ctx);
    Z3_ast right = Z3_mk_true(ctx);

    Z3_ast Xi[H_t];
    for (int i = 0; i < N; i++)
    {
        // Retrieve every x_e,i with i settled (the column i)
        for (int e = 0; e < H_t; e++) Xi[e] = X[e][i];

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
 * @brief Creates a formula that ensures that
 * each homogeneous component has one and only one parent
 * except the root.
 * 
 * @param ctx The solver context.
 * @param P The p_{j,j'} set.
 * @param L The l_{j,h} set.
 * @param C_H The number of homogeneous components.
 * @return The final formula phi2.
 */
Z3_ast phi2(Z3_context ctx, Z3_ast **P, Z3_ast **L, int C_H) {
    Z3_ast phi2 = Z3_mk_true(ctx);
    for (int j = 0; j < C_H; j++)
    {
        Z3_ast args_and[2] = {atLeast(ctx, P[j], C_H), atMost(ctx, P[j], C_H)};
        Z3_ast args_or[2] = {L[j][0], Z3_mk_and(ctx, 2, args_and)};
        Z3_ast args_phi2[2] = {phi2, Z3_mk_or(ctx, 2, args_or)};
        phi2 = Z3_mk_and(ctx, 2, args_phi2);
    }
    return phi2;
}

/**
 * @brief Creates a formula that ensures that
 * each homogeneous components has exactly one level.
 * 
 * @param ctx The solver context.
 * @param L The l_{i,h} set. 
 * @param C_H The number of homogeneous components.
 * @return The final formula phi3.
 */
Z3_ast phi3(Z3_context ctx, Z3_ast **L, int C_H) {
    Z3_ast phi3 = Z3_mk_true(ctx);

    for (int i = 0; i < C_H; i++) 
    {
        Z3_ast args_and[2] = {atLeast(ctx, L[i], C_H), atMost(ctx, L[i], C_H)};
        Z3_ast args_phi3[2] = {phi3, Z3_mk_and(ctx, 2, args_and)};
        phi3 = Z3_mk_and(ctx, 2, args_phi3);
    }
    return phi3;
}

/**
 * @brief Creates a formula that ensures that
 * the tree has a depth stricly higher than cost.
 * 
 * @param ctx The solver context.
 * @param L The l_{i,h} set.
 * @param C_H The number of homogeneous components.
 * @param cost 
 * @return The final formula phi4.
 */
Z3_ast phi4(Z3_context ctx, Z3_ast **L, int C_H, int cost) {
    Z3_ast phi4 = Z3_mk_false(ctx);

    for(int i = 0; i < C_H; i++)
    {
        for(int j = cost; j < C_H; j++) {
            Z3_ast args_or[2] = {phi4, L[i][j]};
            phi4 = Z3_mk_or(ctx, 2, args_or);
        }
    }
    return phi4;
} 

/**
 * @brief Creates a formula that ensures that, for a given couple of homogeneous components,
 * there exists at least one edge between them, and one of them is a translator.
 * 
 * @param ctx The solver context.
 * @param graph The graph as an EdgeConGraph.
 * @param X The x_e,i set.
 * @param n The number of vertices.
 * @param H_t The number of heterogeneous edges.
 * @param N The number of translators.
 * @param j1 The (index of the) first homogeneous component
 * @param j2 The (index of the) second homogeneous component
 * @return The final phi5 formula.
 */
Z3_ast phi5(Z3_context ctx, EdgeConGraph graph, Z3_ast **X, int n, int H_t, int N, int j1, int j2) {
    Z3_ast phi5 = Z3_mk_false(ctx);

    int e = 0;
    for (int u = 0; u < n; u++)
            for (int v = n - 1; v > u; v--)
                if (isEdgeHeterogeneous(graph, u, v))
                {
                    if (isNodeInComponent(graph, u, j1) && isNodeInComponent(graph, v, j2))
                    {
                        Z3_ast args_and[2] = {atLeast(ctx, X[e], N), atMost(ctx, X[e], N)};
                        Z3_ast args_phi5[2] = {phi5, Z3_mk_and(ctx, 2, args_and)};
                        phi5 = Z3_mk_or(ctx, 2, args_phi5);
                    }
                    e++;
                }
    
    return phi5;
}

/**
 * @brief Creates a formula that ensures that, for a given couple of homogeneous components,
 * if X_j1 has a level of h, then X_j2 has a level of h-1.
 * 
 * @param ctx The solver context.
 * @param L The l_{j,h} set. 
 * @param C_H The number of homogeneous components.
 * @param j1 The (index of the) first homogeneous component
 * @param j2 The (index of the) second homogeneous component
 * @return The final formula phi6.
 */
Z3_ast phi6(Z3_context ctx, Z3_ast **L, int C_H, int j1, int j2) {
    Z3_ast phi6 = Z3_mk_false(ctx);
    
    for (int h = 1; h < C_H; h++)
    {
        Z3_ast args_or[2] = {phi6, Z3_mk_implies(ctx, L[j1][h], L[j2][h-1])};
        phi6 = Z3_mk_or(ctx, 2, args_or);
    }
        
    return phi6;
}

/**
 * @brief Creates a formula that ensures that, for a given couple of homogeneous components,
 * if X_j2 is parent of X_j1 then phi5 and phi6 are guaranteed.
 * 
 * @param ctx The solver context.
 * @param graph A EdgeConGraph.
 * @param X The x_{e,i} set. 
 * @param n The number of vertices.
 * @param H_t The number of heterogeneous components.
 * @param N The number of translators.
 * @param P The p_{j,j'} set.
 * @param L The l_{j,h} set.
 * @param C_H The number of homogeneous components.
 * @param j1 The (index of the) first homogeneous component
 * @param j2 The (index of the) second homogeneous component
 * @return The final formula phi7.
 * @pre graph must be an initialized EdgeConGraph with computed connected
 * components.
 */
Z3_ast phi7(Z3_context ctx, EdgeConGraph graph, Z3_ast **X, int n, int H_t, int N, Z3_ast **P,  Z3_ast **L, int C_H, int j1, int j2) {
    Z3_ast args_and[2] = {phi5(ctx, graph, X, n, H_t, N, j1, j2),
                          phi6(ctx, L, C_H, j1, j2)};
    return Z3_mk_implies(ctx, P[j1][j2], Z3_mk_and(ctx, 2, args_and));
}

/**
 * @brief Generates a SAT formula satisfiable if and only if there is a set of
 * translators of cost @p cost such that the graph admits a valid path between
 * any two nodes.
 *
 * @param ctx The solver context.
 * @param graph A EdgeConGraph.
 * @param cost The cost of the translator set.
 * @return The formula.
 * @pre graph must be an initialized EdgeConGraph with computed connected
 * components.
 */
Z3_ast EdgeConReduction(Z3_context ctx, EdgeConGraph graph, int cost) {
    // Z3_ast X[5];
    // for (int i = 0; i < 5; i++)
    //     X[i] = Z3_mk_true(ctx);
    // return atMost(ctx, X, 5);
    // ----------------------------

    // Initialize useful variables from graph
    int n = orderG(getGraph(graph));
    int H_t = getNumHeteregeneousEdges(graph);
    int C_H = getNumComponents(graph);
    int N = C_H - 1;
    printf("n = %d, H_t = %d, C_H = %d, N = %d\n", n, H_t, C_H, N);

    // Initialize formula variables
    Z3_ast **X = makeX(ctx, graph, n, H_t, N);
    Z3_ast **P = makeP(ctx, graph, C_H);
    Z3_ast **L = makeL(ctx, graph, C_H);

    // phi1
    // Z3_ast varphi1 = phi1(ctx, X, H_t, N);
    // printf("phi1:\n%s\n", Z3_ast_to_string(ctx, varphi1));
    // Z3_lbool isvarPhi1Sat = isFormulaSat(ctx, varphi1);
    // if (isvarPhi1Sat == Z3_L_TRUE)
    //     printf("phi1 is SAT\n");
    // else
    //     printf("phi1 is not SAT\n");

    // phi2
    // Z3_ast varphi2 = phi2(ctx, P, L, C_H);
    // printf("phi2:\n%s\n", Z3_ast_to_string(ctx, varphi2));
    // Z3_lbool isvarPhi2Sat = isFormulaSat(ctx, varphi2);
    // if (isvarPhi2Sat == Z3_L_TRUE)
    //     printf("phi2 is SAT\n");
    // else
    //     printf("phi2 is not SAT\n");

    // phi3
    // Z3_ast varphi3 = phi3(ctx, L, C_H);
    // printf("phi3:\n%s\n", Z3_ast_to_string(ctx, varphi3));
    // Z3_lbool isvarPhi3Sat = isFormulaSat(ctx, varphi3);
    // if (isvarPhi3Sat == Z3_L_TRUE)
    //     printf("phi3 is SAT\n");
    // else
    //     printf("phi3 is not SAT\n");

    // phi4
    // Z3_ast varphi4 = phi4(ctx, L, C_H, cost);
    // printf("phi4:\n%s\n", Z3_ast_to_string(ctx, varphi4));
    // Z3_lbool isvarPhi4Sat = isFormulaSat(ctx, varphi4);
    // if (isvarPhi4Sat == Z3_L_TRUE)
    //     printf("phi4 is SAT\n");
    // else
    //     printf("phi4 is not SAT\n");

    // phi5
    // int j1 = 0;
    // int j2 = 1;
    // Z3_ast varphi5 = phi5(ctx, graph, X, n, H_t, N, j1, j2);
    // printf("phi5:\n%s\n", Z3_ast_to_string(ctx, varphi5));
    // Z3_lbool isvarPhi5Sat = isFormulaSat(ctx, varphi5);
    // if (isvarPhi5Sat == Z3_L_TRUE)
    //     printf("phi5 is SAT\n");
    // else
    //     printf("phi5 is not SAT\n");

    // phi6
    // int j1 = 0;
    // int j2 = 1;
    // Z3_ast varphi6 = phi6(ctx, L, C_H, j1, j2);
    // printf("phi6:\n%s\n", Z3_ast_to_string(ctx, varphi6));
    // Z3_lbool isvarPhi6Sat = isFormulaSat(ctx, varphi6);
    // if (isvarPhi6Sat == Z3_L_TRUE)
    //     printf("phi6 is SAT\n");
    // else
    //     printf("phi6 is not SAT\n");

    // phi7
    // int j1 = 0;
    // int j2 = 1;
    // Z3_ast varphi7 = phi7(ctx, graph, X, n, H_t, N, P, L, C_H, j1, j2);
    // printf("phi7:\n%s\n", Z3_ast_to_string(ctx, varphi7));
    // Z3_lbool isvarPhi7Sat = isFormulaSat(ctx, varphi7);
    // if (isvarPhi7Sat == Z3_L_TRUE)
    //     printf("phi7 is SAT\n");
    // else
    //     printf("phi7 is not SAT\n");

    // Free and return
    freeX(X, H_t);
    freeP(P, C_H);
    freeL(L, C_H);
    return Z3_mk_true(ctx); // TODO…
}

/**
 * @brief Gets the translator set from a model and adds it to the EdgeConGraph.
 * It also computes the homogeneous components taking into account the
 * translators (this number should be one).
 *
 * @param ctx The solver context.
 * @param model A variable assignment.
 * @param graph A EdgeConGraph.
 *
 * @pre @p model must be a valid model.
 * @pre @p graph must be a valid EdgeConGraph with no translators (or at least
 * the original number of homogeneous components).
 */
void getTranslatorSetFromModel(Z3_context ctx, Z3_model model,
                               EdgeConGraph graph) {
    return;
}
