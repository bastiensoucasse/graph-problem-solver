/**
 * @file main.c
 * @author Vincent Penelle (vincent.penelle@u-bordeaux.fr)
 * @brief  File putting the whole program together.
 * @version 1
 * @date 2019-07-22
 * 
 * @copyright Creative Commons.
 * 
 */

#include "Graph.h"
#include "Parsing.h"
#include "EdgeConReduction.h"
#include "Z3Tools.h"
#include "Parser.h"
#include "EdgeConGraph.h"
#include "EdgeConResolution.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

void printModel(Z3_context ctx, Z3_model model, EdgeConGraph biGraph, int numComponent)
{
    printf("\nPrinting tree over components (for understanding your formula) -- refer to display of the graph with -v to see which number corresponds to which component:\n");

    for (int child = 0; child < numComponent; child++)
    {
        for (int parent = 0; parent < numComponent; parent++)
        {
            Z3_ast var = getVariableParent(ctx, child, parent);
            if (valueOfVarInModel(ctx, model, var))
                printf("%d is a child of %d\n", child, parent);
        }
    }

    for (int comp = 0; comp < numComponent; comp++)
    {
        for (int level = 0; level < numComponent; level++)
        {
            Z3_ast var = getVariableLevelInSpanningTree(ctx, level, comp);
            if (valueOfVarInModel(ctx, model, var))
                printf("%d is of level %d\n", comp, level);
        }
    }

    printf("\n");
}

void usage()
{
    printf("Use: graphProblemSolver [options] file\n");
    printf(" file should contain a colored graph in dot format.\n The program decides if there exists a translator set for the graph in input.\n Can apply a brute force algorithm or a reduction to SAT. In the latter case, a weight for the expected solution must be given.\n Can display the result both on the command line or in dot format.\n");
    printf("Options: \n");
    printf(" -h         Displays this help\n");
    printf(" -v         Activate verbose mode (displays parsed graphs)\n");
    printf(" -B         Solves the problem using the brute force algorithm\n");
    printf(" -R COST    Solves the problem using a reduction and determines if for all possible translator sets, all nodes can communicate with cost at most COST\n");
    printf(" -F         Displays the formula computed (obviously not in this version, but you should really display it in your code). Only active if -R is active\n");
    printf(" -t         Displays the translator set found [if not present, only displays the existence of the set].\n");
    printf(" -M         If there is a solution to the reduction, displays the tree obtained over the homogeneous components. Only has an effect if -R is present\n");
    printf(" -f         Writes the result with colors in a .dot file. See next option for the name. These files will be produced in the folder 'sol'.\n");
    printf(" -o NAME    Writes the output graph in \"NAME_Brute.dot\" or \"NAME_SAT.dot\" depending of the algorithm used and the formula in \"NAME.formul\". [if not present: \"result_SAT.dot\", \"result_Brute.dot\" and \"result.formul\"]\n");
}

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        usage();
        return 0;
    }

    bool verbose = false;
    bool displayTerminal = false;
    bool outputFile = false;
    bool printformula = false;
    bool bruteForce = false;
    bool reduction = false;
    bool displayModel = false;
    int size = 0;
    char *solutionName = "default";
    char *realArgs[argc];
    int numArgs = 0;

    int option;

    while ((option = getopt(argc, argv, ":hvFBMGR:tfo:")) != -1)
    {
        switch (option)
        {
        case 'h':
            usage();
            return EXIT_SUCCESS;
        case 'v':
            verbose = true;
            break;
        case 'B':
            bruteForce = true;
            break;
        case 'R':
            reduction = true;
            size = atoi(optarg);
            break;
        case 'F':
            //printf("Don't insist, I'm not showing you the solution of the assignment yet!\n");
            printformula = true;
            break;
        case 't':
            displayTerminal = true;
            break;
        case 'M':
            displayModel = true;
            break;
        case 'f':
            outputFile = true;
            break;
        case 'o':
            solutionName = optarg;
            break;
        case '?':
            printf("unknown option: %c\n", optopt);
            break;
        }
    }

    if (argc - optind < 1)
    {
        printf("No argument given. Exiting.\n");
        return 0;
    }

    Graph graph = getGraphFromFile(argv[optind]);

    if (verbose)
        printGraph(graph);

    EdgeConGraph biGraph = initializeGraph(graph);

    if (verbose)
    {
        printf("Connected Components: \n");
        for (int component = 0; component < getNumComponents(biGraph); component++)
        {
            for (int node = 0; node < orderG(graph); node++)
            {
                if (isNodeInComponent(biGraph, node, component))
                    printf("- %s ", getNodeName(graph, node));
            }
            printf("\n");
        }
    }

    if (bruteForce)
    {
        printf("\n*******************\n*** Brute Force ***\n*******************\n\n");
        clock_t start = clock();
        int res = BruteForceEdgeCon(biGraph);
        double end = (double)(clock() - start) / CLOCKS_PER_SEC;
        if (res >= 0)
        {
            printf("Brute force computed the solution in %g seconds: All possible assignations of translators allow nodes to communicate with at most %d translators on the path\n", end, res);
            if (displayTerminal || outputFile)
                printf("A translator set reaching that bound has been computed\n");
            if (displayTerminal)
                printTranslator(biGraph);
            int length = strlen(solutionName) + 12;
            char nameFile[length];
            snprintf(nameFile, length, "%s_Brute", solutionName);
            if (outputFile)
            {
                createDotOfEdgeConGraph(biGraph, nameFile);
                printf("Solution printed in sol/%s.dot.\n", nameFile);
            }
        }
        else
            printf("No solution found by Brute Force in %g seconds\n", end);
        resetTranslator(biGraph);
    }

    if (reduction)
    {
        printf("\n************************\n*** Reduction to SAT ***\n************************\n\n");

        if (size <= 0)
        {
            printf("No weight given, or weight given less than 0, I refuse to compute the formula for it!\n");
        }
        else
        {
            Z3_context ctx = makeContext();

            clock_t start = clock();

            Z3_ast formula;
            formula = EdgeConReduction(ctx, biGraph, size);

            clock_t timeFormula = clock();

            printf("formula computed in %g seconds\n", (double)(timeFormula - start) / CLOCKS_PER_SEC);

            if (printformula)
            {
#ifndef SUBJECT
                struct stat st = {0};
                if (stat("./sol", &st) == -1)
                    mkdir("./sol", 0777);
                int length = strlen(solutionName) + 13;
                char nameFile[length];
                snprintf(nameFile, length, "sol/%s.formula", solutionName);
                FILE *file = fopen(nameFile, "w");
                fprintf(file, "%s\n", Z3_ast_to_string(ctx, formula));
                fclose(file);
                printf("Formula printed in sol/%s.formula\n", solutionName);
#else
                printf("Nah, I'm not displaying the formula in the given executable\n");
#endif
            }

            Z3_model model;
            Z3_lbool isSat = solveFormula(ctx, formula, &model);

            clock_t timeSat = clock();

            printf("solution computed in %g seconds\n", (double)(timeSat - timeFormula) / CLOCKS_PER_SEC);

            switch (isSat)
            {
            case Z3_L_FALSE:
                printf("All possible translator sets allow all nodes to communicate with cost at most %d.\n", size);
                break;

            case Z3_L_UNDEF:
                printf("Not able to decide if there is a translator set forcing some nodes to communicate with cost bigger than %d.\n", size);
                break;

            case Z3_L_TRUE:
                printf("There is a translator set forcing some node to communicate with cost bigger than %d.\n", size);

                int numComponent = getNumComponents(biGraph);
                if (displayTerminal || outputFile || displayModel)
                    getTranslatorSetFromModel(ctx, model, biGraph);

                if (displayModel)
                    printModel(ctx, model, biGraph, numComponent);

                if (displayTerminal)
                {
                    printTranslator(biGraph);
                }

                if (outputFile)
                {
                    int length = strlen(solutionName) + 12;
                    char nameFile[length];
                    snprintf(nameFile, length, "%s_Sat", solutionName);
                    createDotOfEdgeConGraph(biGraph, nameFile);
                    printf("Solution printed in sol/%s.dot.\n", nameFile);
                }

                break;
            }

            Z3_del_context(ctx);
        }
    }

    deleteEdgeConGraph(biGraph);

    deleteGraph(graph);

    return 0;
}
