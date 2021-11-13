#include "Parsing.h"
#include "Parser.h"
#include "Lexer.h"
#include "GraphListToGraph.h"

int yyparse(GraphList *expression, yyscan_t scanner);
 

/**
 * @brief Parses a string and return the GraphList described by it.
 * 
 * @param expr A string in graphviz format.
 * @return GraphList The parsed GraphList.
 */
GraphList getGraphList(const char *expr)
{
    GraphList expression;
    yyscan_t scanner;
    YY_BUFFER_STATE state;

    expression.nodes=NULL;
    expression.edges=NULL;
 
    if (yylex_init(&scanner)) {
        /* could not initialize */
        printf("Error initialization\n");
        return expression;
    }
 
    state = yy_scan_string(expr, scanner);

    if (yyparse(&expression, scanner)) {
        /* error parsing */
        printf("Error parsing\n");
        return expression;
    }

    yy_delete_buffer(state, scanner);
 
    yylex_destroy(scanner);
 
    return expression;
}

/**
 * @brief Parses a file and return the GraphList described by it.
 * 
 * @param toRead A file in graphviz format.
 * @return GraphList The parsed GraphList.
 */
GraphList getGraphListFromFile(FILE *toRead)
{
    GraphList expression;
    yyscan_t scanner;
    YY_BUFFER_STATE state;

    expression.nodes=NULL;
    expression.edges=NULL;
 
    if (yylex_init(&scanner)) {
        /* could not initialize */
        printf("Error initialization\n");
        return expression;
    }

    state = yy_create_buffer(toRead, YY_BUF_SIZE, scanner);
    yy_switch_to_buffer(state,scanner);

    if (yyparse(&expression, scanner)) {
        /* error parsing */
        printf("Error parsing\n");
        return expression;
    }

    yy_delete_buffer(state, scanner);

    yylex_destroy(scanner);

	fclose(toRead);

    return expression;
}

Graph getGraphFromFile(char *toRead){
    FILE* file = fopen(toRead,"r");
    if(file == NULL){
        printf("file %s does not exist. Exiting.\n",toRead);
        exit(-1);
    }
    GraphList e = getGraphListFromFile(file);
    Graph graph = createGraph(e);
    deleteExpression(e.edges);
    deleteNodeList(e.nodes);
    return graph;
}
