/**
 * @file Parsing.h
 * @author Vincent Penelle (vincent.penelle@u-bordeaux.fr)
 * @brief Contains function to convert graphviz files to the GraphList structure. Needs first to compile Lexer.l and Parser.y with bison.
 * @version 1
 * @date 2019-07-31
 * 
 * @copyright Creative Commons
 * 
 */


#ifndef COCA_PARSING_H_
#define COCA_PARSING_H_

#include "GraphList.h"
#include "Graph.h"


/**
 * @brief Parses a file and return the Graph described by it. If the file with the name given in argument does not exists, it displays an error message and exits the program.
 * 
 * @param toRead the name of a file in graphviz format.
 * @return GraphList The parsed GraphList.
 * @pre @p toRead must be an existing file in graphviz format.
 */
Graph getGraphFromFile(char *toRead);


#endif