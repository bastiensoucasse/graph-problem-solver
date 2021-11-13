%{
/**
 * @file Lexer.l
 * @author Vincent Penelle (vincent.penelle@u-bordeaux.fr)
 * @brief  Lexer for a graphviz parser, intended to serve for a master's project at University of Bordeaux. Adapted from gvizparse by Nikolaos Kavvadias (https://github.com/nkkav/gvizparse v1.0.0).
 * @version 1
 * @date 2019-07-24
 * 
 * @copyright Creative Commons.
 * 
 */
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include "GraphList.h"
#include "Parser.h"

 
 #define MAXSIZE  1024
char any_name_s[MAXSIZE];

%}

/* %option outfile="Lexer.c" header-file="Lexer.h"  //for normal make.*/
%option warn 
 
%option reentrant noyywrap never-interactive nounistd
%option bison-bridge


ws	  [ \t\n]
anum  [a-zA-Z0-9_][a-zA-Z0-9_.]*
digit [0-9]
A	  [aA]
B	  [bB]
C	  [cC]
D	  [dD]
E	  [eE]
F	  [fF]
G	  [gG]
H	  [hH]
I	  [iI]
J	  [jJ]
K	  [kK]
L	  [lL]
M	  [mM]
N	  [nN]
O	  [oO]
P	  [pP]
Q	  [qQ]
R	  [rR]
S	  [sS]
T	  [tT]
U	  [uU]
V	  [vV]
W	  [wW]
X	  [xX]
Y	  [yY]
Z	  [zZ]

%%
"//".*          { }
\"(\\.|[^\\"])*\"	{ strcpy(any_name_s,yytext);
                  yylval->name = any_name_s;
                  return(T_STRING); }
{ws}+		        ;
"["             { return(T_LBRACKET); }
"]"             { return(T_RBRACKET); }
"("             { return(T_LPAREN); }
")"             { return(T_RPAREN); }
"{"             { return(T_LBRACE); }
"}"             { return(T_RBRACE); }
","             { return(T_COMMA); }
":"             { return(T_COLON); }
";"             { return(T_SEMI); }
"->"            { return(T_DEDGE); }
"--"            { return(T_UEDGE); }
"="             { return(T_EQ); }
{D}{I}{G}{R}{A}{P}{H}     { return(T_DIGRAPH); }
{G}{R}{A}{P}{H}           { return(T_GRAPH); }
{S}{U}{B}{G}{R}{A}{P}{H}  { return(T_SUBGRAPH); }
{A}{T}          { return(T_AT); }
{S}{T}{R}{I}{C}{T}        { return(T_STRICT); }
{N}{O}{D}{E}    { return(T_NODE); }
{E}{D}{G}{E}    { return(T_EDGE); }
{anum}          { //any_name_s = malloc(sizeof(yytext)+1);
                  strcpy(any_name_s, yytext);
                  yylval->name = any_name_s;
                  return(T_ID); }

%%


