/*
 * @file Graph.c
 * @author Vincent Penelle (vincent.penelle@u-bordeaux.fr)
 * @brief  Structure to store a graph statically, and to access to its informations easily.
  		   Includes unary automata informations (initial and final nodes).
		   Does not contain acces functions for now.
 * @version 0.5
 * @date 2018-11-18, 2019-07-23
 * 
 * @copyright Creative Commons.
 * 
 */


#include "Graph.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

void printGraph(Graph graph){
	//Colors
	printf("There are %d colors: ",graph.numColor);
	for(int i = 0; i<graph.numColor;i++) printf("%s(%d) ",graph.colorNames[i],i);
	printf("\nnodes:\n");
	for(int i = 0; i<graph.numNodes;i++) printf("%d : %s [%s], ",i,graph.nodes[i],graph.colorNames[graph.color[i]]);

	//If uncolored
	/*	
	printf("\nnodes:\n");
	for(int i = 0; i<graph.numNodes;i++) printf("%d : %s, ",i,graph.nodes[i]);
	*/

	//Automata.
	//printf("\ninitial:\n");
	//for(int i = 0; i<graph.numNodes;i++) printf("%d ",graph.initial[i]);
	//printf("\nfinal:\n");
	//for(int i = 0; i<graph.numNodes;i++) printf("%d ",graph.final[i]);


	printf("\nEdges:\n");
	for(int i = 0; i<graph.numNodes;i++){
		for(int j = 0; j<graph.numNodes;j++){
			printf("%d ",graph.edges[i*graph.numNodes+j]);
		}
		printf("\n");
	}
}

Graph copyGraph(Graph graph){
	Graph copy;
	copy.numNodes = graph.numNodes;
	copy.numEdges = graph.numEdges;
	copy.nodes = (char**)malloc(copy.numNodes*sizeof(char*));
	copy.color = (int*)malloc(copy.numNodes*sizeof(int));
	for(int i = 0; i < copy.numNodes; i++){
		copy.nodes[i] = (char*)malloc(strlen(graph.nodes[i])*sizeof(char));
		strcpy(copy.nodes[i],graph.nodes[i]);
		copy.color[i] = graph.color[i];
	}
	copy.edges = (bool*)malloc(copy.numNodes*copy.numNodes*sizeof(bool));

	for(int i = 0; i < copy.numNodes*copy.numNodes; i++) copy.edges[i] = graph.edges[i];
	copy.numColor = graph.numColor;
	copy.colorNames = (char**)malloc(copy.numColor*sizeof(char*));
	for(int i = 0; i < copy.numColor; i++){
		copy.colorNames[i] = (char*)malloc(strlen(graph.colorNames[i])*sizeof(char));
		strcpy(copy.colorNames[i],graph.colorNames[i]);
	}

	copy.initial = (bool*)malloc(copy.numNodes*sizeof(bool));
	for(int i = 0; i < copy.numNodes; i++) copy.initial[i] = graph.initial[i];
	copy.final = (bool*)malloc(copy.numNodes*sizeof(bool));
	for(int i = 0; i < copy.numNodes; i++) copy.final[i] = graph.final[i];
	return copy;
}

void deleteGraph(Graph graph){
	if(graph.edges!=NULL) free(graph.edges);
	if(graph.nodes!=NULL){
		for(int i = 0; i<graph.numNodes; i++) {
			if(graph.nodes[i]!=NULL) free(graph.nodes[i]);
		}
		free(graph.nodes);
	}
	//Pour les automates.
	if(graph.initial!=NULL) free(graph.initial);
	if(graph.final!=NULL) free(graph.final);

	//Couleurs
	if(graph.color != NULL) free(graph.color);
	
	if(graph.colorNames != NULL){
		for(int i = 0; i<graph.numColor; i++){
			if(graph.colorNames[i]!=NULL) free(graph.colorNames[i]);
		}
		free(graph.colorNames);
	}
	graph.numEdges=0;
	graph.numNodes=0;
}

int orderG(Graph graph){
	return graph.numNodes;
}

int sizeG(Graph graph){
	return graph.numEdges;
}

bool isEdge(Graph graph, int source, int target){
	return graph.edges[(source*graph.numNodes)+target];
}

bool isSource(Graph graph, int node){
	return graph.initial[node];
}

bool isTarget(Graph graph, int node){
	return graph.final[node];
}

char* getNodeName(Graph graph, int node){
	return graph.nodes[node];
}

int getColor(Graph graph, int node){
	return graph.color[node];
}

char* getColorString(Graph graph, int col){
	return graph.colorNames[col];
}

int getNumColor(Graph graph){
	return graph.numColor;
}
