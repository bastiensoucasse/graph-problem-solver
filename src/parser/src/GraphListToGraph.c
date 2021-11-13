#include "GraphListToGraph.h"
#include "EdgeList.h"
#include "NodeList.h"
#include <stdlib.h>
#include <string.h>

/*
 * @brief Auxilary function to determine the index of a node name.
 * 
 * @param list the array of node to search in.
 * @param size the size of list.
 * @param target the node to search
 * @return int the index of target in list.
 */
int findNode(char **list,int size,char *target){
	for(int i = 0; i<size;i++){
		if(strcmp(list[i],target)==0) return i;
	}
	return -1;
}

Graph createGraph(GraphList source){
	Graph res;

	res.numNodes=0;
	res.numEdges=0;
	SNodeList *explore = source.nodes;

	int count = 0;
	while(explore!=NULL){
		count++;
		explore = explore->next;
	}
	res.numNodes=count;

	//printf("nodes: %d\n",count);

	res.edges = (bool *)malloc(res.numNodes*res.numNodes*sizeof(bool));
	res.nodes = (char **)malloc(res.numNodes*sizeof(char*));

	//Ajout pour les automates.
	res.initial = (bool *)malloc(res.numNodes*sizeof(bool));
	res.final = (bool *)malloc(res.numNodes*sizeof(bool));

	count=0;
	explore = source.nodes;
	//Couleurs
	int numCol = 0;
	char **colorTab = (char**)malloc(res.numNodes*sizeof(char*));


	while(explore != NULL){
		res.nodes[count] = (char *)malloc((strlen(explore->node)+1)*sizeof(char));
		strcpy(res.nodes[count],explore->node);

		//Pour les automates.
		res.initial[count] = explore->initial;
		res.final[count] = explore->final;

		//Couleurs
		colorTab[count] = explore->color;
		if(colorTab[count] == NULL) colorTab[count] = "";
		numCol++;
		for(int i = 0; i < count; i++){
			if(strcmp(colorTab[i],colorTab[count])==0) {numCol--;break;}
		}

		count++;
		explore = explore->next;
	}


	//Couleurs
	res.numColor = numCol;
	res.colorNames = (char**)malloc((numCol+4)*sizeof(char*));
	res.color = (int*)malloc(res.numNodes*sizeof(int));
	int maxCurrentCol=0;
	for(int i = 0;i<res.numNodes;i++){
		int currentCol=0;
		while(currentCol<maxCurrentCol){
			if(strcmp(colorTab[i],res.colorNames[currentCol])==0) { break;}
			currentCol++;
		}
		res.color[i] = currentCol;
		if(maxCurrentCol == currentCol){
			maxCurrentCol++;
			res.colorNames[currentCol] = (char*)malloc((strlen(colorTab[i])+1)*sizeof(char));
			strcpy(res.colorNames[currentCol],colorTab[i]);
		}
	}

	free(colorTab);

	/*for(int i = 0; i<res.numNodes;i++){
		printf("%s,",res.nodes[i]);
	}
	printf("\n");
	*/

	for(int i = 0; i<res.numNodes;i++)
		for(int j = 0; j<res.numNodes;j++)
			res.edges[i*res.numNodes+j] = false;

	SEdgeList *exploreBis = source.edges;
	while(exploreBis != NULL){
		int n1,n2;
		n1=findNode(res.nodes,res.numNodes,exploreBis->node1);
		n2=findNode(res.nodes,res.numNodes,exploreBis->node2);
		res.edges[n1*res.numNodes+n2] = true;
		if(!source.directed) res.edges[n2*res.numNodes+n1] = true;
		exploreBis = exploreBis->next;
		res.numEdges++;
	}

	return res;
}