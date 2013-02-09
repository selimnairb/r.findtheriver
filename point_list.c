#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "point_list.h"

PointList_t *createList(int col, int row) {
	PointList_t *head = (PointList_t *)malloc(sizeof(PointList_t));
	head->col = col;
	head->row = row;
	head->next = NULL;
	return head;
}

PointList_t *appendPoint(PointList_t * const head, int col, int row) {
	assert( NULL != head);

	PointList_t *tmp = head;
	while ( NULL != tmp->next ) tmp = tmp->next;
	tmp->next = (PointList_t *)malloc(sizeof(PointList_t));
	tmp->next->col = col;
	tmp->next->row = row;
	tmp->next->next = NULL;
	return tmp->next;
}

void destroyList(PointList_t *head) {
	if ( NULL != head ) {
		PointList_t *curr = head;
		PointList_t *next = curr->next;
		while ( NULL != next ) {
			free(curr);
			curr = next;
			next = curr->next;
		}
		free(curr);
	}
}

PointList_t *findNearestPoint(PointList_t * const head, int col, int row) {

	PointList_t *nearest = NULL;
	double tmpDistance, minDistance = HUGE_VAL;

	PointList_t *curr = head;
	while ( NULL != curr ) {
		tmpDistance = sqrt( pow(col - curr->col, 2 ) + pow(row - curr->row, 2 ) );
		if ( tmpDistance < minDistance ) {
			minDistance = tmpDistance;
			nearest = curr;
		}
		curr = curr->next;
	}
	return nearest;
}

void printList(FILE *fp, PointList_t * const head, const char* const sep) {
	PointList_t *curr = head;
	while ( NULL != curr ) {
		fprintf(fp, "%d%s%d\n", curr->col, sep, curr->row);
		curr = curr->next;
	}
}
