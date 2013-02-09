
#ifndef __POINT_LIST_H__
#define __POINT_LIST_H__

#include <stdio.h>

typedef struct PointNode {
	int col;
	int row;
	struct PointNode *next;
} PointList_t;

PointList_t *createList(int col, int row);
PointList_t *appendPoint(PointList_t * const head, int col, int row);
void destroyList(PointList_t *head);
PointList_t *findNearestPoint(PointList_t * const head, int col, int row);
void printList(FILE *fp, PointList_t * const head, const char* const sep);
#endif
