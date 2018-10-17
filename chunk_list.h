#ifndef LINKEDLIST_HEADER
#define LINKEDLIST_HEADER

#include <stdio.h>
#include <stdlib.h>
#include "dedupdef.h"
//#include "iterator.h"

typedef struct node{
  chunk_t * data;
  struct node * next;
  int allocated;
  char used;
} Node;

typedef struct list {
  //The first node
  Node * head;
  //The last non-empty node
  Node * tail;
  int length;
  //Node * empty;
} List;

List * emptylist();
void add(chunk_t * elem, List * list);
List ** split(int n, List * list);
List ** split_mod(int n, List * list);
List * merge(List * l1, List * l2);
List ** zip_split(int n, List ** lists);

#endif
