#include<stdlib.h>
#include "chunk_list.h"
//#include "dedupdef.h"

typedef struct iterator{
  List * list;
  Node * index;
}Iterator;

Iterator * init_iterator(List * list);
chunk_t * next(Iterator * iter);
Node * next_node(Iterator * iter);
void reset(Iterator * iter);
int hasNext(Iterator * iter);
void destroy_iterator(Iterator * iter);
