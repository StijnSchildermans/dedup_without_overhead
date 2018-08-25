#include "iterator.h"

Iterator * init_iterator(List * list){
  Iterator * iter = malloc(sizeof(Iterator));
  iter->list = list;
  iter->index = NULL;
  return iter;
}

Node * next_node(Iterator * iter){
  //printf("In next_node\n");
  Node * n = iter->index;
  if (n == NULL) {
    n = iter->list->head;
    if (n == NULL || n->data==NULL) return NULL;
    else{
      iter->index = n;
      return n;
    }
  }
  Node * nn = n->next;
  if(nn == NULL || nn->data==NULL) return NULL;
  iter->index = nn;
  //printf("Uit next_node\n");
  return nn;
}
chunk_t * next(Iterator * iter){
  //printf("Voor next_node\n");
  Node * n = next_node(iter);
  //printf("Next_node = %p\n",n);
  if (n == NULL || n->data == NULL) return NULL;
  return n->data;
};
void reset(Iterator * iter){
    iter->index = NULL;
}
int hasNext(Iterator * iter){
  return (iter->index != iter->list->tail);
}
void destroy_iterator(Iterator * iter){
  free(iter);
}
