#include<stdlib.h>
#include<string.h>
#include "chunk_list.h"
#include "iterator.h"

void createNNodes(int n, List * list){
  Node * newNodes = malloc(n * sizeof(Node));

  for (int i = 0; i < n-1; i++){
    newNodes[i].allocated = 0;
    newNodes[i].data = NULL;
    newNodes[i].next = &newNodes[i+1];
    newNodes[i].used = 1;
  }
  newNodes[0].allocated = n;
  newNodes[n-1].next = NULL;
  //No non-empty elemnts
  if(list->tail == NULL) {
    //No empty elements either
    if (list->head == NULL) list->head = newNodes;
    //Only empty elements
    else {
      Node * h = list->head;
      while(h->next != NULL) h = h->next;
      h->next = newNodes;
    }
  }
  else list->tail->next = newNodes;
}

void createNodes(List * list){
  int n;
  if (list->length < 16) n = 16;
  else if (list->length < 1024) n = list->length;
  else n = 1024;
  createNNodes(n,list);
}

List * emptylist(){
  List * list = malloc(sizeof(List));
  list->head = NULL;
  list->tail = NULL;
  list->length = 0;
  return list;
}

void add(chunk_t * elem, List * list){
  if (list->head == NULL) createNodes(list);
  if (list->head->data == NULL || list->tail == NULL){
    list->head->data = elem;
    list->tail = list->head;
  }
  else{
    if(list->tail->next == NULL) createNodes(list);
    list->tail = list->tail->next;
    list->tail->data = elem;
  }
  list->length++;
}
void add_node(Node * node, List * list){
  //Empty list
  if(list->head == NULL) {
    node->next = NULL;
    list->head = node;
    list->tail = node;

  }
  //List with only empty nodes
  else if (list->head->data == NULL){
    node->next = list->head;
    list->head = node;
    list->tail = node;
  }
  else{
    node->next = list->tail->next;
    list->tail->next = node;
    list->tail = node;
  }
  list->length++;
}
//Find the number of memory allocations for a list.
int numAllocs(List * list){
  int len = list->length;
  if (list == NULL) return 0;
  else if (len < 17) return 1;
  else if (len < 33) return 2;
  else if (len < 65) return 3;
  else if (len < 129) return 4;
  else if (len < 257) return 5;
  else if (len < 513) return 6;
  else if (len < 1025) return 7;
  else return 7 + ((len - 1024) / 1024) + ((len % 1024 == 0) ? 0:1);
}
//Splits a list in n sublists of sequential elements.
List ** split(int n, List * list){
  int size = (list->length)/n;
  List ** lists = malloc(n * sizeof(List*));
  for (int q = 0; q<n;q++) lists[q] = emptylist();
  Node * buffer;
  buffer = list->head;

  for(int i = 0; i< list->length; i++){
    int l;
    if(i/size < n) l = (i/size);
    else l = (n-1);
    Node * nn = buffer->next;
    add_node(buffer,lists[l]);
    buffer = nn;
  }
  free(list);

  return lists;
}
//Splits a list in n sublists with each m'th element of the sublist being the (n*m)'th element of the original list
List ** split_mod(int n, List * list){
  List ** lists = malloc(n * sizeof(List*));
  for (int q = 0; q<n;q++) lists[q] = emptylist();
  int i = 0;
  Node * buffer = list->head;

  int len = list->length;
  for(int j = 0; j < len; j++){
    Node * nn = buffer->next;
    add_node(buffer,lists[i]);
    buffer = nn;

    if(i == (n - 1)) i = 0;
    else i++;
  }
  return lists;
}
void merge_empty(List * l1, List * l2){
  Node * fempty = NULL;
  Node * lempty = NULL;

  //Initialize fempty;
  if (l1 == NULL) l1 = emptylist();
  //Only empty elements in l1
  else if(l1->head != NULL && l1->tail == NULL) fempty = l1->head;
  //No empty elements in l1
  if (l1->head == NULL || l1->tail->next == NULL){
    if (l2 == NULL || l2->head == NULL) return;
    else if(l2->tail == NULL) fempty = l2->head;
    else if(l2->tail->next == NULL) return;
    else fempty = l2->tail->next;
  }
  //All other cases
  else fempty = l1->tail->next;

  //Initialize lempty
  lempty = fempty;
  while(lempty->next != NULL) lempty = lempty->next;

  //merge lempty and first empty element of l2 if necessary
  if(l2 != NULL && l2->tail != NULL && l2->tail->next != NULL && fempty != l2->tail->next) lempty->next = l2->tail->next;
  else if(l2->head != NULL && l2->tail == NULL && fempty != l2->head) lempty->next = l2->head;

  //Do what is necessary to return l1 with the merged empty sections;
  if(l1->head == NULL || l1->tail == NULL) l1->head = fempty;
  else l1->tail->next = fempty;

  //Remove empty nodes from l2 if necessary
  if (l2 != NULL && l2->tail != NULL) l2->tail->next = NULL;
  else if (l2 != NULL && l2->head != NULL && l2->tail == NULL) l2->head = NULL;
}

List * merge(List * l1, List * l2){
  if (l1 == NULL) return l2;
  else if (l1->head == NULL){
    free(l1);
    return l2;
  }
  else if (l2 == NULL) return l1;
  else if (l2->head == NULL){
    free(l2);
    return l1;
  }
  merge_empty(l1,l2);
  l2->tail->next = l1->tail->next;
  l1->tail->next = l2->head;
  l1->tail = l2->tail;
  l1->length += l2->length;
  free(l2);
  return l1;
}

//Zips n lists that were split using split_mod.
List ** zip_split(int n, List ** lists){
  List ** output = malloc(n*sizeof(List *));
  for (int i = 0; i < n; i++){
    output[i] = emptylist();
    merge_empty(output[i],lists[i]);
  }
  Node * buffers[n];
  int i;
  for(i=0;i<n;i++)buffers[i] = lists[i]->head;
  int len = lists[0]->length;
  int out_list = 0;
  int count = 0;

  for (i = 0; i<len;i++){
    for (int j = 0; j<n;j++){
      if(buffers[j]!= NULL){
        Node * nnn = buffers[j];
        buffers[j] = nnn->next;
        add_node(nnn,output[out_list]);
        count++;
        if((out_list < n-1) && (count >= len && ((j < n-1 && buffers[j+1] != NULL && buffers[j+1]->data->sequence.l1num != nnn->data->sequence.l1num)
              || (j == n-1 && buffers[0] != NULL && buffers[0]->data->sequence.l1num != nnn->data->sequence.l1num)))){
          out_list++;
          count = 0;
        }
      }
    }
  }
  for(i=0;i<n;i++) free(lists[i]);
  free(lists);
  return output;
}
