#include<stdlib.h>
#include<string.h>
#include "chunk_list.h"
#include "iterator.h"


/*Node * createnode(chunk_t * data){
  Node * newNode = malloc(sizeof(Node));
  newNode->data = data;
  newNode->next = NULL;
  return newNode;
}*/

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
  //printf("Uit createNodes\n");
  //list->empty = newNodes;
}

void createNodes(List * list){
  int n;
  //printf("Added nodes\n");
  //printf("In createNodes\n");
  if (list->length < 16) n = 16;
  else if (list->length < 1024) n = list->length;
  else n = 1024;
  //printf("N = %d\n",n);
  createNNodes(n,list);
}

List * emptylist(){
  List * list = malloc(sizeof(List));
  list->head = NULL;
  list->tail = NULL;
  list->length = 0;
  //createNodes(list);
  return list;
}
List * nonemptylist(int n){
  List * list = emptylist();
  createNNodes(n,list);
  return list;
}

void add(chunk_t * elem, List * list){
  //printf("In add\n");
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
  //printf("Uit add\n");
}
void add_node(Node * node, List * list){
  //Empty list
  if(list->head == NULL) {
    //createNodes(list);
    //node->next = list->head->next;
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
void add_nodes(Node * head, Node * tail, List * list){
  int len = 1;
  Node * n = head;
  while (n != tail){
    n = n->next;
    len++;
  }
  //List is empty
  if (list->head == NULL) {
    list->head = head;
    list->tail = tail;
  }
  //List only has empty nodes
  else if (list->head->data == NULL){
    tail->next = list->head;
    list->head = head;
    list->tail = tail;
  }
  else{
    tail->next = list->tail->next;
    list->tail->next = head;
    list->tail = tail;
  }
  list->length += len;
}
/*int length(List * list){
  int i = 1;
  Iterator * iter = init_iterator(list);
  while(next_node(iter)!=list->tail) i++;
  destroy_iterator(iter);
  return i;
}*/

chunk_t * get(int index, List * list){
  int i = 0;
  Node * n;
  if (list->head == NULL) return NULL;
  else n = list->head;
  while (i<index){
    n = n->next;
    if (n==NULL) return NULL;
    i++;
  }
  return n->data;
}

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

void destroy(List * list){
  if (list == NULL) return;
  Node * current = list->head;
  Node * next = current;
  Node * allocated_nodes[numAllocs(list)];
  int allocated_index = 0;
  while(current != NULL){
    next = current->next;
    //if (current->data != NULL) free(current->data);
    if(current->allocated > 0){
      allocated_nodes[allocated_index] = current;
      allocated_index++;
    }
    current->used--;
    current = next;
  }
  for (int i = 0; i< allocated_index;i++){
    //Node * allocated_node = allocated_nodes[i];
    //char doubles = 0;
    //printf("Voor loop\n");
    //for (int j = 0; j<allocated_node->allocated; j++){
      //printf("Aantal allocated nodes: %d\n",allocated_node->allocated);
      //if (allocated_node[j].used > 0) doubles = 1;
    //}
    //printf("Na loop\n");
    //if (!doubles){
      //printf("Voor free\n");
      free(allocated_nodes[i]);
      //printf("Na free\n");
    //}
  }
  //printf("Voor free list\n");
  //free(list);
  //printf("Na free list\n");
}
void destroy_soft(List * list){
  Node * current = list->head;
  Node * next = current;
  //printf("Lengte lijst:%d\n",list->length);
  while(current != NULL){
    next = current->next;
    //if(current->allocated) {
      free(current);
      //printf("Node vrijgegeven\n");
    //}
    current = next;
    //printf("Node vrijgegeven\n");
  }
  free(list);
}
/*void destroy_empty(List * list){
  if (list->tail == NULL) return;
  Node * current = list->tail->next;
  Node * next = current;
  while(current != NULL){
    next = current->next;
    free(current);
    current = next;
  }
  free(list);
}*/
//Splits a list in n sublists of sequential elements.
List ** split(int n, List * list){
  //printf("In split\n");
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
  //Iterator * iter = init_iterator(list);
  //while(hasNext(iter)){
    //add(next(iter),lists[l]);
    //i++;
  //}//
  //destroy_soft(list);
  //destroy_iterator(iter);
  //destroy_empty(list);
  //printf("Uit split\n");
  free(list);

  return lists;
}
void compare(List * l1, List* l2){
  Iterator * iter1 = init_iterator(l1);
  Iterator* iter2;
  while(hasNext(iter1)){
      iter2 = init_iterator(l2);
      while(hasNext(iter2)){
        if(next(iter2) == next(iter1)) printf("Found double chunk!\n");
      }
      destroy_iterator(iter2);
  }
}
//Splits a list in n sublists
// with each m'th element of the sublist being the (n*m)'th e
//lement of the original list
List ** split_mod(int n, List * list){
  List ** lists = malloc(n * sizeof(List*));
  for (int q = 0; q<n;q++) lists[q] = emptylist();
  int i = 0;
  Node * buffer = list->head;

  int len = list->length;
  for(int j = 0; j < len; j++){
    //if(buffer == NULL)printf("Buffer = NULL, len = %d, j = %d\n",len,j);
    Node * nn = buffer->next;
    add_node(buffer,lists[i]);
    buffer = nn/*->next*/;

    if(i == (n - 1)) i = 0;
    else i++;
  }
  //for (int q = 0; q < n; q++) lists[q]->tail->next = NULL;
  //destroy_empty(list);
  return lists;
}
void check_empty(List * list){
  Node * n = list->tail->next;
  while (n!= NULL){
    if (n->data != NULL) printf("Data is niet NULL, l1num = %d\n",n->data->sequence.l1num);
    n = n->next;
  }
}
//Adds the empty nodes of l2 to the end of l1
void merge_empty(List * l1, List * l2){
  Node * fempty = NULL;
  Node * lempty = NULL;

  //Initialize fempty;
  if (l1 == NULL) l1 = emptylist();
  //Only empty elements in l1
  else if(l1->head != NULL && l1->tail == NULL) fempty = l1->head;
  //No empty elements in l1
  if (l1->head == NULL || l1->tail->next == NULL){
    //printf("l1->head == NULL\n");
    if (l2 == NULL || l2->head == NULL) return;
    else if(l2->tail == NULL) fempty = l2->head;
    else if(l2->tail->next == NULL) return;
    else {
       fempty = l2->tail->next;
       //printf("fempty = l2->tail->next\n");
     }
  }
  //All other cases
  else fempty = l1->tail->next;

  //printf("Initialized fempty: %p\n",fempty);

  //Initialize lempty
  lempty = fempty;
  while(lempty->next != NULL) lempty = lempty->next;

  //printf("Initialized lempty\n");

  //merge lempty and first empty element of l2 if necessary
  if(l2 != NULL && l2->tail != NULL && l2->tail->next != NULL && fempty != l2->tail->next){
    lempty->next = l2->tail->next;
  }
  else if(l2->head != NULL && l2->tail == NULL && fempty != l2->head) lempty->next = l2->head;

  //printf("Merged empty parts of l1 and l2\n");

  //Do what is necessary to return l1 with the merged empty sections;
  if(l1->head == NULL || l1->tail == NULL) l1->head = fempty;
  else l1->tail->next = fempty;

  //Remove empty nodes from l2 if necessary
  if (l2 != NULL && l2->tail != NULL) l2->tail->next = NULL;
  else if (l2 != NULL && l2->head != NULL && l2->tail == NULL) l2->head = NULL;
  //check_empty(l1);
}

void check_sequence(List * list){
  Iterator * iter = init_iterator(list);

  int l1num = 0;
  int l2num = 0;

  while(hasNext(iter)){
    chunk_t* c = next(iter);
    if (c->sequence.l1num < l1num || (c->sequence.l2num < l2num && c->sequence.l1num == l1num) /*|| c->sequence.l2num == l2num*/){
      printf("Error in sequence: old l1num = %d, old l2num = %d, new l1num = %d, new l2num = %d\n",l1num, l2num, c->sequence.l1num,c->sequence.l2num);
    }
    l1num = c->sequence.l1num;
    l2num = c->sequence.l2num;
  }

  destroy_iterator(iter);
}

void check_sequence_exact(List * list){
  Iterator * iter = init_iterator(list);

  chunk_t * c = next(iter);
  int l1num = c->sequence.l1num;
  int l2num = c->sequence.l2num;


  while(hasNext(iter)){
    c = next(iter);
    if (!((c->sequence.l1num == l1num && c->sequence.l2num == l2num+1)
       || (c->sequence.l1num == l1num+1 && c->sequence.l2num == 0))){
      printf("Error in sequence: old l1num = %d, old l2num = %d, new l1num = %d, new l2num = %d\n",l1num, l2num, c->sequence.l1num,c->sequence.l2num);
    }
    l1num = c->sequence.l1num;
    l2num = c->sequence.l2num;
  }

  destroy_iterator(iter);
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
  //free(l2);
  //Node * n = l1->head;
  /*for (int i = 1; i <= l1->length; i++){
    if(n == NULL) printf("Node %d is NULL\n",i);
    n = n->next;
  }*/
  free(l2);
  return l1;
}
//Zips n lists that were split using split_mod.
List ** zip_split(int n, List ** lists){
  List ** output = malloc(n*sizeof(List *));
  for (int i = 0; i < n; i++){
    output[i] = emptylist();
    //printf("Voor merge_empty\n");
    merge_empty(output[i],lists[i]);
    //printf("Na merge_empty\n");
  }
  Node * buffers[n];
  int i;
  for(i=0;i<n;i++)buffers[i] = lists[i]->head;
  int len = lists[0]->length;
  int out_list = 0;
  int count = 0;

  //printf("Voor loops\n");

  for (i = 0; i<len;i++){
    for (int j = 0; j<n;j++){
      if(buffers[j]!= NULL){
        Node * nnn = buffers[j];
        buffers[j] = nnn->next;
        add_node(nnn,output[out_list]);
        count++;
        //printf("Voor belachelijke if, i= %d, j = %d\n",i,j);
        if((out_list < n-1) && (count >= len && ((j < n-1 && buffers[j+1] != NULL && buffers[j+1]->data->sequence.l1num != nnn->data->sequence.l1num)
              || (j == n-1 && buffers[0] != NULL && buffers[0]->data->sequence.l1num != nnn->data->sequence.l1num)))){
        //if (count == len){
          out_list++;
          count = 0;
        }
        //printf("Na belachelijke if\n");
        //buffers[j] = buffers[j]->next;
      }
    }
  }
  //for (i = 0; i<n;i++)check_sequence_exact(output[i]);
  /*len = 0;
  for (i = 0; i<n;i++){
    Node * n = output[i]->head;
    len++;
    while(n->next != NULL && n->next->data != NULL){
      n = n->next;
      len++;
    }
  }
  printf("Totaal aantal chunks na zip_split: %d\n",len);*/
  //printf("Na loops\n");
  for(i=0;i<n;i++){
    free(lists[i]);
    //printf("Lengte lijst compress: %d\n", output[i]->length);
    //printf("I: %d, Eerste l1num: %d, laatste l1num: %d\n",i,output[i]->head->data->sequence.l1num, output[i]->tail->data->sequence.l1num);
    //printf("Eerste  element = %p, grootte uncompressed data = %d\n",output[i]->head, output[i]->head->data->uncompressed_data.n);
    //printf("Laatste element = %p, uncompressed data = %p, compressed data = %p\n",output[i]->tail, output[i]->tail->data->uncompressed_data.ptr, output[i]->tail->data->compressed_data);
  }

  free(lists);
  //free(buffers);
  //check_sequence(output);
  return output;
}

/*
* Takes 2 lists of chunks that are sorted by l1num and l2number
* Returns a new list that is sorted by l1num and l2num
* Cleans up the input lists
*/
List * sorted_merge(List * l1, List * l2){
  //check_sequence(l1);
  //check_sequence(l2);
  //printf("L1 en l2 ok\n");
    //printf("In sorted merge\n");
      //If one of the lists is empty, simply return the other list
      //Clean up the empty list
      if(l1 == NULL) return l2;
      else if(l1->head == NULL){
        free(l1);
        return l2;
      }
      else if (l2 == NULL) return l1;
      else if (l2->head == NULL){
        free(l2);
        return l1;
      }
      //Create output list and temporary variables
      List * output = emptylist();
      //printf("Voor merge_empty\n");
      merge_empty(output,l1);
      //printf("Merged l1\n");
      merge_empty(output,l2);
      //printf("Na merge_empty\n");
      //Iterator * iter1 = init_iterator(l1);
      //Iterator * iter2 = init_iterator(l2);
      Node * node1 = l1->head;//next_node(iter1);
      Node * node2 = l2->head;//next_node(iter2);
      chunk_t * chunk1 = node1->data;
      chunk_t * chunk2 = node2->data;
      Node * buf1 = node1->next;
      Node * buf2 = node2->next;

      //While both lists are not empty, 'zip' the lists by comparing elements
      int cont = 1;
      while(cont){
        if (output->tail == node1){
          //printf("Output->tail == node1\n");
          node1 = buf1;
          buf1 = buf1->next;//next_node(iter1);
          chunk1 = node1->data;
        }
        else if (output->tail == node2){
          //printf("Output->tail == node2\n");
          node2 = buf2;
          buf2 = buf2->next; //next_node(iter2);
          chunk2 = node2->data;
        }
        //printf("Chunk 1: %p, chunk2: %p\n",chunk1,chunk2);
        //Instead of moving chunks, we move entire nodes for efficiency.
        if(chunk1->sequence.l1num < chunk2->sequence.l1num
          || (chunk1->sequence.l1num == chunk2->sequence.l1num
          && chunk1->sequence.l2num < chunk2->sequence.l2num)) {
            //printf("Adding node1\n");
            add_node(node1,output);
            if(buf1 == NULL || buf1->data == NULL){
              cont = 0;
              add_node(node2,output);
            }
          }
        else {
          //printf("Adding node2\n");
          add_node(node2,output);
          //printf("Added node2\n");
          //printf("Buf2: %p\n",buf2);
          if(buf2 == NULL || buf2->data == NULL){
            //printf("Buf2 == NULL\n");
            cont = 0;
            add_node(node1,output);
          }
        }
        //printf("Added node\n");
      }
      //printf("Voorbij loop\n");
      //Add the entire remainder of the non-empty list at once.
      if (buf1 != NULL && buf1->data != NULL) add_nodes(buf1, l1->tail, output);
      else if (buf2 != NULL && buf2->data != NULL) add_nodes(buf2, l2->tail, output);

      //destroy_iterator(iter2);
      //destroy_iterator(iter1);
      free(l1);
      free(l2);
      //printf("Uit sorted merge\n");
      //check_sequence(output);
      //printf("Output sorted_merge ok\n");
      return output;
  }
