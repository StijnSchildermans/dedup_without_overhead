/*
 * Decoder for dedup files
 *
 * Copyright 2010 Princeton University.
 * All rights reserved.
 *
 * Originally written by Minlan Yu.
 * Largely rewritten by Christian Bienia.
 */

/*
 * The pipeline model for Encode is Fragment->FragmentRefine->Deduplicate->Compress->Reorder
 * Each stage has basically three steps:
 * 1. fetch a group of items from the queue
 * 2. process the items
 * 3. put them in the queue for the next stage
 */

#include <assert.h>
#include <strings.h>
#include <math.h>
#include <limits.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>
#include <unistd.h>
#include <string.h>

#include "util.h"
#include "dedupdef.h"
#include "encoder.h"
#include "debug.h"
#include "hashtable.h"
#include "config.h"
#include "rabin.h"
#include "mbuffer.h"
#include "chunk_list.h"
#include "iterator.h"
#include "thpool.h"


#ifdef ENABLE_PTHREADS
//#include "queue.h"
#include "binheap.h"
#include "tree.h"
#endif //ENABLE_PTHREADS

#ifdef ENABLE_GZIP_COMPRESSION
#include <zlib.h>
#endif //ENABLE_GZIP_COMPRESSION

#ifdef ENABLE_BZIP2_COMPRESSION
#include <bzlib.h>
#endif //ENABLE_BZIP2_COMPRESSION

#ifdef ENABLE_PTHREADS
#include <pthread.h>
#endif //ENABLE_PTHREADS

#ifdef ENABLE_PARSEC_HOOKS
#include <hooks.h>
#endif //ENABLE_PARSEC_HOOKS


#define INITIAL_SEARCH_TREE_SIZE 4096


//The configuration block defined in main
config_t * conf;

//Hash table data structure & utility functions
struct hashtable *cache;

static unsigned int hash_from_key_fn( void *k ) {
  //NOTE: sha1 sum is integer-aligned
  return ((unsigned int *)k)[0];
}

static int keys_equal_fn ( void *key1, void *key2 ) {
  return (memcmp(key1, key2, SHA1_LEN) == 0);
}


#ifdef ENABLE_STATISTICS
//Keep track of block granularity with 2^CHUNK_GRANULARITY_POW resolution (for statistics)
#define CHUNK_GRANULARITY_POW (7)
//Number of blocks to distinguish, CHUNK_MAX_NUM * 2^CHUNK_GRANULARITY_POW is biggest block being recognized (for statistics)
#define CHUNK_MAX_NUM (8*32)
//Map a chunk size to a statistics array slot
#define CHUNK_SIZE_TO_SLOT(s) ( ((s)>>(CHUNK_GRANULARITY_POW)) >= (CHUNK_MAX_NUM) ? (CHUNK_MAX_NUM)-1 : ((s)>>(CHUNK_GRANULARITY_POW)) )
//Get the average size of a chunk from a statistics array slot
#define SLOT_TO_CHUNK_SIZE(s) ( (s)*(1<<(CHUNK_GRANULARITY_POW)) + (1<<((CHUNK_GRANULARITY_POW)-1)) )


//Deduplication statistics (only used if ENABLE_STATISTICS is defined)
typedef struct {
  /* Cumulative sizes */
  size_t total_input; //Total size of input in bytes
  size_t total_dedup; //Total size of input without duplicate blocks (after global compression) in bytes
  size_t total_compressed; //Total size of input stream after local compression in bytes
  size_t total_output; //Total size of output in bytes (with overhead) in bytes

  /* Size distribution & other properties */
  unsigned int nChunks[CHUNK_MAX_NUM]; //Coarse-granular size distribution of data chunks
  unsigned int nDuplicates; //Total number of duplicate blocks
} stats_t;

//Arguments to pass to each thread
struct thread_args {
  //thread id, unique within a thread pool (i.e. unique for a pipeline stage)
  int tid;
  //number of queues available, first and last pipeline stage only
  int nqueues;
  //file descriptor, first pipeline stage only
  int fd;
  //List of chunks
  List * list;

  //char ** compressed_data;
  Compressed_data * compressed_data;

  List ** list_addr;
  //input file buffer, first pipeline stage & preloading only
  struct {
    void *buffer;
    size_t size;
  } input_file;

  stats_t * stats;
};

//Initialize a statistics record
static void init_stats(stats_t *s) {
  int i;

  assert(s!=NULL);
  s->total_input = 0;
  s->total_dedup = 0;
  s->total_compressed = 0;
  s->total_output = 0;

  for(i=0; i<CHUNK_MAX_NUM; i++) {
    s->nChunks[i] = 0;
  }
  s->nDuplicates = 0;
}

#ifdef ENABLE_PTHREADS

//Merge two statistics records: s1=s1+s2
static void merge_stats(stats_t *s1, stats_t *s2) {
  int i;

  assert(s1!=NULL);
  assert(s2!=NULL);
  s1->total_input += s2->total_input;
  s1->total_dedup += s2->total_dedup;
  s1->total_compressed += s2->total_compressed;
  s1->total_output += s2->total_output;

  for(i=0; i<CHUNK_MAX_NUM; i++) {
    s1->nChunks[i] += s2->nChunks[i];
  }
  s1->nDuplicates += s2->nDuplicates;
}
#endif //ENABLE_PTHREADS

//Print statistics
static void print_stats(stats_t *s) {
  const unsigned int unit_str_size = 7; //elements in unit_str array
  const char *unit_str[] = {"Bytes", "KB", "MB", "GB", "TB", "PB", "EB"};
  unsigned int unit_idx = 0;
  size_t unit_div = 1;

  assert(s!=NULL);

  //determine most suitable unit to use
  for(unit_idx=0; unit_idx<unit_str_size; unit_idx++) {
    unsigned int unit_div_next = unit_div * 1024;

    if(s->total_input / unit_div_next <= 0) break;
    if(s->total_dedup / unit_div_next <= 0) break;
    if(s->total_compressed / unit_div_next <= 0) break;
    if(s->total_output / unit_div_next <= 0) break;

    unit_div = unit_div_next;
  }

  printf("Total input size:              %14.2f %s\n", (float)(s->total_input)/(float)(unit_div), unit_str[unit_idx]);
  printf("Total output size:             %14.2f %s\n", (float)(s->total_output)/(float)(unit_div), unit_str[unit_idx]);
  printf("Effective compression factor:  %14.2fx\n", (float)(s->total_input)/(float)(s->total_output));
  printf("\n");

  //Total number of chunks
  unsigned int i;
  unsigned int nTotalChunks=0;
  for(i=0; i<CHUNK_MAX_NUM; i++) nTotalChunks+= s->nChunks[i];

  //Average size of chunks
  float mean_size = 0.0;
  for(i=0; i<CHUNK_MAX_NUM; i++) mean_size += (float)(SLOT_TO_CHUNK_SIZE(i)) * (float)(s->nChunks[i]);
  mean_size = mean_size / (float)nTotalChunks;

  //Variance of chunk size
  float var_size = 0.0;
  for(i=0; i<CHUNK_MAX_NUM; i++) var_size += (mean_size - (float)(SLOT_TO_CHUNK_SIZE(i))) *
                                             (mean_size - (float)(SLOT_TO_CHUNK_SIZE(i))) *
                                             (float)(s->nChunks[i]);

  printf("Total number of chunks: %d, Duplicate chunks: %d\n",nTotalChunks,s->nDuplicates);

  printf("Mean data chunk size:          %14.2f %s (stddev: %.2f %s)\n", mean_size / 1024.0, "KB", sqrtf(var_size) / 1024.0, "KB");
  printf("Amount of duplicate chunks:    %14.2f%%\n", 100.0*(float)(s->nDuplicates)/(float)(nTotalChunks));
  printf("Data size after deduplication: %14.2f %s (compression factor: %.2fx)\n", (float)(s->total_dedup)/(float)(unit_div), unit_str[unit_idx], (float)(s->total_input)/(float)(s->total_dedup));
  printf("Data size after compression:   %14.2f %s (compression factor: %.2fx)\n", (float)(s->total_compressed)/(float)(unit_div), unit_str[unit_idx], (float)(s->total_dedup)/(float)(s->total_compressed));
  printf("Output overhead:               %14.2f%%\n", 100.0*(float)(s->total_output-s->total_compressed)/(float)(s->total_output));
}

//variable with global statistics
stats_t stats;
#endif //ENABLE_STATISTICS


//Simple write utility function
/*static int write_file(int fd, u_char type, u_long len, u_char * content) {
  if (xwrite(fd, &type, sizeof(type)) < 0){
    perror("xwrite:");
    EXIT_TRACE("xwrite type fails\n");
    return -1;
  }
  if (xwrite(fd, &len, sizeof(len)) < 0){
    EXIT_TRACE("xwrite content fails\n");
  }
  if (xwrite(fd, content, len) < 0){
    EXIT_TRACE("xwrite content fails\n");
  }
  return 0;
}*/

/*
 * Helper function that creates and initializes the output file
 * Takes the file name to use as input and returns the file handle
 * The output file can be used to write chunks without any further steps
 */
static int create_output_file(char *outfile) {
  int fd;

  //Create output file
  fd = open(outfile, O_CREAT|O_TRUNC|O_WRONLY|O_TRUNC, S_IRGRP | S_IWUSR | S_IRUSR | S_IROTH);
  if (fd < 0) {
    EXIT_TRACE("Cannot open output file.");
  }

  //Write header
  if (write_header(fd, conf->compress_type)) {
    EXIT_TRACE("Cannot write output file header.\n");
  }
  return fd;
}

/*
 * Helper function that writes a chunk to an output file depending on
 * its state. The function will write the SHA1 sum if the chunk has
 * already been written before, or it will write the compressed data
 * of the chunk if it has not been written yet.
 *
 * This function will block if the compressed data is not available yet.
 * This function might update the state of the chunk if there are any changes.
 */
//NOTE: The serial version relies on the fact that chunks are processed in-order,
//      which means if it reaches the function it is guaranteed all data is ready.
/*static void write_chunk_to_file(int fd, chunk_t *chunk) {
  assert(chunk!=NULL);

  if(!chunk->header.isDuplicate) {
    //Unique chunk, data has not been written yet, do so now
    //printf("Unieke chunk\n");
    write_file(fd, TYPE_COMPRESS, chunk->compressed_data.n, chunk->compressed_data.ptr);
    //printf("Chunk is geschreven\n");
    mbuffer_free(&chunk->compressed_data);
    //printf("Buffer is vrij\n");
  } else {
    //printf("Dubbele chunk\n");
    //Duplicate chunk, data has been written to file before, just write SHA1
    write_file(fd, TYPE_FINGERPRINT, SHA1_LEN, (unsigned char *)(chunk->sha1));
    //printf("Dubbele chunk is gechreven\n");
  }
}*/
//#endif //ENABLE_PTHREADS

int rf_win;
int rf_win_dataprocess;

/*
 * Computational kernel of compression stage
 * Actions performed: Compress a data chunk
 */
void sub_Compress(chunk_t *chunk) {
    //size_t n;
    int r;

    //if(chunk->header.isDuplicate) return;

    assert(chunk!=NULL);
    switch (conf->compress_type) {
      case COMPRESS_NONE:
        //printf("No compression\n");
        //Simply duplicate the data
        //r = mbuffer_create(&chunk->compressed_data, n);
        /*if(r != 0) {
          EXIT_TRACE("Creation of compression buffer failed.\n");
        }*/
        //copy the block
        chunk->compressed_data.n = chunk->uncompressed_data.n;
        memcpy(chunk->compressed_data.ptr, chunk->uncompressed_data.ptr, chunk->uncompressed_data.n);
        break;
#ifdef ENABLE_GZIP_COMPRESSION
      case COMPRESS_GZIP:
        //printf("GZIP compression\n");
        //Gzip compression buffer must be at least 0.1% larger than source buffer plus 12 bytes
        //n = chunk->uncompressed_data.n + (chunk->uncompressed_data.n >> 9) + 12;
        //r = mbuffer_create(&chunk->compressed_data, n);
        /*if(r != 0) {
          EXIT_TRACE("Creation of compression buffer failed.\n");
        }*/
        //compress the block
        //printf("Voor compress functie\n");
        /*printf("Compressed data ptr = %p\n", chunk->compressed_data.ptr);
        printf("Compressed data n = %d\n", chunk->compressed_data.n);
        printf("Uncompressed data ptr = %p\n", chunk->uncompressed_data.ptr);
        printf("Uncompressed data n = %d\n", chunk->uncompressed_data.n);*/

        r = compress(chunk->compressed_data.ptr, &chunk->compressed_data.n, chunk->uncompressed_data.ptr, chunk->uncompressed_data.n);
        //printf("Na compress functie\n");
        //memcpy(chunk->compressed_data.ptr, chunk->uncompressed_data.ptr, chunk->uncompressed_data.n);
        //chunk->compressed_data.n = chunk->uncompressed_data.n;
        if (r != Z_OK) {
          EXIT_TRACE("Compression failed. Error code: %d\n",r);
        }
        //Shrink buffer to actual size
        /*if(n < chunk->compressed_data.n) {
          r = mbuffer_realloc(&chunk->compressed_data, n);
          assert(r == 0);
        }*/
        break;
#endif //ENABLE_GZIP_COMPRESSION
#ifdef ENABLE_BZIP2_COMPRESSION
      case COMPRESS_BZIP2:
        //printf("BZIP2 compression");
        //Bzip compression buffer must be at least 1% larger than source buffer plus 600 bytes
        n = chunk->uncompressed_data.n + (chunk->uncompressed_data.n >> 6) + 600;
        r = mbuffer_create(&chunk->compressed_data, n);
        if(r != 0) {
          EXIT_TRACE("Creation of compression buffer failed.\n");
        }
        //compress the block
        unsigned int int_n = n;
        r = BZ2_bzBuffToBuffCompress(chunk->compressed_data.ptr, &int_n, chunk->uncompressed_data.ptr, chunk->uncompressed_data.n, 9, 0, 30);
        n = int_n;
        if (r != BZ_OK) {
          EXIT_TRACE("Compression failed\n");
        }
        //Shrink buffer to actual size
        if(n < chunk->compressed_data.n) {
          r = mbuffer_realloc(&chunk->compressed_data, n);
          assert(r == 0);
        }
        break;
    #endif //ENABLE_BZIP2_COMPRESSION
      default:
        EXIT_TRACE("Compression type not implemented.\n");
        break;
    }
    //printf("Voor mbuffer_free\n");
    mbuffer_free(&chunk->uncompressed_data);
    //printf("Na mbuffer_free\n");

#ifdef ENABLE_PTHREADS
    chunk->header.state = CHUNK_STATE_COMPRESSED;
#endif //ENABLE_PTHREADS

}

/*
 * Pipeline stage function of compression stage
 *
 * Actions performed:
 *  - Dequeue items from compression queue
 *  - Execute compression kernel for each item
 *  - Enqueue each item into send queue
 */
//#ifdef ENABLE_PTHREADS
void Compress(void * targs) {

  struct thread_args *args = (struct thread_args *)targs;
  List * list = args->list;

  #ifdef ENABLE_STATISTICS
  /*stats_t *thread_stats = malloc(sizeof(stats_t));
  if(thread_stats == NULL) EXIT_TRACE("Memory allocation failed.\n");*/
  stats_t * thread_stats = args->stats;
  init_stats(thread_stats);
  #endif //ENABLE_STATISTICS

//printf("In compress\n");

//Allocate memory for compressed data buffers

int total_chunks = 0;
int duplicate_chunks = 0;
size_t total_size = 0;
void * mbuffers;
//size_t total_compressed_size = 0;
chunk_t * chunk_refs[1000];

int write_buffers_index = 0;

//int total = 0;


Iterator * iter = init_iterator(list);

  //printf("Lengte lijst: %d\n",list->length);


while(hasNext(iter)){
  //printf("Blijkbaar is er nog een element..\n");
  //Take next chunk
  chunk_t * c = next(iter);
  //printf("L1: %d, L2: %d\n",c->sequence.l1num, c->sequence.l2num);

  chunk_refs[total_chunks] = c;
  total_chunks++;
  //printf("Total chunks: %d\n", total_chunks);



  // If chunk is unique, update counter to reserve memory for compressed buffer.
  if(c->header.isDuplicate) duplicate_chunks++;
  else{
    thread_stats->total_dedup += c->uncompressed_data.n;
    size_t * size = &c->compressed_data.n;
    //if (&c->uncompressed_data == NULL) printf("Dikke zever\n");
    if(conf->compress_type == COMPRESS_NONE) *size =  c->uncompressed_data.n;
    else *size = c->uncompressed_data.n + (c->uncompressed_data.n >> 9) + 12;
    //printf("Size compressed data = %d\n", *size);
    total_size += *size;
  }


  //If we found 1000 chunks or found the last chunk, process the batch.
  if(total_chunks == 1000 || !hasNext(iter)){
    int index = 0;
    //printf("Total size: %d\n", total_size);
    mbuffers = malloc(total_size);
    total_size = 0;
    //printf("Total chunks: %d\n",total_chunks);
    //total = 0;
    for(int i = 0; i<total_chunks; i++){

      //printf("In compressie\n");
      chunk_t * chunk = chunk_refs[i];
      if(!chunk->header.isDuplicate){
        chunk->compressed_data.ptr = mbuffers + index;
        index += chunk->compressed_data.n;


        /*char * ptr = &mbuffers[index] + sizeof(mcb_t);
        mbuffer_t * m = &(chunk->compressed_data);
        m->mcb = (mcb_t*) (&mbuffers[index]);
        m->ptr = ptr;
        m->mcb->i = 1;
        m->mcb->ptr = ptr;
        #ifdef ENABLE_MBUFFER_CHECK
          m->check_flag = MBUFFER_CHECK_MAGIC;
        #endif
        index += sizeof(mcb_t) + m->n;*/
        //printf("Compressed data vlak voor sub_compress: %d\n",m->n);

        //printf("Pointer naar mbuffer compressed data: %p\n",ptr);
        //printf("Unique chunks = %d, total = %d\n",unique_chunks,total);
        //printf("Voor sub_compress\n");
        //printf("Waarde i: %d\n",i);
        //printf("Totale waarde vrijgemaakt voor compressed data: %d, waarde ingenomen: %d, grootte volgende: %d\n", total_malloced, index, m->n);
        sub_Compress(chunk);
        //printf("Achter sub_compress\n");
        thread_stats->total_compressed += chunk->compressed_data.n;
        total_size += chunk->compressed_data.n;
      }
      //printf("Aan einde compressie\n");
      //mbuffer_free(&chunk->uncompressed_data);
    }
    //printf("Voorbij compressie\n");
    int write_buffer_size = duplicate_chunks * SHA1_LEN + total_chunks * 9 + total_size;
    char * write_buffer = malloc(write_buffer_size);

    index = 0;
    for(int i = 0; i < total_chunks; i++){
      //printf("In schrijven\n");
      chunk_t * chunk = chunk_refs[i];
      if (chunk->header.isDuplicate){
        thread_stats->nDuplicates++;
        write_buffer[index] = 0;
        *((u_long*) (&write_buffer[index+1])) = SHA1_LEN;
        index += 9;
        memcpy(write_buffer + index, &chunk->sha1, SHA1_LEN);
        index += SHA1_LEN;
        //mbuffer_free(&chunk->uncompressed_data);
        //free(chunk);
      }
      else{
        write_buffer[index] = 1;
        *((u_long*) (&write_buffer[index+1])) = chunk->compressed_data.n;
        index += 9;
        memcpy(write_buffer + index, chunk->compressed_data.ptr, chunk->compressed_data.n);
        //if(chunk->compressed_data.n>100000) printf("Chunk van meer dan 100 kB\n");
        index += chunk->compressed_data.n;
        //mbuffer_free(&chunk->compressd_data);
        //free(chunk);
      }
      //free(chunk);
    }
    //free(mbuffers);
    //write_buffers[write_buffers_index] = write_buffer;
    args->compressed_data[write_buffers_index].data = write_buffer;
    args->compressed_data[write_buffers_index].size = write_buffer_size;
    write_buffers_index++;
    total_chunks = 0;
    duplicate_chunks = 0;
    //printf("Voorbij schrijven\n");


  }

   //total++;
   //printf("Processed a total of %d chunks\n",total);
}
  //printf("Voorbij loop\n");
  destroy_iterator(iter);
  //printf("Iterator kaput\n");
  //destroy_soft(list);
  //destroy(list);
  //printf("List kaput\n");
  free(list);

  //printf("Uit compress\n");


  /*#ifdef ENABLE_STATISTICS
    return thread_stats;
  #else
    return NULL;
  #endif //ENABLE_STATISTICS*/

}
//#endif //ENABLE_PTHREADS

 /* Computational kernel of deduplication stage
 *
 * Actions performed:
 *  - Calculate SHA1 signature for each incoming data chunk
 *  - Perform database lookup to determine chunk redundancy status
 *  - On miss add chunk to database
 *  - Returns chunk redundancy status */
int sub_Deduplicate(chunk_t *chunk) {
  int isDuplicate;
  int isFirst = 1;
  chunk_t *entry;

  assert(chunk!=NULL);
  assert(chunk->uncompressed_data.ptr!=NULL);

  SHA1_Digest(chunk->uncompressed_data.ptr, chunk->uncompressed_data.n, (unsigned char *)(chunk->sha1));

  //Query database to determine whether we've seen the data chunk before
#ifdef ENABLE_PTHREADS
  pthread_mutex_t *ht_lock = hashtable_getlock(cache, (void *)(chunk->sha1));
  pthread_mutex_lock(ht_lock);
#endif
  entry = (chunk_t *)hashtable_search(cache, (void *)(chunk->sha1));
  isDuplicate = (entry != NULL);
  if (isDuplicate){
    if (entry->sequence.l1num > chunk->sequence.l1num
      || (entry->sequence.l1num == chunk->sequence.l1num
        && entry->sequence.l2num > chunk->sequence.l2num)){
      isFirst = 1;
      entry->header.isDuplicate = 1;
      chunk->header.isDuplicate = 0;
      entry->compressed_data_ref = chunk;
      mbuffer_free(&entry->uncompressed_data);
      //free(&entry->uncompressed_data);
      if (hashtable_insert(cache, (void *)(chunk->sha1), (void *)chunk) == 0) {
        EXIT_TRACE("hashtable_insert failed");
      }
    }
    else{
      isFirst = 0;
      chunk->header.isDuplicate = 1;
      entry->header.isDuplicate = 0;
      chunk->compressed_data_ref = entry;
      mbuffer_free(&chunk->uncompressed_data);
      //free(&chunk->uncompressed_data);
    }
  }
  else{
    chunk->header.isDuplicate = 0;
    // Cache miss: Create entry in hash table and forward data to compression stage
    #ifdef ENABLE_PTHREADS
      pthread_mutex_init(&chunk->header.lock, NULL);
      pthread_cond_init(&chunk->header.update, NULL);
    #endif
      //NOTE: chunk->compressed_data.buffer will be computed in compression stage
    if (hashtable_insert(cache, (void *)(chunk->sha1), (void *)chunk) == 0) {
        EXIT_TRACE("hashtable_insert failed");
    }
  }
#ifdef ENABLE_PTHREADS
  pthread_mutex_unlock(ht_lock);
#endif

  return (isDuplicate && !isFirst);
}

/* Pipeline stage function of deduplication stage
 *
 * Actions performed:
 *  - Take input data from fragmentation stages
 *  - Execute deduplication kernel for each data chunk
 *  - Route resulting package either to compression stage or to reorder stage, depending on deduplication status */
#ifdef ENABLE_PTHREADS
void Deduplicate(void * targs) {
  struct thread_args *args = (struct thread_args *)targs;
  List * list = args->list;
  Node * node;
  Node * buffer = list->head;

  #ifdef ENABLE_STATISTICS
    /*stats_t *thread_stats = malloc(sizeof(stats_t));
    if(thread_stats == NULL) {
      EXIT_TRACE("Memory allocation failed.\n");
    }
    */
    stats_t* thread_stats = args->stats;
    init_stats(thread_stats);
  #endif //ENABLE_STATISTICS

  //Iterator * iter = init_iterator(list);
  int len = list->length;
  for(int i = 0; i< len; i++) {
    //chunk = next(iter);
    node = buffer;
    buffer = buffer->next;
    //if(node->data == NULL) printf("Data is null, i = %d\n",i);

    assert(node->data!=NULL);

    //Do the processing
    sub_Deduplicate(node->data);

/*#ifdef ENABLE_STATISTICS
    if(isDuplicate) {
      thread_stats->nDuplicates++;
    } else {
      thread_stats->total_dedup += node->data->uncompressed_data.n;
    }
#endif //ENABLE_STATISTICS*/
    //Enqueue chunk either into compression queue or into send queue
    //if(isDuplicate) thread_stats->nDuplicates++;
    //list_index++;
  }
  //destroy_soft(list);
  //destroy_iterator(iter);

/*#ifdef ENABLE_STATISTICS
  return thread_stats;
#else
  return NULL;
#endif //ENABLE_STATISTICS*/
}
#endif //ENABLE_PTHREADS

 /* Pipeline stage function and computational kernel of refinement stage
 *
 * Actions performed:
 *  - Take coarse chunks from fragmentation stage
 *  - Partition data block into smaller chunks with Rabin rolling fingerprints
 *  - Send resulting data chunks to deduplication stage
 *
 * Notes:
 *  - Allocates mbuffers for fine-granular chunks*/
void FragmentRefine(void * targs) {
  //printf("In fragmentRefine\n");
  struct thread_args *args = (struct thread_args *)targs;
  int r;
  List * list = (List *)args->list;

  chunk_t *temp;
  chunk_t *chunk;
  u32int * rabintab = malloc(256*sizeof rabintab[0]);
  u32int * rabinwintab = malloc(256*sizeof rabintab[0]);
  if(rabintab == NULL || rabinwintab == NULL) EXIT_TRACE("Memory allocation failed.\n");

#ifdef ENABLE_STATISTICS
  /*stats_t *thread_stats = malloc(sizeof(stats_t));
  if(thread_stats == NULL) {
    EXIT_TRACE("Memory allocation failed.\n");
  }
  */

  stats_t *thread_stats = args->stats;
  init_stats(thread_stats);
#endif //ENABLE_STATISTICS

  int chcount = 0;
  List * refined = emptylist();
  Iterator * iter = init_iterator(list);
  while (hasNext(iter)) {
    //printf("Aan begin loop\n");
    //printf("Voor keuze nieuwe chunk\n");
    chunk = next(iter);
    //printf("Nieuwe chunk = %p\n", chunk);
    assert(chunk!=NULL);

    rabininit(rf_win, rabintab, rabinwintab);

    int split;
    chcount = 0;
    do {
      //printf("Aan begin loop\n");
      //Find next anchor with Rabin fingerprint
      int offset = rabinseg(chunk->uncompressed_data.ptr, chunk->uncompressed_data.n, rf_win, rabintab, rabinwintab);
      //Can we split the buffer?
      if(offset < chunk->uncompressed_data.n) {
        //Allocate a new chunk and create a new memory buffer
        temp = (chunk_t *)malloc(sizeof(chunk_t));
        if(temp==NULL) EXIT_TRACE("Memory allocation failed.\n");
        temp->header.state = chunk->header.state;
        temp->sequence.l1num = chunk->sequence.l1num;

        //split it into two pieces
        r = mbuffer_split(&chunk->uncompressed_data, &temp->uncompressed_data, offset);
        if(r!=0) EXIT_TRACE("Unable to split memory buffer in refinement stage.\n");

        //Set correct state and sequence numbers
        chunk->sequence.l2num = chcount;
        chunk->isLastL2Chunk = FALSE;
        chcount++;

        #ifdef ENABLE_STATISTICS
          //update statistics
          thread_stats->nChunks[CHUNK_SIZE_TO_SLOT(chunk->uncompressed_data.n)]++;
        #endif //ENABLE_STATISTICS

        //put it into send buffer
        add(chunk, refined);
        //prepare for next iteration
        chunk = temp;
        split = 1;
      } else {
        //End of buffer reached, don't split but simply enqueue it
        //Set correct state and sequence numbers
        chunk->sequence.l2num = chcount;
        chunk->isLastL2Chunk = TRUE;

        #ifdef ENABLE_STATISTICS
          //update statistics
          thread_stats->nChunks[CHUNK_SIZE_TO_SLOT(chunk->uncompressed_data.n)]++;
        #endif //ENABLE_STATISTICS

        add(chunk, refined);
        //prepare for next iteration
        chunk = NULL;
        split = 0;
      }
      //printf("Aan einde loop\n");
    } while(split);
    //printf("Aan einde loop\n");
  }
  //printf("Voorbij loop\n");
  //printf("Lengte refined: %d\n", refined->length);
  *(args->list_addr) = refined;
  //printf("Lijst adressen aangepast\n");
  //destroy_soft(list);
  //printf("Lijst vrijgegeven\n");

  free(rabintab);
  free(rabinwintab);
  destroy_iterator(iter);
/*
#ifdef ENABLE_STATISTICS
  return thread_stats;
#else
  return NULL;
#endif //ENABLE_STATISTICS
//printf("Uit fragmentRefine\n");*/
}

/*
 * Pipeline stage function of fragmentation stage
 *
 * Actions performed:
 *  - Read data from file (or preloading buffer)
 *  - Perform coarse-grained chunking
 *  - Send coarse chunks to refinement stages for further processing
 *
 * Notes:
 * This pipeline stage is a bottleneck because it is inherently serial. We
 * therefore perform only coarse chunking and pass on the data block as fast
 * as possible so that there are no delays that might decrease scalability.
 * With very large numbers of threads this stage will not be able to keep up
 * which will eventually limit scalability. A solution to this is to increase
 * the size of coarse-grained chunks with a comparable increase in total
 * input size.
 */
#ifdef ENABLE_PTHREADS
List * Fragment(void * targs){
  //printf("In fragment\n");
  struct thread_args *args = (struct thread_args *)targs;
  size_t preloading_buffer_seek = 0;
  int fd = args->fd;
  int r;
  sequence_number_t anchorcount = 0;

  List * list = emptylist();

  chunk_t *temp = NULL;
  chunk_t *chunk = NULL;
  u32int * rabintab = malloc(256*sizeof rabintab[0]);
  u32int * rabinwintab = malloc(256*sizeof rabintab[0]);
  if(rabintab == NULL || rabinwintab == NULL) {
    EXIT_TRACE("Memory allocation failed.\n");
  }

  rf_win_dataprocess = 0;
  rabininit(rf_win_dataprocess, rabintab, rabinwintab);

  //Sanity check
  if(MAXBUF < 8 * ANCHOR_JUMP) {
    printf("WARNING: I/O buffer size is very small. Performance degraded.\n");
    fflush(NULL);
  }

  //read from input file / buffer
  while (1) {
    //printf("Aan het begin loop\n");
    size_t bytes_left; //amount of data left over in last_mbuffer from previous iteration

    //Check how much data left over from previous iteration resp. create an initial chunk
    //printf("Voor nakijken temp\n");
    if(temp != NULL) {
      bytes_left = temp->uncompressed_data.n;
    } else {
      bytes_left = 0;
    }
    //printf("Na nakijken temp\n");

    //Make sure that system supports new buffer size
    if(MAXBUF+bytes_left > SSIZE_MAX) {
      EXIT_TRACE("Input buffer size exceeds system maximum.\n");
    }
    //Allocate a new chunk and create a new memory buffer
    chunk = (chunk_t*)malloc(sizeof(chunk_t));
    if(chunk==NULL) EXIT_TRACE("Memory allocation failed.\n");
    mbuffer_create(&chunk->uncompressed_data, MAXBUF+bytes_left);

    if(bytes_left > 0) {
      //printf("In bytes left > 0\n");
      //printf("Chunk = %p\n",chunk);
      //printf("Bytes lef")
      //FIXME: Short-circuit this if no more data available

      //"Extension" of existing buffer, copy sequence number and left over data to beginning of new buffer
      chunk->header.state = CHUNK_STATE_UNCOMPRESSED;
      chunk->sequence.l1num = temp->sequence.l1num;

      //printf("Voor memcpy\n");
      //NOTE: We cannot safely extend the current memory region because it has already been given to another thread
      memcpy(chunk->uncompressed_data.ptr, temp->uncompressed_data.ptr, temp->uncompressed_data.n);
      mbuffer_free(&temp->uncompressed_data);
      //printf("Achter mbuffer_free\n");
      //printf("Temp: %p\n",temp);
      free(temp);
      //printf("Freed temp\n");
      temp = NULL;
      //printf("Uit bytes left > 0\n");
    } else {
      //brand new mbuffer, increment sequence number
      chunk->header.state = CHUNK_STATE_UNCOMPRESSED;
      chunk->sequence.l1num = anchorcount;
      anchorcount++;
    }
    //printf("Voor inlezen\n");
    //Read data until buffer full
    size_t bytes_read=0;
    if(conf->preloading) {
      //printf("Aan het inlezen\n");
      size_t max_read = MIN(MAXBUF, args->input_file.size-preloading_buffer_seek);
      memcpy(chunk->uncompressed_data.ptr+bytes_left, args->input_file.buffer+preloading_buffer_seek, max_read);
      bytes_read = max_read;
      preloading_buffer_seek += max_read;
      //printf("Uit inlezen\n");
    } else {
      //printf("Aan het inlezen\n");
      while(bytes_read < MAXBUF) {
        int r = read(fd, chunk->uncompressed_data.ptr+bytes_left+bytes_read, MAXBUF-bytes_read);
        if(r<0) switch(errno) {
          case EAGAIN:
            EXIT_TRACE("I/O error: No data available\n");break;
          case EBADF:
            EXIT_TRACE("I/O error: Invalid file descriptor\n");break;
          case EFAULT:
            EXIT_TRACE("I/O error: Buffer out of range\n");break;
          case EINTR:
            EXIT_TRACE("I/O error: Interruption\n");break;
          case EINVAL:
            EXIT_TRACE("I/O error: Unable to read from file descriptor\n");break;
          case EIO:
            EXIT_TRACE("I/O error: Generic I/O error\n");break;
          case EISDIR:
            EXIT_TRACE("I/O error: Cannot read from a directory\n");break;
          default:
            EXIT_TRACE("I/O error: Unrecognized error\n");break;
        }
        if(r==0) break;
        bytes_read += r;
      }
      //printf("Uit inlezen\n");
    }
    //No data left over from last iteration and also nothing new read in, simply clean up and quit
    if(bytes_left + bytes_read == 0) {
      mbuffer_free(&chunk->uncompressed_data);
      #ifdef ENABLE_MBUFFER_CHECK
          m->check_flag=0;
      #endif
      free(chunk);
      //printf("Freed chunk: %p\n",chunk);
      chunk = NULL;
      break;
    }
    //Shrink buffer to actual size
    if(bytes_left+bytes_read < chunk->uncompressed_data.n) {
      r = mbuffer_realloc(&chunk->uncompressed_data, bytes_left+bytes_read);
      assert(r == 0);
    }
    //Check whether any new data was read in, enqueue last chunk if not
    if(bytes_read == 0) {
      //printf("Voor add\n");
      add(chunk, list);
      //printf("Achter add\n");
      //NOTE: No need to empty a full send_buf, we will break now and pass everything on to the queue
      break;
    }
    //partition input block into large, coarse-granular chunks
    int split;
    do {
      //printf("Aan het begin split loop\n");
      split = 0;
      //Try to split the buffer at least ANCHOR_JUMP bytes away from its beginning
      if(ANCHOR_JUMP < chunk->uncompressed_data.n) {
        int offset = rabinseg(chunk->uncompressed_data.ptr + ANCHOR_JUMP, chunk->uncompressed_data.n - ANCHOR_JUMP, rf_win_dataprocess, rabintab, rabinwintab);
        //Did we find a split location?
        if(offset == 0) {
          //Split found at the very beginning of the buffer (should never happen due to technical limitations)
          assert(0);
          split = 0;
        } else if(offset + ANCHOR_JUMP < chunk->uncompressed_data.n) {
          //Split found somewhere in the middle of the buffer
          //Allocate a new chunk and create a new memory buffer
          temp = (chunk_t *)malloc(sizeof(chunk_t));
          if(temp==NULL) EXIT_TRACE("Memory allocation failed.\n");

          int size = offset + ANCHOR_JUMP;

          mbuffer_create(&temp->uncompressed_data,size);
          memcpy(temp->uncompressed_data.ptr, chunk->uncompressed_data.ptr ,size);

          //split it into two pieces
          //r = mbuffer_split(&chunk->uncompressed_data, &temp->uncompressed_data, offset + ANCHOR_JUMP);

          void * p = chunk->uncompressed_data.ptr + size;
          int p_n = chunk->uncompressed_data.n - size;
          mcb_t * p_mcb = chunk->uncompressed_data.mcb;

          chunk->uncompressed_data.ptr = temp->uncompressed_data.ptr;
          chunk->uncompressed_data.n = size;
          chunk->uncompressed_data.mcb = temp->uncompressed_data.mcb;

          temp->uncompressed_data.ptr = p;
          temp->uncompressed_data.n = p_n;
          temp->uncompressed_data.mcb = p_mcb;

          #ifdef ENABLE_MBUFFER_CHECK
            m2->check_flag=MBUFFER_CHECK_MAGIC;
          #endif

          //if(r!=0) EXIT_TRACE("Unable to split memory buffer in fragmentation stage.\n");
          temp->header.state = CHUNK_STATE_UNCOMPRESSED;
          temp->sequence.l1num = anchorcount;
          anchorcount++;

          //put it into send buffer
          //printf("Voor add \n");
          add(chunk, list);
          //printf("Achter add\n");
          //prepare for next iteration
          chunk = temp;
          temp = NULL;
          split = 1;
        } else {
          //Due to technical limitations we can't distinguish the cases "no split" and "split at end of buffer"
          //This will result in some unnecessary (and unlikely) work but yields the correct result eventually.
          temp = chunk;
          chunk = NULL;
          split = 0;
        }
      } else {
        //NOTE: We don't process the stub, instead we try to read in more data so we might be able to find a proper split.
        //      Only once the end of the file is reached do we get a genuine stub which will be enqueued right after the read operation.
        temp = chunk;
        chunk = NULL;
        split = 0;
      }
      //printf("Aan het einde split loop\n");
    } while(split);
    //printf("Aan het einde loop\n");
  }
  //printf("Uit loop");
  free(rabintab);
  free(rabinwintab);
  //printf("Uit fragment\n");

  return list;
}
#endif //ENABLE_PTHREADS

/*
 * Pipeline stage function of reorder stage
 *
 * Actions performed:
 *  - Receive chunks from compression and deduplication stage
 *  - Check sequence number of each chunk to determine correct order
 *  - Cache chunks that arrive out-of-order until predecessors are available
 *  - Write chunks in-order to file (or preloading buffer)
 *
 * Notes:
 *  - This function blocks if the compression stage has not finished supplying
 *    the compressed data for a duplicate chunk.
 */
#ifdef ENABLE_PTHREADS
void Write(Compressed_data ** data, int * counts) {


  int fd = 0;
  //chunk_t *chunk;
  fd = create_output_file(conf->outfile);
  /*
  Iterator * iter = init_iterator(list);
  int size = 0;
  //int chunks = 0;
  while(hasNext(iter)) {
    //chunks++;
    //write as many chunks as possible, current chunk is next in sequence
    chunk = next(iter);
    //printf("Chunk: %p\n",chunk);
    //write_chunk_to_file(fd, chunk);
    if(chunk->header.isDuplicate) size += (9 + SHA1_LEN);
    else size += (9 + chunk->compressed_data.n);
  }
  //printf("Size = %d\n",size);
  //printf("Chunks: %d\n",chunks);
  //int chunks = 0;
  //int duplicate = 0;
  unsigned char * toWrite = (unsigned char *) malloc(size);
  int index = 0;
  //int type_fingerprint = TYPE_FINGERPRINT;
  //int type_compress = TYPE_COMPRESS;
  //u_long sha1_len = SHA1_LEN;
  //printf("Grootte type: %d\n",sizeof(TYPE_FINGERPRINT));
  //printf("Grootte int: %d\n",sizeof(type_fingerprint));
  reset(iter);
  while(hasNext(iter)){
    //chunks++;
    chunk = next(iter);
    if(chunk->header.isDuplicate){
      //duplicate++;
      //for(int i = 0; i<sizeof(TYPE_FINGERPRINT);i++) toWrite+index+i = 0;
      toWrite[index] = 0;
      index++;
      //memcpy(toWrite + index, &sha1_len, 8);
      *((u_long*) (&toWrite[index])) = SHA1_LEN;
      index += 8;
      memcpy(toWrite + index, &chunk->sha1, SHA1_LEN);
      index += SHA1_LEN;
    }
    else {
      toWrite[index] = 1;
      index ++;
      //memcpy(toWrite + index, &chunk->compressed_data.n, sizeof(size_t));
      *((u_long*) (&toWrite[index])) = chunk->compressed_data.n;
      index += 8;
      memcpy(toWrite + index, chunk->compressed_data.ptr, chunk->compressed_data.n);
      //if(chunk->compressed_data.n>100000) printf("Chunk van meer dan 100 kB\n");
      index += chunk->compressed_data.n;
      //mbuffer_free(&chunk->compressed_data);
    }
      mbuffer_free(&chunk->uncompressed_data);
      free(chunk);
    //mbuffer_destroy
    //free(chunk);
  }*/
  //printf("Uiteindelijke waarde index: %d\n",index);
  //printf("Chunks en duplicates volgens write: %d, %d\n",chunks,duplicate);
  //printf("Chunks: %d\n",chunks);

  //for (int i = 0; i < conf->nthreads; i++) printf("Counts[%d] = %d\n",i,counts[i]);

  for (int i = 0; i<conf->nthreads; i++){
    for (int j = 0; j < counts[i];j++){
      //printf("Writing %d %d\n",i,j);
      xwrite(fd, data[i][j].data, data[i][j].size);
      free(data[i][j].data);
      //free(&(data[i][j]);
    }
    free(data[i]);
  }
  //destroy(list);
  //destroy_iterator(iter);
  close(fd);
  //return NULL;*/
}
#endif //ENABLE_PTHREADS



//Allocate all mbuffers at once to reduce malloc calls from threads
/*void assign_compression_buffers(int duplicates, int size, List ** lists){
  char * mbuffers = (char *)malloc(duplicates * sizeof(mcb_t) + size);

  int index = 0;
  for (int i = 0; i < conf->nthreads; i++){
    Iterator * iter = init_iterator(lists[i]);
    while(hasNext(iter)){
      chunk_t * c = next(iter);
      if(!c->header.isDuplicate){
        int buffer_size;
        if (conf->compress_type == COMPRESS_NONE) buffer_size = c->uncompressed_data.n;
        else buffer_size = c->uncompressed_data.n + (c->uncompressed_data.n >> 9) + 12;

        //void * ptr = malloc(buffer_size);
        char * ptr = &mbuffers[index] + sizeof(mcb_t);
        mbuffer_t * m = &c->compressed_data;
        m->mcb = (mcb_t*) (&mbuffers[index]);
        m->ptr = ptr;
        m->n = buffer_size;
        m->mcb->i = 1;
        m->mcb->ptr = ptr;
        #ifdef ENABLE_MBUFFER_CHECK
          m->check_flag=MBUFFER_CHECK_MAGIC;
        #endif
        index += sizeof(mcb_t) + buffer_size;
      }
    }
    destroy_iterator(iter);
  }
}*/

/*--------------------------------------------------------------------------*/
/* Encode
 * Compress an input stream
 *
 * Arguments:
 *   conf:    Configuration parameters
 *
 */
void Encode(config_t * _conf) {
  struct stat filestat;
  int32 fd;

  conf = _conf;

#ifdef ENABLE_STATISTICS
  init_stats(&stats);
#endif

  //Create chunk cache
  cache = hashtable_create(65536, hash_from_key_fn, keys_equal_fn, FALSE);
  if(cache == NULL) {
    printf("ERROR: Out of memory\n");
    exit(1);
  }

#ifdef ENABLE_PTHREADS
  struct thread_args data_process_args;
#else
  struct thread_args generic_args;
#endif //ENABLE_PTHREADS

  /* src file stat */
  if (stat(conf->infile, &filestat) < 0)
      EXIT_TRACE("stat() %s failed: %s\n", conf->infile, strerror(errno));

  if (!S_ISREG(filestat.st_mode))
    EXIT_TRACE("not a normal file: %s\n", conf->infile);
#ifdef ENABLE_STATISTICS
  stats.total_input = filestat.st_size;
#endif //ENABLE_STATISTICS

  /* src file open */
  if((fd = open(conf->infile, O_RDONLY | O_LARGEFILE)) < 0)
    EXIT_TRACE("%s file open error %s\n", conf->infile, strerror(errno));

  //Load entire file into memory if requested by user
  void *preloading_buffer = NULL;
  if(conf->preloading) {
    size_t bytes_read=0;
    int r;

    preloading_buffer = malloc(filestat.st_size);
    if(preloading_buffer == NULL)
      EXIT_TRACE("Error allocating memory for input buffer.\n");

    //Read data until buffer full
    while(bytes_read < filestat.st_size) {
      r = read(fd, preloading_buffer+bytes_read, filestat.st_size-bytes_read);
      if(r<0) switch(errno) {
        case EAGAIN:
          EXIT_TRACE("I/O error: No data available\n");break;
        case EBADF:
          EXIT_TRACE("I/O error: Invalid file descriptor\n");break;
        case EFAULT:
          EXIT_TRACE("I/O error: Buffer out of range\n");break;
        case EINTR:
          EXIT_TRACE("I/O error: Interruption\n");break;
        case EINVAL:
          EXIT_TRACE("I/O error: Unable to read from file descriptor\n");break;
        case EIO:
          EXIT_TRACE("I/O error: Generic I/O error\n");break;
        case EISDIR:
          EXIT_TRACE("I/O error: Cannot read from a directory\n");break;
        default:
          EXIT_TRACE("I/O error: Unrecognized error\n");break;
      }
      if(r==0) break;
      bytes_read += r;
    }
#ifdef ENABLE_PTHREADS
    data_process_args.input_file.size = filestat.st_size;
    data_process_args.input_file.buffer = preloading_buffer;
#else
    generic_args.input_file.size = filestat.st_size;
    generic_args.input_file.buffer = preloading_buffer;
#endif //ENABLE_PTHREADS
  }

  /* Variables for 3 thread pools and 2 pipeline stage threads.
   * The first and the last stage are serial (mostly I/O).
   */
  /*pthread_t threads_anchor[MAX_THREADS],
    threads_chunk[MAX_THREADS],
    threads_compress[MAX_THREADS];*/

  data_process_args.tid = 0;
  data_process_args.fd = fd;

#ifdef ENABLE_PARSEC_HOOKS
    __parsec_roi_begin();
#endif


  threadpool pool = thpool_init(14);

  printf("Entering fragmentation stage\n");
  //fflush(stdout);
  int threadCount = conf->nthreads;
  List * fragmented = Fragment(&data_process_args);
  //clean up after preloading
    if(conf->preloading) free(preloading_buffer);

  printf("Entering refinement stage\n");
  List ** refined = split(threadCount, fragmented);
  //printf("Successfully split list\n");
  int i = 0;
  struct thread_args anchor_thread_args[conf->nthreads];
  stats_t threads_anchor_rv[conf->nthreads];
  for (i = 0; i < conf->nthreads; i ++) {
     anchor_thread_args[i].tid = i;
     anchor_thread_args[i].list = refined[i];
     //printf("Length refined[i]: %d\n",refined[i]->length);

     anchor_thread_args[i].stats = &threads_anchor_rv[i];
     anchor_thread_args[i].list_addr = &(refined[i]);

     thpool_add_work(pool, FragmentRefine, &anchor_thread_args[i]);

     //pthread_create(&threads_anchor[i], NULL, FragmentRefine, &anchor_thread_args[i]);
     //FragmentRefine(&anchor_thread_args[i]);
  }
  //

  thpool_wait(pool);

  /*for (i = 0; i < conf->nthreads; i ++)
    pthread_join(threads_anchor[i], (void **)&threads_anchor_rv[i]);*/

  printf("Entering deduplication stage\n");
  //int qs = 0;
  List * refined_merged = emptylist();
  //printf("Voor merge\n");
  /*while (qs < threadCount){
    //printf("Voor merge. Lengte refined[qs]: %d\n",refined[qs]->length);
    //printf("Lengte refined: %d\n", refined[qs]->length);
    refined_merged = merge(refined_merged,refined[qs]);
    //printf("Na merge\n");
    qs++;
  }*/
  for (i = 0; i < threadCount; i++) refined_merged = merge(refined_merged,refined[i]);
  //printf("Na merge\n");
  /*Node * nnn = refined_merged->head;
  for (int i = 1; i <= refined_merged->length; i++){
    //if(nnn == NULL) printf("Node %d is NULL\n",i);
    //else printf("%d",i);
    nnn = nnn->next;
  }*/
  //printf("Voor split_mod. Totaal aantal chunks: %d\n",refined_merged->length);
  List ** dedup = split_mod(threadCount, refined_merged);
  //printf("Na split_mod\n");
  struct thread_args chunk_thread_args[conf->nthreads];
  stats_t threads_chunk_rv[conf->nthreads];
  for (i = 0; i < conf->nthreads; i ++) {
    chunk_thread_args[i].tid = i;
    chunk_thread_args[i].list = dedup[i];
    //printf("Assigning stats\n");
    chunk_thread_args[i].stats = &threads_chunk_rv[i];
    //printf("Checking sequence to_dedup\n");
    //check_sequence(to_dedup[i]);
    //printf("Sequence to_dedup OK\n");
    //printf("Creating dedup task\n");
    thpool_add_work(pool, Deduplicate, &chunk_thread_args[i]);
    //printf("Created dedup task\n");
    //pthread_create(&threads_chunk[i], NULL, Deduplicate, &chunk_thread_args[i]);

    //Deduplicate(&chunk_thread_args[i]);
  }
  //printf("Created dedup threads\n");
  thpool_wait(pool);
  /*for (i = 0; i < conf->nthreads; i ++)
    pthread_join(threads_chunk[i], (void **)&threads_chunk_rv[i]);*/

    //return;

  printf("Entering compression stage\n");

  //Gather all useful stats by traversing the dedup lists once
  /*int chunks = 0;
  int duplicates = 0;
  int uncompressed_data_size = 0;
  for (i = 0; i < threadCount; i++){
    chunks+=dedup[i]->length;
    Iterator * iter = init_iterator(dedup[i]);
    while(hasNext(iter)){
      chunk_t * c = next(iter);
      if(!c->header.isDuplicate) duplicates++;
      uncompressed_data_size+=c->uncompressed_data.n;
    }
    destroy_iterator(iter);
  }*/
  //return;

  //List ** test = malloc(threadCount * sizeof(List *));
  //for (int i = 0; i<threadCount; i++) test[i] = dedup[i];
  //assign_compression_buffers(duplicates, uncompressed_data_size, dedup);
  //printf("Before zip_split\n");
  List ** compress = zip_split(threadCount, dedup);
  //printf("After zip_split\n");
  Compressed_data * total_compressed_data[threadCount];

  int buffer_counts[threadCount];

  //char ** write_buffers = malloc(sizeof(char *) * (list->length/100 + ((list->length % 100 == 0) ? 1:0)));

  struct thread_args compress_thread_args[conf->nthreads];
  stats_t threads_compress_rv[conf->nthreads];
  for (i = 0; i < conf->nthreads; i ++) {
    compress_thread_args[i].tid = i;
    compress_thread_args[i].list = compress[i];
    int write_buffer_count = compress[i]->length/100 + ((compress[i]->length % 100 == 0) ? 0:1);
    //printf("Lengte buffer: %d\n", write_buffer_count);
    total_compressed_data[i] = malloc(write_buffer_count * sizeof(Compressed_data));
    buffer_counts[i] = write_buffer_count;
    compress_thread_args[i].compressed_data = total_compressed_data[i];
    compress_thread_args[i].stats = &threads_compress_rv[i];
    //printf("Checking sequence compressed before compression\n");
    //check_sequence(compressed[i]);
    //printf("Sequence compressed before compression OK\n");
    thpool_add_work(pool, Compress, &compress_thread_args[i]);

    //pthread_create(&threads_compress[i], NULL, Compress, &compress_thread_args[i]);
    //Compress(&compress_thread_args[i]);
  }
  thpool_wait(pool);
  thpool_destroy(pool);

  /*for (i = 0; i < conf->nthreads; i ++){
    //printf("Joining thread\n");
    pthread_join(threads_compress[i], (void **)&threads_compress_rv[i]);
    //printf("Joined thread\n");
  }*/


  //while(1) printf("Voorbij join\n");
  //return;
  printf("Entering writing stage\n");
  /*List ** dedupped_merged = malloc(threadCount * sizeof(List *));
  for (i = 0; i < threadCount; i++){
    //printf("Checking sequence compressed\n");
    //check_sequence(compressed[i]);
    //printf("Sequence compressed OK\n");
    //printf("Checking sequence dedup_qs\n");
    //check_sequence(dedup_qs[i]->send);
    //printf("Sequence dedup_qs OK\n");
     dedupped_merged[i] = sorted_merge(compressed[i], dedup_qs[i]->send);
   }*/
   //printf("Voorbij sorted_merge\n");

  //printf("Voorbij zip\n");



  Write(total_compressed_data, buffer_counts);


#ifdef ENABLE_PARSEC_HOOKS
  __parsec_roi_end();
#endif

#ifdef ENABLE_STATISTICS
  //Merge everything into global `stats' structure
  for(i=0; i<conf->nthreads; i++) {
    merge_stats(&stats, &threads_anchor_rv[i]);
    //free(threads_anchor_rv[i]);
  }
  for(i=0; i<conf->nthreads; i++) {
    merge_stats(&stats, &threads_chunk_rv[i]);
    //free(threads_chunk_rv[i]);
  }
  for(i=0; i<conf->nthreads; i++) {
    merge_stats(&stats, &threads_compress_rv[i]);
    //free(threads_compress_rv[i]);
  }
#endif //ENABLE_STATISTICS

  // clean up with the src file
  if (conf->infile != NULL)
    close(fd);

  //free(refined);

  //free(refined_merged);

  //Send queues are freed by merge function.
  //for (int i = 0; i < threadCount; i++) free(dedup_qs[i]);
  //free(dedup_qs);


  //Sublists are freed by sorted_merge.
  //free(compressed);

  //free(dedupped_merged);

  //destroy(to_reorder);



  hashtable_destroy(cache, FALSE);

#ifdef ENABLE_STATISTICS
   //dest file stat
  if (stat(conf->outfile, &filestat) < 0)
      EXIT_TRACE("stat() %s failed: %s\n", conf->outfile, strerror(errno));
  stats.total_output = filestat.st_size;

  //Analyze and print statistics
  if(conf->verbose) print_stats(&stats);
#endif //ENABLE_STATISTICS


}
