#include "chunk_list.h"

#ifndef _ENCODER_H_
#define _ENCODER_H_ 1



#endif /* !_ENCODER_H_ */

typedef struct{
  size_t size;
  char * data;
} Compressed_data;

void Encode(config_t * conf);
