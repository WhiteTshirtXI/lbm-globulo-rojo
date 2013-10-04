#ifndef _MEMORY_H_
#define _MEMORY_H_

void alloc_memory_GPU(float **data_d, size_t size);
void free_memory_GPU(float *data_d);

void send_data_to_GPU(float *data, float *data_d, size_t size);
void retrieve_data_from_GPU(float *data, float *data_d, size_t size);


#endif
