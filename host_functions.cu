#include "declarations.cuh"

//allocate memory on host
void hostInitialize(DataArray& host)
{
	int size = total * sizeof(float);
	host.a = (float*) malloc(size);

	if (host.a == NULL) Error("Host memory allocation failed\n");
}

//free memory on host
void hostClose(DataArray& host)
{
	free(host.a);
}
