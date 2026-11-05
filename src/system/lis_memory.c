/* Copyright (C) 2005 The Scalable Software Infrastructure Project. All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:
   1. Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
   2. Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
   3. Neither the name of the project nor the names of its contributors 
      may be used to endorse or promote products derived from this software 
      without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE SCALABLE SOFTWARE INFRASTRUCTURE PROJECT
   ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
   TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
   PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE SCALABLE SOFTWARE INFRASTRUCTURE
   PROJECT BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
   OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
*/

#ifdef HAVE_CONFIG_H
	#include "lis_config.h"
#else
#ifdef HAVE_CONFIG_WIN_H
	#include "lis_config_win.h"
#endif
#endif
/*#undef malloc*/

#include <stdlib.h>
#ifdef HAVE_MALLOC_H
        #include <malloc.h>
#endif
#include <string.h>
#include <stdarg.h>
#include "lislib.h"

/************************************************
 * *rpl_malloc
 * *lis_malloc
 * *lis_calloc
 * *lis_realloc
 * lis_free
 * lis_free_mat
 * lis_free2
 * lis_free_all
 * lis_is_malloc
 * lis_free2
 * lis_is_malloc
 ************************************************/

#ifdef _DEBUG
	#define USE_MALLOC_TAG	1
#endif

#if 1
typedef struct _malloc_address {
	struct _malloc_address *next, *prev;
	void *address, *real_address;
	size_t size;
	#ifdef USE_MALLOC_TAG
		char *tag;
	#endif
} malloc_address;

malloc_address malloc_address_top = 
  {&malloc_address_top, &malloc_address_top, NULL, NULL};
#else
#define MALLOC_ADDRESS_SIZE	1021

typedef struct _malloc_address {
  struct _malloc_address *next, *prev;
  void *address, *real_address;
  size_t size;
} malloc_address;

malloc_address malloc_address_top = 
  {&malloc_address_top, &malloc_address_top, NULL, NULL};

malloc_address *malloc_address_table;
#endif

/*
void *malloc ();
 
void *rpl_malloc( size_t n )
{
	if( n == 0 )
	n = 1;
	return malloc( n );
}
*/


/*
 * Tracked allocations use 16-byte alignment. Publish a pointer in the
 * allocation list only after both the storage and its metadata exist.
 */
static void *lis_alloc_tracked(size_t size, char *tag, int clear)
{
	union {
		void *ptr;
		unsigned long l;
	} addr;
	malloc_address *ma;
	void *real_address;

	if( size > (size_t)-1 - 16 ) return NULL;
	real_address = malloc(size + 16);
	if( real_address==NULL ) return NULL;

	ma = (malloc_address *)malloc(sizeof(malloc_address));
	if( ma==NULL )
	{
		free(real_address);
		return NULL;
	}
	#ifdef USE_MALLOC_TAG
	{
		const char *name = tag ? tag : "";
		ma->tag = (char *)malloc(strlen(name)+1);
		if( ma->tag==NULL )
		{
			free(ma);
			free(real_address);
			return NULL;
		}
		strcpy(ma->tag,name);
	}
	#else
	(void)tag;
	#endif

	addr.ptr = real_address;
	addr.l += 0xf;
	addr.l &= ~0xf;
	ma->address = addr.ptr;
	ma->real_address = real_address;
	ma->size = size;
	if( clear ) memset(addr.ptr,0,size);

	ma->next = &malloc_address_top;
	ma->prev = malloc_address_top.prev;
	ma->prev->next = ma;
	ma->next->prev = ma;
	return addr.ptr;
}

void *lis_malloc( size_t size, char *tag )
{
	return lis_alloc_tracked(size,tag,0);
}

void *lis_calloc( size_t size, char *tag )
{
	return lis_alloc_tracked(size,tag,1);
}

void *lis_realloc( void *p, size_t size )
{
	union {
		void *ptr;
		unsigned long l;
	} addr;
	malloc_address *ma;
	void *real_address;
	size_t copy_size;

	ma = malloc_address_top.next;
	while (ma->address)
	{
		if (ma->address == p) break;
		ma = ma->next;
	}
	if (ma->address)
	{
		/* Preserve the original tracked allocation on failure. */
		if( size > (size_t)-1 - 16 ) return NULL;
		real_address = malloc(size + 16);
		if( real_address==NULL ) return NULL;
		addr.ptr = real_address;
		addr.l += 0xf;
		addr.l &= ~0xf;

		copy_size = ma->size < size ? ma->size : size;
		memcpy(addr.ptr,p,copy_size);
		free(ma->real_address);

		ma->address      = addr.ptr;
		ma->real_address = real_address;
		ma->size         = size;
		return addr.ptr;
	}
	/* Preserve standard realloc behaviour for untracked pointers. */
	return realloc(p,size);
}

void lis_free(void *p)
{
  malloc_address *ma;

  ma = malloc_address_top.next;
  while (ma->address) {
    if (ma->address == p) break;
    ma = ma->next;
  }
  if (ma->address) {
    ma->next->prev = ma->prev;
    ma->prev->next = ma->next;
    free(ma->real_address);
#ifdef USE_MALLOC_TAG
    free(ma->tag);
#endif
    free(ma); 
    return;
  }
  free(p);
  return;
}

void lis_free_mat(LIS_MATRIX A)
{
	malloc_address *ma;
	malloc_address *ma_bak;
	LIS_INT i;

	ma = malloc_address_top.next;
	for(i=0;i<A->n;i++)
	{
		if( A->w_index[i]==NULL ) continue;
		while (ma->address) {
			if (ma->address == A->w_index[i]) break;
			ma = ma->next;
		}
		if (ma->address) {
			ma->next->prev = ma->prev;
			ma->prev->next = ma->next;
			ma_bak         = ma->next;
			free(ma->real_address);
			#ifdef USE_MALLOC_TAG
				free(ma->tag);
			#endif
			free(ma);
			ma = ma_bak;
		}
		while (ma->address) {
			if (ma->address == A->w_value[i]) break;
			ma = ma->next;
		}
		if (ma->address) {
			ma->next->prev = ma->prev;
			ma->prev->next = ma->next;
			ma_bak         = ma->next;
			free(ma->real_address);
			#ifdef USE_MALLOC_TAG
				free(ma->tag);
			#endif
			free(ma);
			ma = ma_bak;
		}
	}
	return;
}

void lis_free2(LIS_INT n, ...)
{
	va_list vvlist;
	LIS_INT	i;
	char *p;

	va_start( vvlist, n );
	for(i=0;i<n;i++)
	{
		p = va_arg(vvlist, char *);
		if( p )
		{
			lis_free(p);
		}
	}
	va_end(vvlist);
}

void lis_free_all()
{
	malloc_address *ma;
	malloc_address *ma_bak;

	ma = malloc_address_top.next;
	while (ma->address) {
#ifdef USE_MALLOC_TAG
		lis_printf(LIS_COMM_WORLD,"memory leak: address = %x size=%d (%s)\n",ma->address,ma->size,ma->tag);
#endif
		ma_bak = ma->next;
		/*  if (ma->real_address) free(ma->real_address); */
		ma->real_address = NULL;
#ifdef USE_MALLOC_TAG
		if (ma->tag) lis_free(ma->tag);
		ma->tag = NULL;
#endif
		ma = ma_bak;
	}
	malloc_address_top.next         = &malloc_address_top;
	malloc_address_top.prev         = &malloc_address_top;
	malloc_address_top.address      = NULL;
	malloc_address_top.real_address = NULL;
	malloc_address_top.size         = 0;
	return;
}

LIS_INT lis_is_malloc( void *p )
{
	malloc_address *ma;

	ma = malloc_address_top.next;
	while (ma->address)
	{
	    if (ma->address == p) break;
		ma = ma->next;
	}
	if (ma->address)
	{
		return LIS_TRUE;
	}
	return LIS_FALSE;
}
