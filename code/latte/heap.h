/* heap.cpp -- Binary heaps

   Copyright 2003 Utz-Uwe Haus
   Copyright 2004, 2006, 2007 Matthias Koeppe

   This file is part of LattE.
   
   LattE is free software; you can redistribute it and/or modify it
   under the terms of the version 2 of the GNU General Public License
   as published by the Free Software Foundation.

   LattE is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with LattE; if not, write to the Free Software Foundation,
   Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
*/

// $GYWOPT: misc/heap.h,v 1.6 2004/12/09 12:49:18 durzinsk Exp $
 
#ifndef  HEAP_H_
#define  HEAP_H_     /*+ To stop multiple inclusions. +*/

#include <cstddef>

/* heaps have one of the 4 default or a custom comparisons available */
enum HEAP_TYPE {
     HEAP_INT_MIN = 1,
     HEAP_INT_MAX,
     HEAP_DBL_MIN,
     HEAP_DBL_MAX,
     HEAP_CUSTOM
 };

/* If you use a custom heap comparison, you must specify a comparison function.
 * It allows you to build a MIN-first heap. Comparison as for strcmp: <0, 0, >0
 * depending on whether d1<d2, d1==d2 or d1>d2
 */
typedef int (*heap_cmp_fun) (const void* d1, const void* d2);

struct heap_component_int {
    int    weight;
    void  *data;
};
struct heap_component_dbl {
    double weight;
    void  *data;
};
struct heap_component_cust {
    void  *weight;
    void  *data;
};

struct heap {
    enum HEAP_TYPE type;
    size_t size;
    size_t next_pos;
    
    heap_cmp_fun custom_compare;

    bool frozen;

    union {
	/* depending on TYPE one of the following self-extending structs is active */
	struct heap_component_int  int_heap[1];
	struct heap_component_dbl  dbl_heap[1];
	struct heap_component_cust cust_heap[1];
    } u;
};

/* allocate a new heap of specify type. If type == HEAP_CUSTOM, you must specify
   the comparison function. Otherwise that argument is ignored */

struct heap * 
heap_alloc(enum HEAP_TYPE type, size_t size, heap_cmp_fun custom_compare_fun);

/* freeze a heap for mass-insertions without immediate heapification */ 
void
heap_freeze(struct heap *h);

/* thaw a frozen heap and heapify it immediately */
void
heap_thaw(struct heap *h);

/* retrieve top element without removing it from the heap */
void *
heap_top(struct heap *h);

/* retrieve top element */
void *
heap_pop(struct heap *h);

/* insert, 3 versions depending on heap type */
void
heap_insert_intweight(struct heap *h, void *data, int weight);
void
heap_insert_dblweight(struct heap *h, void *data, double weight);
void
heap_insert_custweight(struct heap *h, void *data, void *weight);

/* change weight of top element */ 
void
heap_change_top_intweight(struct heap *h, int weight);

/* if you just want to store an elementary data type in the void* use */
#define INT_TO_POINTER(i)   ((void*) (i))
#define POINTER_TO_INT(p)   ((int)   (p))
#define UINT_TO_POINTER(ui) ((void*) (ui))
#define POINTER_TO_UINT(p)  ((unsigned int) (p))

/* clear a heap for re-use. Just forgets, does not free any of your data */
void
heap_clear(struct heap *h);

/* get rid of the heap H */
void
heap_destroy(struct heap *h);

bool
heap_empty_p(const struct heap *h);
bool
heap_nonempty_p(const struct heap *h);

size_t
heap_num_fill(const struct heap *h);

#endif  /* HEAP_H_ */
