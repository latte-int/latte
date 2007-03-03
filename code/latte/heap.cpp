/* heap.cpp -- Binary heaps

   Copyright 2003 Utz-Uwe Haus
   Copuright 2004, 2006, 2007 Matthias Koeppe

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

// $GYWOPT: misc/heap.c,v 1.8 2004/08/16 17:08:54 mkoeppe Exp $

#include <cassert>
#include <cstdlib>
#include "heap.h"

/******************* local helpers ******************/

/* all kinds of percolating stuff */
static void
percolate_down_int_min(struct heap *h, size_t pos)
{
    size_t child;
    size_t tmp_pos;

    tmp_pos = pos;
    child = (tmp_pos << 1) + 1;
    while(child < h->next_pos) {
	/* compare with neighbors */
	if(child < h->next_pos - 1 /* still have two children */
	   && h->u.int_heap[child+1].weight < h->u.int_heap[child].weight)
	    child += 1;
	
	if(h->u.int_heap[tmp_pos].weight > h->u.int_heap[child].weight) {
	    struct heap_component_int xxx = h->u.int_heap[tmp_pos];
	    h->u.int_heap[tmp_pos]        = h->u.int_heap[child];
	    h->u.int_heap[child]          = xxx;
	    tmp_pos = child;
	    child = (tmp_pos<<1) + 1;
	} else
	    break;
    }
}
static void
percolate_down_int_max(struct heap *h, size_t pos)
{
    size_t child;
    size_t tmp_pos;

    tmp_pos = pos;
    child = (tmp_pos << 1) + 1;
    while(child < h->next_pos) {
	/* compare with neighbors */
	if(child < h->next_pos - 1 /* still have two children */
	   && h->u.int_heap[child+1].weight > h->u.int_heap[child].weight)
	    child += 1;
	
	if(h->u.int_heap[tmp_pos].weight < h->u.int_heap[child].weight) {
	    struct heap_component_int xxx = h->u.int_heap[tmp_pos];
	    h->u.int_heap[tmp_pos]        = h->u.int_heap[child];
	    h->u.int_heap[child]          = xxx;
	    tmp_pos = child;
	    child = (tmp_pos<<1) + 1;
	} else
	    break;
    }
}
static void
percolate_down_dbl_min(struct heap *h, size_t pos)
{
    size_t child;
    size_t tmp_pos;

    tmp_pos = pos;
    child = (tmp_pos << 1) + 1;
    while(child < h->next_pos) {
	/* compare with neighbors */
	if(child < h->next_pos - 1 /* still have two children */
	   && h->u.dbl_heap[child+1].weight < h->u.dbl_heap[child].weight)
	    child += 1;
	
	if(h->u.dbl_heap[tmp_pos].weight > h->u.dbl_heap[child].weight) {
	    struct heap_component_dbl xxx = h->u.dbl_heap[tmp_pos];
	    h->u.dbl_heap[tmp_pos]        = h->u.dbl_heap[child];
	    h->u.dbl_heap[child]          = xxx;
	    tmp_pos = child;
	    child = (tmp_pos<<1) + 1;
	} else
	    break;
    }
}
static void
percolate_down_dbl_max(struct heap *h, size_t pos)
{
    size_t child;
    size_t tmp_pos;

    tmp_pos = pos;
    child = (tmp_pos << 1) + 1;
    while(child < h->next_pos) {
	/* compare with neighbors */
	if(child < h->next_pos - 1 /* still have two children */
	   && h->u.dbl_heap[child+1].weight > h->u.dbl_heap[child].weight)
	    child += 1;
	
	if(h->u.dbl_heap[tmp_pos].weight < h->u.dbl_heap[child].weight) {
	    struct heap_component_dbl xxx = h->u.dbl_heap[tmp_pos];
	    h->u.dbl_heap[tmp_pos]        = h->u.dbl_heap[child];
	    h->u.dbl_heap[child]          = xxx;
	    tmp_pos = child;
	    child = (tmp_pos<<1) + 1;
	} else
	    break;
    }
}
static void
percolate_down_cust(struct heap *h, size_t pos)
{
    size_t child;
    size_t tmp_pos;

    tmp_pos = pos;
    child = (tmp_pos << 1) + 1;
    while(child < h->next_pos) {
	/* compare with neighbors */
	if(child < h->next_pos - 1 /* still have two children */
	   && h->custom_compare(h->u.cust_heap[child+1].weight,
				h->u.cust_heap[child].weight) < 0)
	    child += 1;
	
	if(h->custom_compare(h->u.cust_heap[tmp_pos].weight,
			     h->u.cust_heap[child].weight) > 0) {
	    struct heap_component_cust xxx = h->u.cust_heap[tmp_pos];
	    h->u.cust_heap[tmp_pos]        = h->u.cust_heap[child];
	    h->u.cust_heap[child]          = xxx;
	    tmp_pos = child;
	    child = (tmp_pos<<1) + 1;
	} else
	    break;
    }
}

/* All kinds of bubbles */
static void
bubble_up_int_min(struct heap *h, size_t pos)
{
    size_t temp_pos = pos;
    while(temp_pos > 0) {
	size_t parent = (temp_pos - 1) >> 1;
	if(h->u.int_heap[temp_pos].weight < h->u.int_heap[parent].weight) {
	    /* hard hat area, heavy structures moving here */
	    struct heap_component_int xxx =  h->u.int_heap[temp_pos];
	    h->u.int_heap[temp_pos] = h->u.int_heap[parent];
	    h->u.int_heap[parent] = xxx;
	    temp_pos = parent;
	} else
	    temp_pos = 0;
    }
}
static void
bubble_up_int_max(struct heap *h, size_t pos)
{
    size_t temp_pos = pos;
    while(temp_pos > 0) {
	size_t parent = (temp_pos - 1) >> 1;
	if(h->u.int_heap[temp_pos].weight > h->u.int_heap[parent].weight) {
	    /* hard hat area, heavy structures moving here */
	    struct heap_component_int xxx =  h->u.int_heap[temp_pos];
	    h->u.int_heap[temp_pos] = h->u.int_heap[parent];
	    h->u.int_heap[parent] = xxx;
	    temp_pos = parent;
	} else
	    temp_pos = 0;
    }
}
static void
bubble_up_dbl_min(struct heap *h, size_t pos)
{
    size_t temp_pos = pos;
    while(temp_pos > 0) {
	size_t parent = (temp_pos - 1) >> 1;
	if(h->u.dbl_heap[temp_pos].weight < h->u.dbl_heap[parent].weight) {
	    /* hard hat area, heavy structures moving here */
	    struct heap_component_dbl xxx =  h->u.dbl_heap[temp_pos];
	    h->u.dbl_heap[temp_pos] = h->u.dbl_heap[parent];
	    h->u.dbl_heap[parent] = xxx;
	    temp_pos = parent;
	} else
	    temp_pos = 0;
    }
}
static void
bubble_up_dbl_max(struct heap *h, size_t pos)
{
    size_t temp_pos = pos;
    while(temp_pos > 0) {
	size_t parent = (temp_pos - 1) >> 1;
	if(h->u.dbl_heap[temp_pos].weight > h->u.dbl_heap[parent].weight) {
	    /* hard hat area, heavy structures moving here */
	    struct heap_component_dbl xxx =  h->u.dbl_heap[temp_pos];
	    h->u.dbl_heap[temp_pos] = h->u.dbl_heap[parent];
	    h->u.dbl_heap[parent] = xxx;
	    temp_pos = parent;
	} else
	    temp_pos = 0;
    }
}
static void
bubble_up_cust(struct heap *h, size_t pos)
{
    size_t temp_pos = pos;
    assert(h->custom_compare!=NULL);

    while(temp_pos > 0) {
	size_t parent = (temp_pos - 1) >> 1;
	if(h->custom_compare(h->u.cust_heap[temp_pos].weight,
			     h->u.cust_heap[parent].weight) < 0) {
	    /* hard hat area, heavy structures moving here */
	    struct heap_component_cust xxx =  h->u.cust_heap[temp_pos];
	    h->u.cust_heap[temp_pos] = h->u.cust_heap[parent];
	    h->u.cust_heap[parent] = xxx;
	    temp_pos = parent;
	} else
	    temp_pos = 0;
    }
}

/* heapify a formerly frozen heap as smart as we can */
static void
heapify(struct heap *h)
{
    /* about half of the elements are suitable as children. */
    if(h->next_pos <= 1)
	return; /* nothing to be done: empty or 1 element */
    else {
	size_t next_insert = ( (h->next_pos-1) /* last element we have */
			       - 1) >> 1;      /* its parent */
	
	switch (h->type) {
	  case HEAP_INT_MIN:
	      do {
		  percolate_down_int_min(h, next_insert);
	      } while (next_insert-- > 0);
	      break;
	  case HEAP_INT_MAX:
	      do {
		  percolate_down_int_max(h, next_insert);
	      } while (next_insert-- > 0);
	      break;
	  case HEAP_DBL_MIN:
	      do {
		  percolate_down_dbl_min(h, next_insert);
	      } while (next_insert-- > 0);
	      break;
	  case HEAP_DBL_MAX:
	      do {
		  percolate_down_dbl_max(h, next_insert);
	      } while (next_insert-- > 0);
	      break;
	  case HEAP_CUSTOM:
	      do {
		  percolate_down_cust(h, next_insert);
	      } while (next_insert-- > 0);
	      break;
	  default:
	      assert(h->type == HEAP_INT_MIN || h->type == HEAP_INT_MAX
			 || h->type == HEAP_DBL_MIN || h->type == HEAP_DBL_MAX
			 || h->type ==HEAP_CUSTOM);
	      return;
	}
    }
}


/* externally visible stuff */

struct heap *
heap_alloc(enum HEAP_TYPE type, size_t size, heap_cmp_fun custom_compare_fun)
{
    struct heap * result;
    size_t alloc_size;
    
    switch (type) {
      case HEAP_INT_MAX:
      case HEAP_INT_MIN:
	  alloc_size = (sizeof(struct heap)
			+ (size-1) * sizeof(struct heap_component_int));
	  break;
      case HEAP_DBL_MAX:
      case HEAP_DBL_MIN:
	  alloc_size = (sizeof(struct heap)
			+ (size-1) * sizeof(struct heap_component_dbl));
	  break;
      case HEAP_CUSTOM:
	  alloc_size = (sizeof(struct heap)
			+ (size-1) * sizeof(struct heap_component_cust));
	  break;
      default:
	  assert(0);
	  return NULL;
    };

    result = (struct heap *) malloc(alloc_size);
    if(result!=NULL) {
	result->type     = type;
	result->size     = size;
	result->next_pos = 0;
	if(type == HEAP_CUSTOM) 
	    result->custom_compare = custom_compare_fun;
	else
	    result->custom_compare = NULL;
	result->frozen   = false;
    }

    return result;
}

/* freeze a heap for mass-insertions without immediate heapification */ 
void
heap_freeze(struct heap *h)
{
    assert(h!=NULL);
    if(h->frozen) {
      assert(0);
    } else
	h->frozen=true;
}

/* thaw a frozen heap and heapify it immediately */
void
heap_thaw(struct heap *h)
{
    assert(h!=NULL);

    if(!h->frozen) {
      assert(0);
    } else {
	heapify(h);
	h->frozen=false;
    }
}

/* retrieve top element without removing it from the heap */
void *
heap_top(struct heap *h)
{
    void *data;
    assert(h!=NULL);
    assert(h->next_pos>0);
    assert(! h->frozen);
    
    switch (h->type) {
      case HEAP_INT_MIN:
      case HEAP_INT_MAX:
	  data = h->u.int_heap[0].data;
	  break;
      case HEAP_DBL_MIN:
      case HEAP_DBL_MAX:
	  data = h->u.dbl_heap[0].data;
	  break;
      case HEAP_CUSTOM:
	  data = h->u.cust_heap[0].data;
	  break;
      default:
	assert(0);
	  return NULL;
    }
    
    return data;
}

/* retrieve top element */
void *
heap_pop(struct heap *h)
{
    void *data;
    assert(h!=NULL);
    assert(h->next_pos>0);
    assert(! h->frozen);
    
    switch (h->type) {
      case HEAP_INT_MIN:
	  data = h->u.int_heap[0].data;
	  h->next_pos--;
	  h->u.int_heap[0].weight = h->u.int_heap[h->next_pos].weight;
	  h->u.int_heap[0].data   = h->u.int_heap[h->next_pos].data;
	  percolate_down_int_min(h, 0);
	  break;
      case HEAP_INT_MAX:
	  data = h->u.int_heap[0].data;
	  h->next_pos--;
	  h->u.int_heap[0].weight = h->u.int_heap[h->next_pos].weight;
	  h->u.int_heap[0].data   = h->u.int_heap[h->next_pos].data;
	  percolate_down_int_max(h, 0);
	  break;
      case HEAP_DBL_MIN:
	  data = h->u.dbl_heap[0].data;
	  h->next_pos--;
	  h->u.dbl_heap[0].weight = h->u.dbl_heap[h->next_pos].weight;
	  h->u.dbl_heap[0].data   = h->u.dbl_heap[h->next_pos].data;
	  percolate_down_dbl_min(h, 0);
	  break;
      case HEAP_DBL_MAX:
	  data = h->u.dbl_heap[0].data;
	  h->next_pos--;
	  h->u.dbl_heap[0].weight = h->u.dbl_heap[h->next_pos].weight;
	  h->u.dbl_heap[0].data   = h->u.dbl_heap[h->next_pos].data;
	  percolate_down_dbl_max(h, 0);
	  break;
      case HEAP_CUSTOM:
	  data = h->u.cust_heap[0].data;
	  h->next_pos--;
	  h->u.cust_heap[0].weight = h->u.cust_heap[h->next_pos].weight;
	  h->u.cust_heap[0].data   = h->u.cust_heap[h->next_pos].data;
	  percolate_down_cust(h, 0);
	  break;
      default:
	assert(0);
	  return NULL;
    }
    
    return data;
}

/* insert, 3 versions depending on heap type */
void
heap_insert_intweight(struct heap *h, void *data, int weight)
{
    assert(h->next_pos < h->size);
    assert(h->type == HEAP_INT_MIN || h->type == HEAP_INT_MAX);
    
    h->u.int_heap[h->next_pos].weight = weight;
    h->u.int_heap[h->next_pos].data   = data;

    if(!h->frozen) {
	if(h->type == HEAP_INT_MIN)
	    bubble_up_int_min(h, h->next_pos);
	else
	    bubble_up_int_max(h, h->next_pos);
    }
    h->next_pos++;
}
    
void
heap_insert_dblweight(struct heap *h, void *data, double weight)
{
    assert(h->next_pos < h->size);
    assert(h->type == HEAP_DBL_MIN || h->type == HEAP_DBL_MAX);
    
    h->u.dbl_heap[h->next_pos].weight = weight;
    h->u.dbl_heap[h->next_pos].data   = data;

    if(!h->frozen) {
	if(h->type == HEAP_DBL_MIN)
	    bubble_up_dbl_min(h, h->next_pos);
	else
	    bubble_up_dbl_max(h, h->next_pos);
    }
    h->next_pos++;
}

void
heap_insert_custweight(struct heap *h, void *data, void *weight)
{
    assert(h->next_pos < h->size);
    assert(h->type == HEAP_CUSTOM);
    
    h->u.cust_heap[h->next_pos].weight = weight;
    h->u.cust_heap[h->next_pos].data   = data;

    if(!h->frozen)
	bubble_up_cust(h, h->next_pos);

    h->next_pos++;
}

/* change weight of top element */ 

void
heap_change_top_intweight(struct heap *h, int weight)
{
    assert(h->next_pos != 0);
    assert(h->type == HEAP_INT_MIN || h->type == HEAP_INT_MAX);

    h->u.int_heap[0].weight = weight;

    if(!h->frozen) {
	if(h->type == HEAP_INT_MIN)
	    percolate_down_int_min(h, 0);
	else
	    percolate_down_int_max(h, 0);
    }
}

/* clear a heap for re-use. Just forgets, does not free any of your data */
void
heap_clear(struct heap *h)
{
    assert(h!=NULL);
    h->next_pos = 0;
}

void
heap_destroy(struct heap *h)
{
    assert(h!=NULL);
    free(h);
}

bool
heap_empty_p(const struct heap *h)
{
    assert(h!=NULL);
    return h->next_pos==0;
}
bool
heap_nonempty_p(const struct heap *h)
{
    assert(h!=NULL);
    return h->next_pos!=0;
}

size_t
heap_num_fill(const struct heap *h)
{
    assert(h!=NULL);
    return h->next_pos;
}
