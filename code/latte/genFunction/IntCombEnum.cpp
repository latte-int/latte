#include "IntCombEnum.h"
#include <stdlib.h> 
#include <string.h> 
#include <iostream>

using namespace std;

/*
 * use initialization list which is faster to set len = l
 */
IntCombEnum::IntCombEnum(int *u, int l) : len(l), cur_col(l - 1),
   upper_bound(u) {
   /*
    * initialize memory for the start position and "next" position
    * set start and end position to be the zero vector.
    */
   prev = new int[len];
   next = new int[len];
   memset(prev, 0, len*sizeof(int));
   memset(next, 0, len*sizeof(int));

   print_debug();
}

IntCombEnum::~IntCombEnum() {
   /* upper_bound is deleted by creator */ 
   delete [] prev;
   delete [] next;
}

/* sets zeros from column lower_col to column upper_col */    
void IntCombEnum::set_zero(int *v, int lower_col, int upper_col) {
   for (int i = lower_col; i <= upper_col; i++) {
      v[i] = 0;
   }
} 

int *IntCombEnum::getNext() {
   /* check if there are any more integer combinations */
   if (is_last()) {
      cout << "IntCombEnum::getNext -- found last integer combination.\n";
      return NULL;
   }
   /* set next initially equal to the prev iteration */
   copy_comb(next, prev);
   if (prev[cur_col] == upper_bound[cur_col]) {
      while (cur_col >= 0) {
         /* sets all cells after cur_col to 0 */
         set_zero(next, cur_col, (len - 1));
         /* decrement cur_col */
         cur_col--;
         if ((cur_col >= 0) && (prev[cur_col] < upper_bound[cur_col])) {
            next[cur_col]++;
            cur_col = (len - 1);
            break;
         } 
      } 
   } else {
      next[cur_col]++;
   }
   copy_comb(prev, next);
   print_debug();

   return next;
}

void IntCombEnum::decrementUpperBound() {
   for (int i = 0; i < len; i++) {
      upper_bound[i]--;
   }
}

void IntCombEnum::copy_comb(int *dest, int *src) {
   for (int i = 0; i < len; i++) {
      dest[i] = src[i];
   }
}

int IntCombEnum::is_last() {
   int is_zero = 1; 
   /*
    * The enumerator is finished when the prev integer combination is all
    * 0s AND the active column is -1.
    */  
   for (int i = 0; i < len; i++) {
      if (prev[i] != 0) {
         is_zero = 0;
      }
   }
   if ((is_zero) && (cur_col < 0)) {
      return (1);
   }
   return (0);
}

void IntCombEnum::print_debug() {
   cout << "IntCombEnum:: Begin print_debug...\n";
   cout << "len = " << len << ",cur_col = " << cur_col << "\n";
   cout << "next = ";
   for (int i = 0; i < len; i++) {
      cout << next[i] << ",";
   }
   cout << "\n";
   cout << "prev = ";
   for (int i = 0; i < len; i++) {
      cout << prev[i] << ",";
   }
   cout << "\n";
   cout << "upper_bound = ";
   for (int i = 0; i < len; i++) {
      cout << upper_bound[i] << ",";
   }
   cout << "\n";
   cout << "IntCombEnum:: End print_debug...\n";
}
