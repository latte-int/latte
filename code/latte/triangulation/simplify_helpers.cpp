/* simplify_helpers.cpp -- Methods to simplify data used to find
   a subspace avoiding completion vector for a boundary triangulation
           
   Copyright 2011       Christof Soeger

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



//make the first entry positive by *(-1) if needed
static void make_first_entry_positive(vec_ZZ &v) {
  for (int i=0; i<v.length(); i++) {
    if (v[i] > 0) break;
    if (v[i] < 0) {
      v*=(-1);
      break;
    }
  }
}

//compute the gcd of all entries
static ZZ vec_gcd(const vec_ZZ &v) {
  ZZ g = ZZ(); //initialized to 0
  for (int i=0; i<v.length(); i++) {
    g = GCD(v[i],g);
    if (g == 1)
      break;
  }
  return g;
}

static ZZ make_coprime(vec_ZZ &v) {
  ZZ g = vec_gcd(v);
  if (g != 0) {
    for (int i=0; i<v.length(); i++) {
      v[i] /= g;
    }
  }
  return g;
}

//works only for vectors of same length
static bool vec_ZZ_is_less(const vec_ZZ &a, const vec_ZZ &b) {
  assert (a.length() == b.length());
  for (int i=0; i<a.length(); i++) {
    if (a[i] < b[i]) return true;
    if (a[i] > b[i]) return false;
  }
  return false;
}

//find the first non-zero position in the vector
static int first_non_zero_pos(const vec_ZZ& v) {
  int i;
  for (i=0; i<v.length() && v[i]==0; i++);
  return i;
}


//remove rows for which a multiple of it is also in the list
void remove_dominated_rows(list<vec_ZZ>& a_list) {
  list<vec_ZZ>::iterator it = a_list.begin();
  list<vec_ZZ>::const_iterator cit;
  int fnz; //first non zero position
  ZZ factor1, factor2;
  while (it!=a_list.end()) {
    fnz = first_non_zero_pos(*it);
    cit = it; ++cit;
    while (cit!=a_list.end() && fnz == first_non_zero_pos(*cit)) {
      factor1 = (*it)[fnz];
      factor2 = (*cit)[fnz];
      if ( factor1 <= factor2 && (*it)*factor2 == (*cit)*factor1 ) {
        fnz = -1; //to mark success
        it = a_list.erase(it);
        break;
      }
	  ++cit;
    }
    if (fnz >= 0) ++it;
  }
}
