
/*                                                                    *
 * Author: Ruriko Yoshida                                             *
 * Date: October 25th, 2003                                           *
 * Update: October 25th, 2003                                         *
 * This is for checking redundant inequalities and hidden equations.  *
 *
 */

#ifndef CHECKREDUNDANT__H
#define CHECKREDUNDANT__H

#include "myheader.h"
#include "ramon.h"

listVector* CheckRedIneq(mat_ZZ S);
listVector* CheckHidEqs(mat_ZZ S, listVector* equ);

#endif
