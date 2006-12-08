;;; latte.el --- Emacs editing helpers for the LattE project
;;; 
;;; Copyright 2006 Matthias Koeppe
;;;
;;; This file is part of LattE.
;;;
;;; LattE is free software; you can redistribute it and/or modify it
;;; under the terms of the version 2 of the GNU General Public License
;;; as published by the Free Software Foundation.
;;;
;;; LattE is distributed in the hope that it will be useful, but
;;; WITHOUT ANY WARRANTY; without even the implied warranty of
;;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
;;; General Public License for more details.
;;;
;;; You should have received a copy of the GNU General Public License
;;; along with LattE; if not, write to the Free Software Foundation,
;;; Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.

(defconst ntl-types
  '("ZZ" "vec_ZZ" "mat_ZZ"
    "RR" "vec_RR" "mat_RR"))

(defconst gmp-types
  '("mpz_class" "mpq_class"))

(defconst lidia-types
  '("bigint" "bigint_matrix"))

(defconst cddlib-types
  '("dd_MatrixPtr" "dd_rowset" "dd_colset" "dd_PolyhedraPtr"
    "dd_ErrorType" "dd_SetFamilyPtr"))

(defconst latte-types
  '("listCone" "listVector"
    "vector" "rationalVector"
    "PtrCone" "Cone"
    "Integer"
    "mpq_vector" "mpz_vector"
    "PointsInParallelepipedGenerator" "IntCombEnum"
    "NotGenericException"
    "BarvinokParameters"
    "Vertex" "Polyhedron"))

(setq c++-font-lock-extra-types
      (cons (regexp-opt (append ntl-types latte-types gmp-types lidia-types cddlib-types))
	    c++-font-lock-extra-types))

(if load-file-name
    (let ((directory (file-name-directory load-file-name)))
      (if (file-exists-p (concat directory "TAGS"))
	  (visit-tags-table directory))))

(add-to-list 'auto-insert-alist 
	     '(("\\.\\([Hh]\\|hh\\|hpp\\)\\'" . "C++ header")
	       "Short description: " 
	       (progn (c++-mode) 
		      "// This is a -*- C++ -*- header file.

/* ") (file-name-nondirectory buffer-file-name) " -- " str "
	       
   Copyright " (substring (current-time-string) -4) " " (user-full-name) "

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

"
	       "#ifndef " (upcase (concat (file-name-nondirectory (substring buffer-file-name 0 (match-beginning 0)))
			       "_" 
			       (substring buffer-file-name (1+ (match-beginning 0)))))
	       "
"
	       "#define " (upcase (concat (file-name-nondirectory (substring buffer-file-name 0 (match-beginning 0)))
			       "_" 
			       (substring buffer-file-name (1+ (match-beginning 0)))))
	       "

" _ "

#endif"))

(add-to-list 'auto-insert-alist 
	     '(("\\.\\(cpp\\)\\'" . "C++ source file")
	       "Short description: " 
	       "/* "
	       (file-name-nondirectory buffer-file-name) " -- " str "
	       
   Copyright " (substring (current-time-string) -4) " " (user-full-name) "

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

"))

