;;; latte.el --- Emacs editing stuff for the LattE project

(defconst ntl-types
  '("ZZ" "vec_ZZ" "mat_ZZ"
    "RR" "vec_RR" "mat_RR"))

(defconst gmp-types
  '("mpz_class" "mpq_class"))

(defconst lidia-types
  '("bigint" "bigint_matrix"))

(defconst cddlib-types
  '("dd_MatrixPtr" "dd_rowset" "dd_colset"))

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
