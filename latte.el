;;; latte.el --- Emacs editing stuff for the LattE project

(defconst ntl-types
  '("ZZ" "vec_ZZ" "mat_ZZ"))

(defconst gmp-types
  '("mpz_class" "mpq_class"))

(defconst latte-types
  '("listCone" "listVector"
    "vector" "rationalVector"
    "Integer"
    "mpq_vector" "mpz_vector"))

(setq c++-font-lock-extra-types
      (cons (regexp-opt (append ntl-types latte-types gmp-types))
	    c++-font-lock-extra-types))
