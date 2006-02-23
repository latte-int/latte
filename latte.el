;;; latte.el --- Emacs editing stuff for the LattE project

(defconst ntl-types
  '("ZZ" "vec_ZZ" "mat_ZZ"))

(defconst latte-types
  '("listCone" "listVector"))

(setq c++-font-lock-extra-types
      (cons (regexp-opt (append ntl-types latte-types))
	    c++-font-lock-extra-types))
