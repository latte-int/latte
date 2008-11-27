#| (asdf:operate 'asdf:load-op :atrop) 
|#

(in-package :cl-user)

(defun make-vector (dim &key (initial-element 0))
  (make-array dim
	      :initial-element initial-element
	      :adjustable nil :fill-pointer nil))

(defun unit-vector (entry dim)
  ;; 0-based
  (let ((vec (make-vector dim :initial-element 0)))
    (setf (aref vec entry) 1)
    vec))

(defun cyclotomic-polytope/prime (p)
  (let ((dim (- p 1)))
    (poly:make-v-poly dim
		      (append (loop for i from 0 upto (- p 2)
				 collect (unit-vector i dim))
			      (list (make-vector dim :initial-element -1)))
		      nil)))

(defun vector-tensor-product (v w)
  (let ((u (make-vector (* (length v) (length w))))
	(fp 0))
    (loop 
       for a across v
       do (loop for b across w
	     do (setf (aref u fp) (* a b))
	     (incf fp))
       finally (return u))))

#| (vector-tensor-product #(1 2 3) #(4 5)) |#

(defun tensor-product (&rest factors)
  (case (length factors)
    (0 (poly:make-v-poly 0 nil nil))
    (1 (first factors))
    (t (let ((rest-tensor-product (apply #'tensor-product (rest factors)))
	     (first-poly (first factors)))
	 (poly:make-v-poly (* (poly:adim first-poly)
			      (poly:adim rest-tensor-product))
			   (loop for v across (poly:poly-points first-poly)
			      append (loop for w across (poly:poly-points rest-tensor-product)
					collect (vector-tensor-product v w)))
			   nil)))))
			   
(defun cyclotomic-polytope/prime-product (&rest primes)
  (assert (apply #'/= primes))
  (let ((factors (mapcar #'cyclotomic-polytope/prime primes)))
    (apply #'tensor-product factors)))

(in-package :poly)

;; This wants to go into ATROP

(defun print-homogenized-vector (dst-stream homog vec)
  (loop for r of-type rational across vec
     initially (format dst-stream "~&~A " homog)
     do (format dst-stream "~A " r)
     finally (format dst-stream "~%")))

(defun write-poly-latte-vrep-format (p &key (dst *standard-output*))
  "Write out the v-poly P in LattE vrep format.
Keyword argument DST specifies output stream (default: standard-output) or a pathname or string."
  (declare (type v-poly p))
  (with-open-file-or-stream (dst-stream dst :direction :output)
    (format dst-stream "~4D ~4D~%" 
	    (+ (length (poly-points p))
	       (length (poly-rays p)))
	    (1+ (adim p)))
    (loop for point across (poly-points p)
       do (print-homogenized-vector dst-stream 1 point))
    (loop for ray across (poly-rays p)
       do (print-homogenized-vector dst-stream 0 ray))))

(in-package :cl-user)

(defun m-prime (cyclotomic-polytope)
  ;; Add the zero point as an additional generator.
  ;; The lifting into the next dimension is implicit in the vrep format.
  (poly:make-v-poly (poly:adim cyclotomic-polytope)
		    (cons (make-vector (poly:adim cyclotomic-polytope)
				       :initial-element 0)
			  (coerce (poly:poly-points cyclotomic-polytope) 'list))
		    nil))

(defun write-cyclo-mprime-file (&rest primes)
  (with-open-file (f (format nil "/home/mkoeppe/w/latte/EXAMPLES/beck/cyclo-mprime~{-~A~}" primes)
		     :direction :output :if-exists :supersede)
    (poly::write-poly-latte-vrep-format (m-prime (apply #'cyclotomic-polytope/prime-product primes)) 
					:dst f)))  

#|
(write-cyclo-mprime-file 3 5)
(write-cyclo-mprime-file 2 3 5)
(write-cyclo-mprime-file 2 3 7)
(write-cyclo-mprime-file 3 5 7)


|#