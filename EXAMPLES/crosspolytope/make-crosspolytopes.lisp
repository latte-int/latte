;; Common Lisp program to create cross polytope input files 
;; (V-representation)

(defun make-cross-polytope (dimension)
  (with-open-file (s (format nil "cross-polytope-~D.vrep"
			     dimension)
		     :direction :output
		     :if-exists :supersede)
    (format s "~D ~D~%"
	    (* 2 dimension) (1+ dimension))
    (dotimes (i dimension)
      (format s "1 ")
      (dotimes (j dimension) 
	(format s "~D " (if (= i j) 1 0)))
      (format s "~%")
      (format s "1 ")
      (dotimes (j dimension) 
	(format s "~D " (if (= i j) -1 0)))
      (format s "~%"))))

(loop for dimension from 2 to 30
   do (make-cross-polytope dimension))


