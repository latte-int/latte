(in-package :cl-user)

(defun make-random-simplex (dimension digits seed)
  (with-open-file (s (format nil "random-simplex-dim~D-digits~D-~D.vrep"
			     dimension digits seed)
		     :direction :output
		     :if-exists :supersede)
    (let ((*random-state* (excl::make-random-state-from-seed seed))
	  (coeff-bound (expt 10 digits)))
      (format s "~D ~D~%" (1+ dimension) (1+ dimension))
      (dotimes (i (1+ dimension))
	(format s "1 ")
	(dotimes (j dimension)
	  (format s "~D " (random coeff-bound)))
	(format s "~%")))))
	
(defun make-random-simplices (dimension digits)
  (loop for seed from 1 to 10 
     do (make-random-simplex dimension digits seed)))


#||
(make-random-simplices 4 2)

(make-random-simplices 3 2)
||#