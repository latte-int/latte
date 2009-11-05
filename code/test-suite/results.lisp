(require :aserve)
(use-package :net.html.generator)

(defun parse-summary-file (filename)
  (with-open-file (f filename) 
    (loop 
       with alist = '()
       for line = (read-line f nil :eof)
       until (eql line :eof)
       do (multiple-value-bind (ok whole-match name time-string)
	      (excl:match-regexp "\\([^:]*\\):.*Result:.*Time: \\(.*\\) sec"
				 line :return :string)
	    (declare (ignore whole-match))
	    (when ok
	      (let ((time (read-from-string time-string)))
		(push (cons name time) alist))))
       finally (return (nreverse alist)))))

(defun read-dimensions (latte-filename)
  (with-open-file (f latte-filename) 
    (values (read f) (read f))))

(defun make-results-table-html (stream result-file-names)
  (let* ((result-nicknames (loop for result in result-file-names
			     for nick across "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
			      collect nick))
	 (result-alists (mapcar #'parse-summary-file result-file-names))
	 (examples (sort (reduce (lambda (a b) (union a b :test #'string=))
				 result-alists 
				 :key (lambda (alist) (mapcar #'car alist)))
			 #'string<=)))
    (html-stream stream
		 (:p "Key to test runs:")
		 (loop for nick in result-nicknames
		    for file-name in result-file-names
		    do (html (:p)
			     (:princ nick) 
			     (:princ " ")
			     (:tt (:princ file-name))))
		 (:p)
		 ((:table border 1)
		  ((:tr bgcolor "blue")
		   (:td "Example")
		   (:td (:i "m"))
		   (:td (:i "n"))
		   (loop for nick in result-nicknames
		      do (html ((:td align "center") (:princ nick)))))
		  (loop for example in examples
		     do (html (:tr 
			       (:td (:tt (:princ example)))
			       (multiple-value-bind (m n)
				   (read-dimensions (concatenate 'string 
								 *examples-directory*
								 example))
				 (html ((:td align "right") (:princ m))
				       ((:td align "right") (:princ n))))
			       (let* ((times 
				       (loop for alist in result-alists
					  collect (cdr (assoc example alist :test #'string=))))
				      (min-time 
				       (loop for time in times
					  if time
					  minimize time)))
				 (loop for time in times
				    do (html ((:td align "right"
						   :if* (and time (= min-time time))
						   bgcolor 
						   "yellow")
					      (if time (html (:princ (format nil "~,2F" time)))
						  (html "&nbsp;")))))))))))))

;;(parse-summary-file "/home/mkoeppe/w/latte/results/log-2006-04-01-count --exp --maxdet=10-pid2244@zeta/summary")
  

(defparameter *examples-directory* "/home/mkoeppe/w/latte/EXAMPLES/")

(with-open-file (stream "/home/mkoeppe/public_html/latte/latte-benchmark-2007-11-01.html" :direction :output
			:if-exists :supersede)
  (html-stream stream 
	       (:html
		(:head (:title "LattE macchiato benchmarks"))
		(:body
		 (:h1 "LattE macchiato benchmarks")
		 (:p "We provide benchmarks comparing old LattE (version 1.2; shown in column A) and "
		     ((:a href "http://www.math.uni-magdeburg.de/~mkoeppe/latte/")
		      "LattE macchiato")
		     " (version 1.2-mk-0.9; column B; various \"flavors\" shown in the remaining columns).  The running time is given in CPU seconds on an Intel Core Duo running at 2 GHz.")
		 (make-results-table-html stream
					  '("/home/mkoeppe/w/latte/results/log-2007-10-31-latte-1.1-count -pid24980@moose/summary"
					    "/home/mkoeppe/w/latte/results/log-2007-10-31-count -pid5822@moose/summary"
					    "/home/mkoeppe/w/latte/results/log-2007-11-01-count --irr --exp --maxdet=1000-pid7391@moose/summary"
					    "/home/mkoeppe/w/latte/results/log-2007-11-01-count --all-primal --exp --maxdet=1000-pid5352@moose/summary"
					    "/home/mkoeppe/w/latte/results/log-2007-11-01-count --triangulation=4ti2 --dualization=4ti2 --compute-vertex-cones=4ti2-pid8371@moose/summary"
					    "/home/mkoeppe/w/latte/results/log-2007-11-01-count --irr --exp --maxdet=1000 --triangulation=4ti2 --dualization=4ti2 --compute-vertex-cones=4ti2-pid10754@moose/summary"
					    "/home/mkoeppe/w/latte/results/log-2007-11-01-count --all-primal --exp --maxdet=1000 --triangulation=4ti2 --dualization=4ti2 --compute-vertex-cones=4ti2-pid11951@moose/summary"
					    ))))))

			   