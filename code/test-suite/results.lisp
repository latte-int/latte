(require :aserve)
(use-package :net.html.generator)

(defun parse-summary-file (filename)
  (with-open-file (f filename) 
    (loop 
       with alist = '()
       for line = (read-line f nil :eof)
       until (eql line :eof)
       do (multiple-value-bind (ok whole-match name time-string)
	      (excl:match-regexp "\\([^:]*\\):.*GOOD Time: \\(.*\\) sec"
				 line :return :string)
	    (declare (ignore whole-match))
	    (when ok
	      (let ((time (read-from-string time-string)))
		(push (cons name time) alist))))
       finally (return (nreverse alist)))))

(defun make-results-table-html (stream result-file-names)
  (let* ((result-nicknames (loop for result in result-file-names
			     for nick across "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
			      collect nick))
	 (result-alists (mapcar #'parse-summary-file result-file-names))
	 (examples (reduce #'union result-alists 
			   :key (lambda (alist) (mapcar #'car alist)))))
    (html-stream stream
		 (loop for nick in result-nicknames
		    for file-name in result-file-names
		    do (html (:p)
			     (:princ nick) 
			     (:princ " ")
			     (:princ file-name)))
		 ((:table border 1)
		  (:tr
		   (:td "Example")
		   (loop for nick in result-nicknames
		      do (html (:td (:princ nick)))))
		  (loop for example in examples
		     do (html (:tr 
			       (:td (:princ example)))))))))

;;(parse-summary-file "/home/mkoeppe/w/latte/results/log-2006-04-01-count --exp --maxdet=10-pid2244@zeta/summary")
  

  
		
(with-open-file (stream "benchmark.html" :direction :output
			:if-exists :supersede)
  (html-stream stream 
	       (:html
		(:head (:title "LattE macchiato benchmarks"))
		(:body
		 (:h1 "LattE macchiato benchmarks")
		 (make-results-table-html stream
					  '("/home/mkoeppe/w/latte/results/log-2006-04-01-count --exp --maxdet=10-pid2244@zeta/summary")
					  )))))

			   