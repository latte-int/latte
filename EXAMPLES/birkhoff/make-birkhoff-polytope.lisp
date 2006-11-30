#||
 (asdf:oos 'asdf:load-op :gywopt)
||#

(defpackage :birkhoff
  (:use :cl :gywopt :gyw :misc :app))

(in-package :birkhoff)

(defun create-birkhoff-problem (dimension)
  (call-with-modeled-problem (GYW-MIN)
    (lambda (p)
      (permutations-model p dimension))))

(defun make-birkhoff-file (dimension)
  (let ((p (create-birkhoff-problem dimension)))
    (write-prob/latte p (format nil "birkhoff-~D" dimension))
    p))


