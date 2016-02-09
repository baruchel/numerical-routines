; Copyright (c) 2016 Thomas Baruchel
; 
; Permission is hereby granted, free of charge, to any person obtaining a copy
; of this software and associated documentation files (the "Software"), to deal
; in the Software without restriction, including without limitation the rights
; to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
; copies of the Software, and to permit persons to whom the Software is
; furnished to do so, subject to the following conditions:
; 
; The above copyright notice and this permission notice shall be included in
; all copies or substantial portions of the Software.
; 
; THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
; IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
; FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
; AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
; LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
; OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
; SOFTWARE.

; Compute the smallest (integer) coefficient for converting (by multiplication)
; a list of rational numbers to a list of integers; this is the LCM of all
; denominators.
(defmacro coeff-normalize-list-fractions (v)
  `(reduce #'lcm (mapcar #'denominator ,v)))

; Convolution between two series (lists of coefficients); the final size is the
; size of the shortest list
(defun convolution (a b)
  (loop
    for NIL in a
    for y in b
    for z = (list (car b)) then (cons y z)
    collect (loop
              for i in a
              for j in z
              sum (* i j))))
    

; Convolution between one series and one polynomial; the final size is the
; size of the longest list (first argument).
; The polynomial (second argument) MUST have less coefficients than the series.
(defun convolution-poly (a b)
  (loop
    for NIL in a
    for y = b then (if (cdr y) (cdr y) '(0))
    for z = (list (car b)) then (cons (car y) z)
    collect (loop
              for i in a
              for j in z
              sum (* i j))))

; Compute the reciprocal of a series (list of coefficient); the first coefficient
; MUST not be zero.
(defun convolution-reciprocal (l)
  (loop
    for NIL in l
    for m = (list (/ 1 (car l))) then
      (cons (/ (-
        (loop
          for i in (cdr l)
          for j in m
          sum (* i j)))
        (car l)) m)
    finally (return (nreverse m))))
             
; Compute the minimal recurrence vector; returned coefficients are rational
; (though integer coefficients may often been returned, non integer ones can
; also been returned in some cases).
; Warning: an empty sequence returns (1)
(defun  recurrence-vector-raw (v)
  (loop named main
    with z = (floor (/ (length v) 2))
    with l = v
    with q1 = '(0)
    for q2 = '(1)
           then (multiple-value-bind (m qq1)
                  (loop
                    for a on l
                    for b = q1 then (cons 0 b)
                    finally (return (values NIL NIL))
                    do (if (/= 0 (car a)) (return (values a b))))
                  (if qq1
                    (loop named add-and-compute-size
                      for k2 = q2 then (if (cdr k2) (cdr k2) '(0))
                      for k1 = qq1 then (if (cdr k1) (cdr k1) '(0))
                      for o = (cons (+ (/ (car q2) (car m)) (car qq1)) o)
                            then (cons (+ (/ (car k2) (car m)) (car k1)) o)
                      for ss from 1
                      do (if (not (or (cdr k1) (cdr k2)))
                           (loop
                             for oo on o
                             for i downfrom ss
                             do (if (/= 0 (car oo))
                                  (if (> i z)
                                    (return-from main NIL)
                                    (progn
                                      (setf q1 (cons 0 q2))
                                      (setf l (cdr (convolution-reciprocal m)))
                                      (return-from add-and-compute-size
                                        (nreverse oo))))))))
                    (return-from main q2)))))

; Compute the minimal recurrence vector; returned coefficients are integers.
(defun recurrence-vector (v)
  (let* ((l (recurrence-vector-raw v))
         (c (coeff-normalize-list-fractions l)))
    (mapcar #'(lambda (a) (* c a)) l)))

; Compute the generating function for a sequence if g.f is a rational function.
(defun ggf (v)
  (let ((l (recurrence-vector-raw v)))
    (if l
      (let* ((s (nreverse
                  (loop
                    for x on (nreverse (convolution-poly v l))
                    do (if (/= 0 (car x)) (return x)))))
             (c (lcm
                  (coeff-normalize-list-fractions l)
                  (coeff-normalize-list-fractions s))))
        (list
          (mapcar #'(lambda (a) (* c a)) s)
          (mapcar #'(lambda (a) (* c a)) l)))
      NIL)))
