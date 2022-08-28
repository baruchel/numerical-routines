; Copyright (c) 2022 Thomas Baruchel
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

; compute the Carleman matrix for series whose coefficients are in v
; return a list of lists
(defun carleman (v)
  (let ((n (list-length v)))
    (loop repeat n
          with w = (reverse v)
          for u = (cons 1 (make-list (1- n) :initial-element 0))
            then (loop for x on w
                       with s = NIL
                       do (setf s (cons (loop for a in x
                                              for b in u
                                              summing (* a b)) s))
                       finally (return s))
          collect u)))

; compute the Carleman matrix for series whose coefficients are in v
; return an array of lists
(defun carleman-array (v)
  (let* ((n (list-length v)) (a (make-array n)))
    (loop for i below n
          with w = (reverse v)
          for u = (cons 1 (make-list (1- n) :initial-element 0))
            then (loop for x on w
                       with s = NIL
                       do (setf s (cons (loop for a in x
                                              for b in u
                                              summing (* a b)) s))
                       finally (return s))
          do (setf (aref a i) u))
    a))

; compute the Carleman matrix for series whose coefficients are in v
; return an array of arrays
(defun carleman-array2 (v)
  (let* ((n (list-length v)) (a (make-array n))
                             (j (make-array n :initial-element 0)))
    (setf (aref j 0) 1)
    (setf (aref a 0) j)
    (loop for i from 1 below n
          with w = (reverse v)
          do (let ((s (make-array n :initial-element 0)))
               (setf (aref a i) s)
               (loop for x on w
                     for j from (1- n) downto 0
                     do (setf (aref s j) (loop for e in x
                                               for f across (aref a (1- i))
                                               summing (* e f))))))
    a))

(defun carleman-array3 (v)
  (let* ((n (list-length v)) (a (make-array n))
                             (j (make-array n :initial-element 0)))
    (setf (aref j 0) 1)
    (setf (aref a 0) j)
    (loop for i from 1 below n
          for y across a
          with w = (reverse v)
          do (let ((s (make-array n :initial-element 0)))
               (setf (aref a i) s)
               (loop for x on w
                     for j from (1- n) downto 0
                     do (setf (aref s j) (loop for e in x
                                               for f across y
                                               summing (* e f))))))
    a))



; let M be the Carelman matrix of a function having 0 as a fixed point
; (ie. f(0)=0) and f'(0) not in {0, 1} ; now, V(M) is such
; that M = V^(-1) . L . V with L a diagonal matrix of eigenvalues
; Formula (4.16) in "Continuous time evolution form iterated maps and
; Carleman linearization" (Gralewicz and Kowalski)
(defun carleman-diag-right (m d)
  (loop for NIL in m
        for i from 0
        for d1 = 1 then (* d1 d)
        ; transpose matrix m to z and iterate on rows of z
        for z = (cdr (apply #'mapcar #'list m)) then (mapcar #'cdr (cdr z))
        collect (loop with x = (list 1)
                      for y = x then (cdr y)
                      for u in z
                      for d2 = (* d1 d) then (* d2 d)
                      do (setf (cdr y)
                               (list (/ (loop for e in x
                                              for f in u
                                              summing (* e f))
                                        (- d1 d2))))
                      finally (loop repeat i do (setf x (cons 0 x)))
                              (return x))))

; Formula (4.17) - but there seems to be some misprint in the PDF and
; the sign "-" has been added here
; Instead of returning the true matrix:
;
;      (a1 a2 a3 a4 … )
;      (b1 b2 b3 b4 … )
;      (c1 c2 c3 c4 … )
;      (d1 d2 d3 d4 … )
;
; a partial computation of the transpose is returned:
;
;      (a1)
;      (a2 b2)
;      (a3 b3 c3)
;      (a4 b4 c4 … )
;
(defun carleman-diag-left-partial (m d)
  (loop for w in m
        for i from 1
        for y = (list (cdr w)) then (cons (nthcdr i w) y)
        for d1 = 1 then (* d1 d)
        collect (loop
                      for u in y
                      for x = '(1)
                        then (cons 
                               (/ (loop for e in x
                                        for f in u
                                        summing (* e f))
                                  (- d1 d2)) x)
                      for d2 = (/ d1 d) then (/ d2 d)
                      finally (return x))))

; return only row 2 of the left matrix
(defun carleman-diag-left-row2 (m d)
  (cons 0 (mapcar #'cadr (cdr (carleman-diag-left-partial m d)))))

(defun carleman-diag-left (m d)
  (apply #'mapcar #'list
    (loop for w in m
          for i from 1
          for y = (list (cdr w)) then (cons (nthcdr i w) y)
          for d1 = 1 then (* d1 d)
          for r = (cdar m) then (cdr r) ; exactly the required number of 0's !!!
          collect (loop with z = (list 1)
                        for u in y
                        for x = z
                          then (cons 
                                 (/ (loop for e in x
                                          for f in u
                                          summing (* e f))
                                    (- d1 d2)) x)
                        for d2 = (/ d1 d) then (/ d2 d)
                        finally (setf (cdr z) r)
                                (return x)))))


;;; MAXIMA interface
;;; ================
(defun $carleman (v)
  (cons '($matrix simp)
        (mapcar #'(lambda (x) (cons '(mlist simp) x))
                (carleman (cdr v)))))
