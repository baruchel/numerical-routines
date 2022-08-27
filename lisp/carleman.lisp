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



; (defun carleman-diag-right (m)
;   (let* ((n (list-length m))
;          (v (make-array n))
;          (r (make-array n :initial-element 0)))
;     (setf (aref r 0) 1)
;     (setf (aref v 0) r)
;     (loop for j from 1 below (1- n)
;           do (let ((s (make-array n :initial-element 0)))
;                (setf (aref s j) 1)
;                (setf (aref v j) s)
;                (loop for k from (1+ j) below n
;                      do (setf (aref s k)
;                               (loop for i from j below k
;                                     sum (* (aref s i)
;                                            (aref (aref
; 
; 
;
(defun carleman-diag (v)
  (let ((m (carleman v)))
    ; let M be the Carelman matrix of a function having 0 as a fixed point
    ; (ie. f(0)=0) and f'(0) not in {0, 1} ; now, V(M) is such
    ; that M = V^(-1) . L . V with L a diagonal matrix of eigenvalues
    ; Formula (4.16) in "Continuous time evolution form iterated maps and
    ; Carleman linearization" (Gralewicz and Kowalski)
    (loop for NIL in v
          for i from 0
          with d = (cadr v)
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
                                (return x)))))

; Formula (4.17) - but there seems to be some misprint in the PDF and
; the sign "-" has been added here; furthermore, since loops in pari-gp
; only work with increasing values, 'j' has been replaced with 'k-j'
; in order to make the (original) j go from k-1 down to 1.
; The case (original) j = 0 has been discarded as useless.
(defun carleman-diag2 (v)
  (let ((m (carleman v)))
    ; let M be the Carelman matrix of a function having 0 as a fixed point
    ; (ie. f(0)=0) and f'(0) not in {0, 1} ; now, V(M) is such
    ; that M = V^(-1) . L . V with L a diagonal matrix of eigenvalues
    ; Formula (4.16) in "Continuous time evolution form iterated maps and
    ; Carleman linearization" (Gralewicz and Kowalski)
    (loop for NIL in v
          for i from 0
          with w = (reverse m)
          with d = (cadr v)
          for d1 = 1 then (* d1 d)
          collect (loop
                        for x = (list 1)
                          then (cons 
                                 (/ (loop for e in x
                                          for f in (nthcdr j u)
                                          summing (* e f))
                                    (- d1 d2)) x)
                        for u in (nthcdr (- (list-length v) i) w)
                        for j from i downto 0
                        for d2 = (/ d1 d) then (/ d2 d)
                        finally (return x)))))
(defun carleman-diag3 (v)
  (let ((m (carleman v)))
    ; let M be the Carelman matrix of a function having 0 as a fixed point
    ; (ie. f(0)=0) and f'(0) not in {0, 1} ; now, V(M) is such
    ; that M = V^(-1) . L . V with L a diagonal matrix of eigenvalues
    ; Formula (4.16) in "Continuous time evolution form iterated maps and
    ; Carleman linearization" (Gralewicz and Kowalski)
    (loop for w in m
          for i from 0
          for y = (list w) then (cons w y)
          with d = (cadr v)
          for d1 = 1 then (* d1 d)
          collect (loop
                        for u in y
                        for x = (list 1)
                          then (cons 
                                 (/ (loop for e in x
                                          for f in (nthcdr j u)
                                          summing (* e f))
                                    (- d1 d2)) x)
                        for j from i downto 0
                        for d2 = (/ d1 d) then (/ d2 d)
                        finally (return x)))))
(defun carleman-diag4 (v)
  (let ((m (carleman v)))
    ; let M be the Carelman matrix of a function having 0 as a fixed point
    ; (ie. f(0)=0) and f'(0) not in {0, 1} ; now, V(M) is such
    ; that M = V^(-1) . L . V with L a diagonal matrix of eigenvalues
    ; Formula (4.16) in "Continuous time evolution form iterated maps and
    ; Carleman linearization" (Gralewicz and Kowalski)
    (loop for w in m
          for i from 1
          for y = (list (cdr w)) then (cons (nthcdr i w) y)
          with d = (cadr v)
          for d1 = 1 then (* d1 d)
          collect (loop
                        for u in y
                        for x = (list 1)
                          then (cons 
                                 (/ (loop for e in x
                                          for f in u
                                          summing (* e f))
                                    (- d1 d2)) x)
                        for d2 = (/ d1 d) then (/ d2 d)
                        finally (return x)))))


;;; MAXIMA interface
;;; ================
(defun $carleman (v)
  (cons '($matrix simp)
        (mapcar #'(lambda (x) (cons '(mlist simp) x))
                (carleman (cdr v)))))
