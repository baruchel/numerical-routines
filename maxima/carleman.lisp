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
;
; Installation
; ------------
; The functions can be used with or without compilation:
;   * without compilation:
;         load("carleman.lisp")$
;   * with compilation (must be compiled only once):
;         :lisp  (compile-file "carleman.lisp");
;     look for the compiled file like "convolution.o" and from now on:
;         load("carleman.o")$

; Compute the Carleman matrix for series whose coefficients are in v.
; Return a list of lists.
; The vector v contains Maxima objects, but the car '(mlist simp) should
; be removed before calling the function.
(defun carleman (v)
  (let ((n (list-length v)))
    (loop repeat n
          with w = (reverse v)
          for u = (cons 1 (make-list (1- n) :initial-element 0))
            then (loop for x on w
                       with s = NIL
                       do (setf s
                            (cons
                              (addn
                                  (loop for a in x
                                        for b in u
                                        collect (mul a b))) s))
                       finally (return s))
          collect u)))


; Formula (4.17) - but there seems to be some misprint in the PDF and
; the sign "-" has been added here
(defun carleman-diag-left (m d)
  (apply #'mapcar #'list
    (loop for w in m
          for i from 1
          for y = (list (cdr w)) then (cons (nthcdr i w) y)
          for d1 = 1 then (mul d1 d)
          for r = (cdar m) then (cdr r) ; exactly the required number of 0's !!!
          collect (loop with z = (list 1)
                        for u in y
                        for x = z
                          then (cons 
                                 (div
                                   (addn
                                       (loop for e in x
                                             for f in u
                                             collect (mul e f)))
                                   (sub d1 d2)) x)
                        for d2 = (div d1 d) then (div d2 d)
                        finally (setf (cdr z) r)
                                (return x)))))


(defun carleman-diag-middle (m d)
  (loop for NIL in m
        for z = (cdar m) then (cdr z)
        for q = NIL then (cons 0 q)
        for d1 = 1 then (mul d1 d)
        collect (append q (cons d1 z))))


; Formula (4.16) in "Continuous time evolution form iterated maps and
; Carleman linearization" (Gralewicz and Kowalski)
(defun carleman-diag-right (m d)
  (loop for NIL in m
        for d1 = 1 then (mul d1 d)
        ; transpose matrix m to z and iterate on rows of z
        for z = (cdr (apply #'mapcar #'list m)) then (mapcar #'cdr (cdr z))
        for q = NIL then (cons 0 q)
        collect (loop with x = (list 1)
                      for y = x then (cdr y)
                      for u in z
                      for d2 = (mul d1 d) then (mul d2 d)
                      do (push (div
                                 (addn
                                     (loop for e in x
                                           for f in u
                                           collect (mul e f)))
                                 (sub d1 d2)) (cdr y))
                      finally (return (append q x)))))


;;; MAXIMA interface
;;; ================

; Return the Carleman matrix of a function whose Taylor expansion is
; given as a list of coefficients.
(defun $carleman (v)
  (simplifya (cons '($matrix)
        (mapcar #'(lambda (x) (simplify (cons '(mlist) x)))
                (carleman (cdr v))))))

; Let M be the Carleman matrix of a function having 0 as a fixed point
; (ie. f(0)=0) and f'(0) not in {0, 1} ; now, V(M) is such
; that M = V^(-1) . L . V with L a diagonal matrix of eigenvalues.
;
; Return the diagonalized Carleman matrix of a function whose Taylor
; expansion is given as a list of coefficients.
; The first coefficient MUST be 0 (since f(0)=0 is a fixed point).
; The second coefficient MUST be some positive value different from 1.
;
; The function returns a list of three matrices whose product is the
; Carleman matrix.
(defun $carleman_diag (v)
  (let ((m (carleman (cdr v)))
        (d (caddr v)))
    (simplifya (list '(mlist)
      ; left part
      (simplifya
        (cons '($matrix) (mapcar #'(lambda (x) (simplify (cons '(mlist) x)))
                              (carleman-diag-left m d))))
      ; diagonal matrix
      (simplifya
        (cons '($matrix) (mapcar #'(lambda (x) (simplify (cons '(mlist) x)))
                              (carleman-diag-middle m d))))
      ; right part
      (simplifya
        (cons '($matrix) (mapcar #'(lambda (x) (simplify (cons '(mlist) x)))
                              (carleman-diag-right m d))))))))
  
