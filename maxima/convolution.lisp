; Compute the smallest (integer) coefficient for converting (by multiplication)
; a list of rational numbers to a list of integers; this is the LCM of all
; denominators.
 (defmacro coeff-normalize-list-fractions (v)
   `(labels ((rl (m x)
               (if x (rl (lcm m (denominator (car x))) (cdr x)) m)))
      (rl 1 ,v)))


; Convolution between two series (lists of coefficients); the final size is the
; size of the shortest list
(defun convolution (a b)
  (labels ((main (ar br rev comp)
             (if (and ar br)
               (let ((x (cons (car br) rev)))
                 (main (cdr ar) (cdr br) x
                       (cons (sub a x 0) comp)))
               (nreverse comp)))
           (sub (c d s)
             (if d (sub (cdr c) (cdr d) (+ s (* (car c) (car d)))) s)))
    (main a b NIL NIL)))

; Convolution between one series and one polynomial; the final size is the
; size of the longest list (first argument).
; The polynomial (second argument) MUST have less coefficients than the series.
(defun convolution-poly (a b)
  (labels ((main (ar br rev comp)
             (if ar
               (let* ((br2 (if br br '(0)))
                      (x (cons (car br2) rev)))
                 (main (cdr ar) (cdr br2) x
                       (cons (sub a x 0) comp)))
               (nreverse comp)))
           (sub (c d s)
             (if d (sub (cdr c) (cdr d) (+ s (* (car c) (car d)))) s)))
    (main a b NIL NIL)))

; Compute the reciprocal of a series (list of coefficient); the first coefficient
; MUST not be zero.
(defun convolution-reciprocal (l)
  (labels ((main (a m)
             (if a
               (main (cdr a) (cons (/ (sub (cdr l) m 0) (car l)) m))
               (nreverse m)))
           (sub (v m s)
             (if m (sub (cdr v) (cdr m) (- s (* (car v) (car m)))) s)))
    (main (cdr l) (list (/ 1 (car l))))))
             
(defun  recurrence-vector-raw (v)
  (let ((z (floor (/ (length v) 2))))
    (labels ((main (l q1 q2 sz)
               (if (<= sz z)
                 (multiple-value-bind (qq1 m)
                   (labels ((rl (w z)
                              (if (= (car w) 0)
                                (if (cdr w)
                                  (rl (cdr w) (cons 0 z))
                                  (values NIL NIL))
                                (values z w))))
                     (rl l q1))
                   (if qq1
                     (multiple-value-bind (q s)
                       (labels ((rl (k2 k1 o ss)
                                  (if (or k1 k2)
                                    (let ((kk1 (if k1 k1 '(0)))
                                          (kk2 (if k2 k2 '(0))))
                                      (rl (cdr kk2) (cdr kk1)
                                          (cons
                                            (+ (/ (car kk2) (car m)) (car kk1))
                                            o)
                                          (+ 1 ss)))
                                    (purge o ss)))
                                (purge (o ss)
                                  (if (= 0 (car o))
                                    (purge (cdr o) (- ss 1))
                                    (values (nreverse o) ss))))
                         (rl q2 qq1 NIL 0))
                       (main (cdr (convolution-reciprocal m)) (cons 0 q2) q s))
                     q2))
                 NIL)))
      (main v '(0) '(1) 1))))

(defun recurrence-vector (v)
  (let* ((l (recurrence-vector-raw v))
         (c (coeff-normalize-list-fractions l)))
    (mapcar #'(lambda (a) (* c a)) l)))

(defun ggf (v)
  (let ((l (recurrence-vector-raw v)))
    (if l
      (let* ((s (labels ((purge (x) (if (= 0 (car x)) (purge (cdr x)) x) ))
                 (nreverse (purge (nreverse (convolution-poly v l))))))
             (c (lcm
                  (coeff-normalize-list-fractions l)
                  (coeff-normalize-list-fractions s))))
        (list
          (mapcar #'(lambda (a) (* c a)) s)
          (mapcar #'(lambda (a) (* c a)) l)))
      NIL)))


;;; MAXIMA interface
;;; ================
; Take a Maxima list containing either rational or integer numbers and convert
; it into a lisp list
; TODO: probably not the cleanest way to access Maxima rationals
;       (see http://comments.gmane.org/gmane.comp.mathematics.maxima.general/28954 )
(defmacro from-maxima-list (l)
  `(labels ((c (i)
              (if i
                (cons
                  (if
                    (and (listp (car i))
                         (eq 'rat (caaar i)))
                    (apply '/ (cdar i))
                    (car i))
                  (c (cdr i)))
                NIL)))
     (c (cdr ,l))))

(defmacro to-maxima-list (l)
  `(cons '(mlist simp) ,l))

(defmacro to-maxima-polynomial (v x)
  `(if (cdr ,v)
    (labels ((p (i n)
            (if i
              (cons 
                (if (= 1 (car i))
                  (list '(mexpt simp) x n) 
                  (list '(mtimes simp) (car i)
                            (list '(mexpt simp) (quote $x) n)))
                    (p (cdr i) (+ 1 n)))
              NIL)))
      (cons '(mplus simp) (cons (car ,v) (p (cdr ,v) 1))))
    (car ,v)))

(defmacro to-maxima-ratfrac (f x)
  `(list '(mtimes simp)
    (list '(mexpt simp) (to-maxima-polynomial (cadr ,f) x) -1)
    (to-maxima-polynomial (car ,f) x)))

(defun $recvec (v)
  (to-maxima-list (recurrence-vector-raw (from-maxima-list v))))

(defun $recvecn (v)
  (to-maxima-list (recurrence-vector (from-maxima-list v))))

(defun $ggf (v &optional (x (quote $x)))
  (to-maxima-ratfrac (ggf (from-maxima-list v)) x))

;(%i5) recvec(makelist( 1/2^i, i, 12));
;(%o5)                              [- 2, 1]
;(%i6) recvec(makelist( 1/2^i+1, i, 12));
;(%o6)                             [2, - 3, 1]
;(%i7) 
;(%i6) recvec(makelist( fib(i), i, 12));
;(%o6)                             [- 1, 1, 1]
;(%i7) recvec(makelist( 1/2^(12-i), i, 12));
;(%o7)                             [- 1/2, 1]
;(%i10) recvecn(makelist( 1/2^(12-i), i, 12));
;(%o10)                             [- 1, 2]
;(%i64) ggf(makelist( fib(i), i, 12));
;                                        1
;(%o64)                           - -----------
;                                    2    1
;                                   x  + x  - 1
;(%i65) ggf(makelist( fib(i), i, 12),y);
;                                        1
;(%o65)                           - -----------
;                                    2    1
;                                   y  + y  - 1
