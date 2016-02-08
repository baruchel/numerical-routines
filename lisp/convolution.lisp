; Compute the smallest (integer) coefficient for converting (by multiplication)
; a list of rational numbers to a list of integers; this is the LCM of all
; denominators.
(defmacro coeff-normalize-list-fractions (v)
  `(reduce #'lcm (mapcar #'denominator ,v)))

; Convolution between two series (lists of coefficients); the final size is the
; size of the shortest list
(defun convolution (a b)
  (labels ((main (ar br rev comp)
             (if (and ar br)
               (let ((x (cons (car br) rev)))
                 (main (cdr ar) (cdr br) x
                       (cons (loop
                               for i in a
                               for j in x
                               sum (* i j)) comp)))
               (nreverse comp))))
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
                       (cons (loop
                               for i in a
                               for j in x
                               sum (* i j)) comp)))
               (nreverse comp))))
    (main a b NIL NIL)))

; Compute the reciprocal of a series (list of coefficient); the first coefficient
; MUST not be zero.
(defun convolution-reciprocal (l)
  (labels ((main (a m)
             (if a
               (main (cdr a) (cons (/ (-
                                        (loop
                                          for i in (cdr l)
                                          for j in m
                                          sum (* i j)))
                                      (car l)) m))
               (nreverse m))))
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
