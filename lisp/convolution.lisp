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
         (c (labels ((rl (m x)
                       (if x (rl (lcm m (denominator (car x))) (cdr x)) m)))
              (rl 1 l))))
    (mapcar #'(lambda (a) (* c a)) l)))

(defun ggf (v)
  (let ((l (recurrence-vector v)))
    (if l
      (labels ((rl (x) (if (= 0 (car x)) (rl (cdr x)) x))
               (ql (m x) (if x (ql (lcm m (denominator (car x))) (cdr x)) m)))
        (let ((s (nreverse (rl (nreverse (convolution-poly v l)))))
              (c (ql 1 (nreverse (rl (nreverse (convolution-poly v l)))))))
          (list
            (mapcar #'(lambda (a) (* c a)) s)
            (mapcar #'(lambda (a) (* c a)) l))))
      NIL)))
    ;1/2 5/3 8/3 5/3 -23/6 -175/12 -599/24 -895/48 2713/96 24305/192 88969/384 153185/768 -301703/1536 -3346735/3072
; TODO: multiply s by lcm of denominators
