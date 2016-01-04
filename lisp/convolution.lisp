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

; Convolution between two series (lists of coefficients); the final size is the
; size of the longest list (the shortest list is assumed to represent a finite
; sum of coefficients rather than an infinite series).
; The shortest list MUST be the second one.
(defun convolution-full (a b)
  (labels ((main (ar br rev comp)
             (if (or ar br)
               (let* ((ar2 (if ar ar '(0)))
                      (br2 (if br br '(0)))
                      (x (cons (car br2) rev)))
                 (main (cdr ar2) (cdr br2) x
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
  (let* ((l (recurrence-vector v))
         (s (labels ((rl x) (if (= 0 (car x)) (cdr x) x))
              (nreverse (rl (nreverse (convolution-full v l)))))))
    (list s l)))
