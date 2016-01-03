(setf l '(1 1 1 1 1 1 1 1))

(defun ratio-add (a b)
  (let* ((q (lcm (cadr a) (cadr b)))
         (c (/ q (cadr a)))
         (d (/ q (cadr b)))
         (p (+ (* c (car a)) (* d (car b))))
         (g (gcd p q)))
    (list (/ p g) (/ q g))))

(defun ratio-mul (a b)
  (let* ((p (* (car a) (car b)))
         (q (* (cadr a) (cadr b)))
         (g (gcd p q)))
    (list (/ p g) (/ q g))))


; (defun convolution (a b)
;   (labels ((main (ar br rev comp)
;              (if (and ar br)
;                (let ((x (cons (car br) rev)))
;                  (main (cdr ar) (cdr br) x
;                        (cons (sub a x 0) comp)))
;                (nreverse comp)))
;            (sub (c d s)
;              (if d ; d is always shortest than c
;                (sub (cdr c) (cdr d)
;                     (+ s (* (car c) (car d))))
;                s)))
;     (main a b NIL NIL)))

(defmacro make-convolution (a b p m z)
  `(labels ((main (ar br rev comp)
             (if (and ar br)
               (let ((x (cons (car br) rev)))
                 (main (cdr ar) (cdr br) x
                       (cons (sub ,a x ,z) comp)))
               (nreverse comp)))
           (sub (c d s)
             (if d ; d is always shorter than c
               (sub (cdr c) (cdr d)
                    (funcall #',p s (funcall #',m (car c) (car d))))
               s)))
    (main ,a ,b NIL NIL)))

(defun c1 (a b) (make-convolution a b + * 0))
(defun c2 (a b) (make-convolution a b ratio-add ratio-mul '(0 1)))
