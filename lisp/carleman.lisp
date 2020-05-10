; Identity Matrix (as an array NOT AS lists)
(defun eye (n)
  (loop
    for i below n
    with m = (make-array (list n n) :initial-element 0)
    do (setf (aref m i i) 1)
    finally (return m)))

; Convolution of two power series
(defun convolution (a b)
  (loop
    for NIL in a
    for y in b
    for z = (list (car b)) then (cons y z)
    collect (loop
              for i in a
              for j in z
              sum (* i j))))

; Compute iterates of some power series
; Argument must be a list following some rules:
;   * initial coefficient MUST BE 0;
;   * second coefficient is a positive value != 1
; Return a Lambda function
; The lambda is used with a new value for the derivative at 0
; (new value must be a positive value != 1); returns the new coefficients
;  eg. (funcall (C '(0 2 0 2 0 2 0 2)) 3)
(defun C (v)
  (let*
    ((n (length v))
     (l (cadr v))
     (c (make-array (list n n)
                    :initial-contents 
                      (loop
                        for NIL in v
                        for y = (cons 1 (make-list (1- (length v))
                                                   :initial-element 0))
                                then (convolution y v)
                        collect y)))
     (*vi (loop
           with vi* = (eye n)
           for k from 1 below n
           do (loop
                for j from (1- k) downto 1
                do (setf (aref vi* j k)
                         (- (/ (reduce #'+
                                    (loop for u from (1+ j) to k
                                          collect (* (aref vi* u k)
                                                     (aref c j u))))
                            (- (expt l j) (expt l k))))))
           finally (return
                     (loop for i below n collect (aref vi* 1 i)))))
     (*v (loop
           with v* = (eye n)
           for j below (1- n)
           do (loop
                for k from (1+ j) below n
                do (setf (aref v* j k)
                         (/ (reduce #'+
                                    (loop for u from j below k
                                          collect (* (aref v* j u)
                                                     (aref c u k))))
                            (- (expt l j) (expt l k)))))
           finally (return v*))))
    (lambda (d2)
      (let ((d2* (loop
                   for i below n
                   for j = 1 then (* j d2)
                   collect j)))
        (loop
          for i below n
          collect
            (loop for j to i
                  for d in d2*
                  for k in *vi
                  sum (* k d (aref *v j i))))))))
