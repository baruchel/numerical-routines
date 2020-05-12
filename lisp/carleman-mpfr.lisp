(require :sb-mpfr)
(sb-mpfr:set-precision 256) ; 256 bits




; Identity Matrix (as an array NOT AS lists)
(defun eye (n)
  (loop
    for i below n
    with m = (make-array (list n n) :initial-element
                                    (sb-mpfr:coerce 0 'sb-mpfr:mpfr-float))
    do (setf (aref m i i) (sb-mpfr:coerce 1 'sb-mpfr:mpfr-float))
    finally (return m)))

; Compute iterates of some power series
; Argument must be a list following some rules:
;   * list has at least 2 elements;
;   * initial coefficient MUST BE 0;
;   * second coefficient is a positive value != 1
; Return a Lambda function
; The lambda is used with a new value for the derivative at 0
; (new value must be a positive value != 1); returns the new coefficients
;  eg. (funcall (C '(0 2 0 2 0 2 0 2)) 3)
(defun C (v)
  (let*
    ((n (length v))
     (c (loop
          with m = (make-array (list n n) :initial-element
                                    (sb-mpfr:coerce 0 'sb-mpfr:mpfr-float))
          initially (progn
                      (setf (aref m 0 0) (sb-mpfr:coerce 1 'sb-mpfr:mpfr-float))
                      (loop
                        for i below n
                        for j in v
                        do (setf (aref m 1 i) j)))
          for i from 2 below n
          do (loop
               for j from i below n
               with w = (1- i)
               do (setf (aref m i j) (loop
                                       for u from j downto w
                                       for k in v
                                       with s = (sb-mpfr:coerce 0 'sb-mpfr:mpfr-float)
                                       do (setq s
                                            (sb-mpfr:add s
                                              (sb-mpfr:mul (aref m w u) k)))
                                       finally (return s))))
          finally (return m)))
     (*vi (loop
           with vi* = (eye n)
           for k from 1 below n
           do (loop
                for j from (1- k) downto 1
                do (setf (aref vi* j k)
                         (sb-mpfr:negate (sb-mpfr:div
                                        (loop for u from (1+ j) to k
                                              with s = (sb-mpfr:coerce 0 'sb-mpfr:mpfr-float)
                                              do (setq s
                                                   (sb-mpfr:add s
                                                     (sb-mpfr:mul (aref vi* u k)
                                                                  (aref c j u))))
                                              finally (return s))
                            (sb-mpfr:sub (aref c j j) (aref c k k))))))
           finally (return
                     (loop for i below n collect (aref vi* 1 i)))))
     (*v (loop
           with v* = (eye n)
           for j below (1- n)
           do (loop
                for k from (1+ j) below n
                do (setf (aref v* j k)
                         (sb-mpfr:div (loop for u from j below k
                                            with s = (sb-mpfr:coerce 0 'sb-mpfr:mpfr-float)
                                            do (setq s
                                                 (sb-mpfr:add s
                                                   (sb-mpfr:mul (aref v* j u)
                                                                (aref c u k))))
                                            finally (return s))
                            (sb-mpfr:sub (aref c j j) (aref c k k)))))
           finally (return v*))))
    (lambda (d2)
      (let ((d2* (loop
                   for j = (sb-mpfr:coerce 1 'sb-mpfr:mpfr-float)
                           then (sb-mpfr:mul j d2)
                   for k in *vi
                   collect (sb-mpfr:mul k j))))
        (loop
          for i below n
          collect
            (loop for j to i
                  for d in d2*
                  with s = (sb-mpfr:coerce 0 'sb-mpfr:mpfr-float)
                  do (setq s
                       (sb-mpfr:add s (sb-mpfr:mul d (aref *v j i))))
                  finally (return s)))))))
