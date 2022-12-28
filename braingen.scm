(import (chicken format)
        (chicken io)
        (chicken process)
        (chicken string)
        (chicken random)
        vector-lib)

;; config ----------------------------------------------------------------------

(define-constant INIT-POPULATION-SIZE 4)
(define-constant TEMP-FILE-PATH "/tmp/braingen.bf")

;; misc ------------------------------------------------------------------------

(define (1+ n) (+ n 1))

(define string->vector (compose list->vector string->list))
(define vector->string (compose list->string vector->list))

;; brainfuck -------------------------------------------------------------------

(define (brainfuck-eval path inputs)
  (define input-str (apply string-append
                           (map (lambda (n)
                                  (string-append "\\x" (number->string n 16)))
                                inputs)))
  (define-values (ip op pid)
    (process (sprintf "printf '~A' | brainfuck ~A" input-str path)))
  (let-values (((pid2 _ status) (process-wait pid)))
    (let ((out (if (zero? status)
                   (string->number (read-line ip) 16)
                   #f)))
      (close-input-port ip)
      (close-output-port op)
      out)))

;; candidate -------------------------------------------------------------------

(define-constant BF-CHARS #(#\> #\< #\+ #\- #\. #\, #\[ #\]))
(define-constant BF-CHAR-CNT 8)

(define-record candidate genes)

(define-record-printer (candidate cand op)
  (fprintf op "#<candidate ~A>" (vector->list (candidate-genes cand))))

(define (random-bf-char)
  (vector-ref BF-CHARS (pseudo-random-integer BF-CHAR-CNT)))

(define (random-candidate gene-len)
  (define genes (vector-unfold (lambda (i) (random-bf-char)) gene-len))
  (make-candidate genes))

(define (candidate-score cand dats)
  (define str (vector->string (candidate-genes cand)))
  (with-output-to-file TEMP-FILE-PATH (lambda () (print str)))
  (define score 1) ; 1 is minimum score.
  (let loop ((ds dats))
    (unless (null? ds)
      (let* ((d (car ds))
             (result (brainfuck-eval TEMP-FILE-PATH (cdr d))))
        (cond ((not result))
              ((= result (car d)) (set! score (1+ score)) (loop (cdr ds)))
              (else (loop (cdr ds)))))))
  score)

;;(define (candidate-combine cand1 cand2)
;;  )

;;(define (candidate-mutate cand)
;;  )

;; population ------------------------------------------------------------------

(define (random-population pop-size max-gene-len)
  (vector-unfold (lambda (i) (random-candidate
                              (pseudo-random-integer (1+ max-gene-len))))
                 pop-size))

;;(define (select-parents pop)
;;  )

;; main ------------------------------------------------------------------------

(define (main dats-path #!optional seed)
  (define dats (map (lambda (l) (map string->number (string-split l)))
                    (with-input-from-file dats-path read-lines)))
  (set-pseudo-random-seed! (if seed seed (random-bytes)))
  (define population (random-population 4 8))
  (print dats)
  )
