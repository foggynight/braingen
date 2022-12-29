(import (chicken format)
        (chicken io)
        (chicken process)
        (chicken random)
        (chicken string)
        vector-lib)

;; config ----------------------------------------------------------------------

(define-constant INIT-POPULATION-SIZE 4)
(define-constant INIT-MAX-GENE-LEN 8)
(define-constant TEMP-FILE-PATH "/tmp/braingen.bf")

;; misc ------------------------------------------------------------------------

(define (1+ n) (+ n 1))

(define string->vector (compose list->vector string->list))
(define vector->string (compose list->string vector->list))

;; brainfuck -------------------------------------------------------------------

;; TODO: Write to output port instead of using `printf`.
(define (brainfuck-eval path inputs)
  (define input
    (apply string-append (map (lambda (n)
                                (string-append "\\x" (number->string n 16)))
                              inputs)))
  (define-values (ip op pid ep)
    (process* (sprintf "printf '~A' | brainfuck ~A" input path)))
  (let-values (((_ _ status) (process-wait pid)))
    (let ((out (if (zero? status)
                   (string->number (read-line ip) 16)
                   #f)))
      (close-input-port ip)
      (close-output-port op)
      (close-input-port ep)
      out)))

;; data ------------------------------------------------------------------------

(define (read-data path)
  (map (lambda (l) (map string->number (string-split l)))
       (with-input-from-file path read-lines)))

;; program ---------------------------------------------------------------------

(define-constant BF-CHARS #(#\> #\< #\+ #\- #\[ #\] #\. #\,))
(define-constant BF-CHAR-CNT 8)

(define-record-type program
  (%make-program genes score)
  program?
  (genes program-genes)
  (score program-score program-score-set!))

(define (make-program genes #!optional score)
  (%make-program genes score))

(define-record-printer (program prog op)
  (fprintf op "#<program ~A ~A>"
           (program-score prog)
           (vector->list (program-genes prog))))

(define (random-bf-char)
  (vector-ref BF-CHARS (pseudo-random-integer BF-CHAR-CNT)))

(define (random-program gene-len)
  (define genes (vector-unfold (lambda (i) (random-bf-char)) gene-len))
  (make-program genes))

(define (score-program! program data)
  (define str (vector->string (program-genes program)))
  (with-output-to-file TEMP-FILE-PATH (lambda () (print str)))
  (define score 1) ; 1 is minimum score.
  (let loop ((ds data))
    (unless (null? ds)
      (let* ((d (car ds))
             (result (brainfuck-eval TEMP-FILE-PATH (cdr d))))
        (when result
          (when (= result (car d))
            (set! score (1+ score)))
          (loop (cdr ds))))))
  (program-score-set! program score))

;;(define (program-combine prog1 prog2)
;;  )

;;(define (program-mutate prog)
;;  )

;; population ------------------------------------------------------------------

(define (random-population pop-size max-gene-len)
  (vector-unfold (lambda (i) (random-program
                              (pseudo-random-integer (1+ max-gene-len))))
                 pop-size))

;;(define (select-parents pop)
;;  )

;; main ------------------------------------------------------------------------

(define (main data-path #!optional seed)
  (define data (read-data data-path))
  (set-pseudo-random-seed! (if seed seed (random-bytes)))
  (define population (random-population INIT-POPULATION-SIZE
                                        INIT-MAX-GENE-LEN))
  (vector-for-each (lambda (i p) (score-program! p data)) population)
  (print population)
  )
