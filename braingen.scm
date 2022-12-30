(import (chicken format)
        (chicken io)
        (chicken process)
        (chicken process signal)
        (chicken process-context)
        (chicken random)
        (chicken sort)
        (chicken string)
        (chicken time)
        simple-exceptions
        vector-lib)

;; misc ------------------------------------------------------------------------

(define (1+ n) (+ n 1))
(define (1- n) (- n 1))
(define (2+ n) (+ n 2))
(define (2* n) (* n 2))
(define // (compose inexact->exact floor /))

;; string ----------------------------------------------------------------------

(define (string-insert str i c)
  (define left (substring str 0 i))
  (define right (substring str i))
  (string-append left (string c) right))

(define (string-delete str i)
  (define left (substring str 0 i))
  (define right (substring str (1+ i)))
  (string-append left right))

(define (string-swap! str i j)
  (define c (string-ref str i))
  (string-set! str i (string-ref str j))
  (string-set! str j c))

;; config ----------------------------------------------------------------------

(define-constant POPULATION-SIZE 100)        ; Size of the population, remains constant.
(define-constant INIT-MAX-GENE-LEN 80)       ; Max size of the genes in the initial population.
(define PARENT-COUNT (// POPULATION-SIZE 2)) ; Number of parents selected each generation.
(define CHILD-COUNT PARENT-COUNT)            ; Number of children spawned each generation.
(define-constant MAX-GENERATIONS -1)         ; Maximum number of generations, -1 for infinite.

(define-constant MAX-RUNTIME-SECS 2)         ; Max seconds Brainfuck programs are allowed to run.

;; Brainfuck programs are written to/read from this file.
(define-constant TEMP-FILE-PATH "/tmp/braingen.bf")

;; brainfuck -------------------------------------------------------------------

;; TODO: Figure out a way to terminate without using `with-exn-handler`.
;; TODO: Write to output port instead of using `printf`.
;; Evaluate the Brainfuck program contained within the file at `path` given a
;; list of numbers as `inputs`. Each number is converted to a character and
;; written to the input stream of the program. If evaluation was successful,
;; returns the value of the accumulator as a number, otherwise returns false.
(define (brainfuck-eval path inputs genes)
  (printf "braingen: evaluating: ~A\r" inputs)
  (with-exn-handler
   (lambda (exn) (print (message exn)) #f)
   (lambda ()
     (define input
       (apply string-append (map (lambda (n)
                                   (string-append "\\x" (number->string n 16)))
                                 inputs)))
     (define start-time (current-seconds))
     (define-values (ip op pid ep)
       (process* (sprintf "printf '~A' | brainfuck ~A" input path)))
     (define out
       (let loop ((kill #f)) ; kill when over time
         (if kill
             (begin (printf "braingen: infinite loop detected: ~A\n" genes)
                    (process-signal pid)
                    #f)
             (let-values (((pid _ status) (process-wait pid #t)))
               (if (zero? pid) ; process has not terminated
                   (loop (> (- (current-seconds) start-time)
                            MAX-RUNTIME-SECS))
                   (if (zero? status)
                       (string->number (read-line ip) 16)
                       #f))))))
     (close-input-port ip)
     (close-output-port op)
     (close-input-port ep)
     (lambda () out))))

;; data ------------------------------------------------------------------------

;; Read data from the file at `path`, data is returned in the format:
;; ((RESULT INPUT*)*), that is, a list of clauses, where each clause is a list
;; containing a result followed by zero or more inputs.
(define (read-data path)
  (map (lambda (l) (map string->number (string-split l)))
       (with-input-from-file path read-lines)))

;; program ---------------------------------------------------------------------

(define-constant BF-CHARS "<>-+[],.")
(define-constant BF-CHAR-CNT 8)

(define-record-type program
  (%make-program genes score)
  program?
  (genes program-genes program-genes-set!)
  (score program-score program-score-set!))

(define (make-program genes #!optional score)
  (%make-program genes score))

(define-record-printer (program prog op)
  (fprintf op "#<program ~A ~S>" (program-score prog) (program-genes prog)))

(define (random-bf-char)
  (string-ref BF-CHARS (pseudo-random-integer BF-CHAR-CNT)))

(define (random-program gene-len)
  (define genes (make-string gene-len))
  (do ((i 0 (1+ i))) ((= i gene-len))
    (string-set! genes i (random-bf-char)))
  (make-program genes))

(define (score-program! prog data)
  (define genes (program-genes prog))
  (with-output-to-file TEMP-FILE-PATH (lambda () (print genes)))
  (define score 1) ; 1 is minimum score.
  (let loop ((ds data))
    (unless (null? ds)
      (let* ((d (car ds))
             (result (brainfuck-eval TEMP-FILE-PATH (cdr d) genes)))
        (when result
          (when (= result (car d))
            (set! score (1+ score)))
          (loop (cdr ds))))))
  (program-score-set! prog score))

;; Combine 2 parent programs to form 2 children. Each parent has a random point
;; chosen within it's bounds, and each half of the parent goes to one of the
;; children.
(define (combine-programs prog1 prog2)
  (define-values (g1 g2) (values (program-genes prog1) (program-genes prog2)))
  (define-values (p1 p2) (values (pseudo-random-integer (string-length g1))
                                 (pseudo-random-integer (string-length g2))))
  (list (make-program (string-append (substring g1 0 p1) (substring g2 p2)))
        (make-program (string-append (substring g2 0 p2) (substring g1 p1)))))

(define (mutate-program! prog)
  (define genes (program-genes prog))
  (define len (string-length genes))
  (define (add-cmd!)
    (program-genes-set! prog
                        (string-insert genes
                                       (pseudo-random-integer len)
                                       (random-bf-char))))
  (define (del-cmd!)
    (program-genes-set! prog (string-delete genes (pseudo-random-integer len))))
  (define (swap-cmd!) (string-swap! genes
                                    (pseudo-random-integer len)
                                    (pseudo-random-integer len)))
  (define (edit-cmd!) (string-set! genes
                                   (pseudo-random-integer len)
                                   (random-bf-char)))
  (case (pseudo-random-integer 4)
    ((0) (add-cmd!))
    ((1) (unless (zero? len) (del-cmd!)))
    ((2) (unless (< len 2) (swap-cmd!)))
    ((3) (unless (zero? len) (edit-cmd!)))))

;; population ------------------------------------------------------------------

(define (random-population pop-size max-gene-len)
  (vector-unfold (lambda (i) (random-program
                              (pseudo-random-integer (1+ max-gene-len))))
                 pop-size))

(define (score-population! pop data)
  (vector-for-each (lambda (i p) (score-program! p data)) pop))

(define (sort-population! pop)
  (sort! pop (lambda (p1 p2)
               (let ((s1 (program-score p1))
                     (s2 (program-score p2)))
                 (or (> s1 s2)
                     (and (= s1 s2)
                          (< (string-length (program-genes p1))
                             (string-length (program-genes p2)))))))))

;; Select `cnt` parents from `pop` using Roulette Wheel Selection.
(define (select-parents/roulette pop cnt)
  (define len (vector-length pop))
  (define total-score (vector-fold (lambda (i a p) (+ a (program-score p)))
                                   0.0 pop))
  ;; cumulative probability distribution
  (define cpd (vector-map (lambda (i prog) (/ (program-score prog) total-score))
                          pop))
  (let ((sum 0))
    (vector-for-each (lambda (i prob)
                       (set! sum (+ sum prob))
                       (vector-set! cpd i sum))
                     cpd))
  (define (select-parent)
    (define ran (pseudo-random-real))
    (let loop ((i 0))
      (cond ((= i len) (error 'select-parents/roulette "invalid CPD" cpd))
            ((<= ran (vector-ref cpd i)) i)
            (else (loop (1+ i))))))
  (define parents (make-vector PARENT-COUNT))
  (do ((i 0 (1+ i)))
      ((= i PARENT-COUNT))
    (vector-set! parents i
                 (vector-ref pop (select-parent))))
  parents)

;; main ------------------------------------------------------------------------

;; Note: At all times, population is kept sorted by score descending, equal
;; score programs are sorted by gene length ascending.

(define REPLACE-OFFSET (- POPULATION-SIZE CHILD-COUNT))

(define (next-generation! pop data)
  ;; Select PARENT-COUNT parents from population.
  (print "braingen: selecting parents...")
  (define parents (select-parents/roulette pop PARENT-COUNT))
  ;; Create children by combining each pair of 2 parents for 2 children.
  (print "braingen: creating offspring...")
  (define children (make-vector PARENT-COUNT))
  (do ((i 0 (2+ i)))
      ((>= i PARENT-COUNT))
    (let ((offspring (combine-programs (vector-ref parents i)
                                       (vector-ref parents (1+ i)))))
      (for-each mutate-program! offspring)
      (for-each (lambda (c) (score-program! c data)) offspring)
      (vector-set! children i (car offspring))
      (vector-set! children (1+ i) (cadr offspring))))
  ;; Replace least score programs in population with children.
  (print "braingen: replacing with children...")
  (do ((i 0 (1+ i)))
      ((= i CHILD-COUNT))
    (vector-set! pop (+ i REPLACE-OFFSET)
                 (vector-ref children i)))
  (print "braingen: sorting population...")
  (sort-population! pop))

(define (main data-path #!optional seed)
  (set-pseudo-random-seed! (if seed seed (random-bytes)))
  (define pop (random-population POPULATION-SIZE INIT-MAX-GENE-LEN))
  (define data (read-data data-path))
  (define target-score (length data))
  (print "braingen: scoring population...")
  (score-population! pop data)
  (print "braingen: sorting population...")
  (sort-population! pop)
  (let loop ((gen 0))
    (printf "[~A]: ~A\n" gen pop)
    (let ((prime (vector-ref pop 0)))
      (if (or (= gen MAX-GENERATIONS) (>= (program-score prime) target-score))
          (printf "BEST PROGRAM (SCORE: ~A, LENGTH: ~A): ~S\n"
                  (program-score prime)
                  (string-length (program-genes prime))
                  (program-genes prime))
          (begin (next-generation! pop data) (loop (1+ gen)))))))

(let ((args (command-line-arguments)))
  (if (= (length args) 1)
      (main (car args))
      (print "braingen: invalid arguments")))
