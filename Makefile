.PHONY: all
all:
	csc -o braingen -O5 -d0 braingen.scm

.PHONY: debug
debug:
	csc -o braingen -d3 braingen.scm
