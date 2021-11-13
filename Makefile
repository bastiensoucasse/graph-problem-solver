# Makefile

FILESPARS	= $(wildcard src/parser/src/*.c)
FILESSRC	= src/main/Graph.c src/main/Z3Tools.c
FILESBICON	= $(wildcard src/EdgeConProblem/*.c)
CC			= gcc
CFLAGS		= -g -Iinclude/main -Iinclude/EdgeConProblem -Isrc/parser/include -Isrc/parser -Isrc/EdgeConProblem
LDLIBS		= -lz3
OBJPARS		= $(FILESPARS:parser/src/%.c=build/%.o)
OBJSRC		= $(FILESSRC:src/main/%.c=build/%.o) $(FILESBICON:src/EdgeConProblem/%.c=build/%.o)
OBJNOTMAIN	= build/Parser.o build/Lexer.o $(OBJPARS) $(OBJSRC) 
OBJ			= $(OBJNOTMAIN) build/main.o

.PHONY: all
all: graphProblemSolver

graphProblemSolver: $(OBJ) 
		$(CC) $(CFLAGS) $(OBJ) $(LDLIBS) -o graphProblemSolver

graphProblemSolver-sol: $(OBJNOTMAIN)
		echo "#define SUBJECT" > toto
		cat toto src/main/main.c > src/main/new-main.c
		$(CC) $(CFLAGS) $(OBJNOTMAIN) src/main/new-main.c $(LDLIBS) -o graphProblemSolver-sol
		rm src/main/new-main.c
		rm toto

build/Lexer.o: src/parser/Lexer.c src/parser/Parser.c
		mkdir -p build
		$(CC) -c $(CFLAGS) $< -o $@

build/Parser.o: src/parser/Parser.c src/parser/Lexer.c
		mkdir -p build
		$(CC) -c $(CFLAGS) $< -o $@

src/parser/Lexer.c:	src/parser/Lexer.l 
		flex --header-file=src/parser/Lexer.h -o src/parser/Lexer.c src/parser/Lexer.l


src/parser/Parser.c:	src/parser/Parser.y src/parser/Lexer.c
		bison --defines=src/parser/Parser.h -o src/parser/Parser.c src/parser/Parser.y

build/%.o:	src/main/%.c 
		mkdir -p build
		$(CC) -c $(CFLAGS) $^ -o $@

build/%.o:	src/EdgeConProblem/%.c 
		mkdir -p build
		$(CC) -c $(CFLAGS) $^ -o $@

build/%.o:	src/parser/src/%.c
		mkdir -p build
		$(CC) -c $(CFLAGS) $^ -o $@

build/graphUsage.o: examples/graphUsage.c 
		mkdir -p build
		$(CC) -c $(CFLAGS) $^ -o $@

graphParser: build/Lexer.o build/Parser.o $(OBJPARS) build/Graph.o build/graphUsage.o
		$(CC) $(CFLAGS) $^ -o $@

build/Z3Example.o: examples/Z3Example.c 
		mkdir -p build
		$(CC) -c $(CFLAGS) $^ -o $@

Z3Example: build/Z3Example.o build/Z3Tools.o
		$(CC) $(CFLAGS) $^ $(LDLIBS) -o $@

.PHONY: doc
doc:
		doxygen doxygen.config
		rm -f doc.html
		ln -s doc/html/files.html doc.html

.PHONY: clean
clean:
		rm -f build/*.o *~ src/parser/Lexer.c src/parser/Lexer.h src/parser/Parser.c src/parser/Parser.h graphProblemSolver graphParser Z3Example doc.html
		rm -rf doc
