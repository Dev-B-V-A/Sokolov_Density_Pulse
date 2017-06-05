CWFLAGS  = -mfpmath=sse -W -Wall -Wextra -Wunused -Wcast-align -Werror -pedantic -pedantic-errors \
-fstack-protector-all -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute \
-Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long -Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual
CPROFILE = 
CDEBUG   = $(CPROFILE)
CC       = gcc -c
LD       = gcc
CFLAGS   = $(DEBUG)
LIBS     = -L./lib -llaspack -lxc -lm
LDFLAGS  = $(DEBUG)
INC	 = -I./include

OBJS = function.o gas_params.o norma.o report.o init.o fill_in.o solver.o scheme.o main.o

all : scheme.out

scheme.out : $(OBJS)
	$(LD) $(LDFLAGS) $(INC) $^ $(LIBS) -o $@ -O3 -lpthread --fast-math
function.o : function.cpp function.h
	$(CC) $(CFLAGS) $< -o $@ -O3 --fast-math
gas_params.o : gas_params.cpp gas_params.h
	$(CC) $(CFLAGS) $< -o $@ -O3 --fast-math
norma.o : norma.cpp norma.h
	$(CC) $(CFLAGS) $< -o $@ -O3 --fast-math
report.o : report.cpp report.h
	$(CC) $(CFLAGS) $< -o $@ -O3 --fast-math
init.o: init.cpp init.h
	$(CC) $(CFLAGS) $< -o $@ -O3 --fast-math
fill_in.o: fill_in.cpp fill_in.h
	$(CC) $(CFLAGS) $(INC) $< -o $@ -O3 --fast-math
solver.o : solver.cpp solver.h
	$(CC) $(CFLAGS) $(INC) $< -o $@ -O3 --fast-math
scheme.o : scheme.cpp scheme.h
	$(CC) $(CFLAGS) $(INC) $< -o $@ -O3 --fast-math
main.o : main.cpp
	$(CC) $(CFLAGS) $(INC) $< -o $@ -O3 --fast-math
clean: 
	rm -f $(OBJS) a.out
clean_rand:
	rm -rf rand.txt
	
