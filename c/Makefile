CC = cc
OPT = -O2 -march=native
CFLAGS = -Wall -std=gnu99 $(OPT)

targets = hpf

.PHONY : default
default : all

.PHONY : all
all : clean $(targets)

hpf : hpf.c ../core/hpfcore.c
	$(CC) -o $@ $^

%.o : %.c
	$(CC) -c $(CFLAGS) $<

.PHONY : clean
clean:
	rm -f $(targets)
