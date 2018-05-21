OPT = -O2 -march=native
CFLAGS = -c -fpic -Wall -std=gnu99 $(OPT)
LDFLAGS = -shared

SOURCES = libhpf.c
TARGET = ../libhpf.so
HEADERS = $(SOURCES:.c=.h)
OBJECTS = $(SOURCES:.c=.o)

.PHONY : all clean
all: $(TARGET)

clean:
	rm -f $(OBJECTS) $(TARGET)

$(TARGET): $(OBJECTS)
	$(CC) $(LDFLAGS) -o $(TARGET) $(OBJECTS)

%.o: %.c
	$(CC) $(CFLAGS) $< -o $@
