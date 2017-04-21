OPT = -O2 -march=native
CFLAGS = -c -fpic -Wall -std=gnu99 $(OPT)
LDFLAGS =

SOURCES = hpf.c ../core/libhpf.c
TARGET = hpf
HEADERS = $(SOURCES:.c=.h)
OBJECTS = $(SOURCES:.c=.o)

.PHONY : all clean
all: $(TARGET)

clean:
	rm -f $(OBJECTS) $(TARGET)

$(TARGET): $(OBJECTS)
	$(CC) $(LDFLAGS) -o $(TARGET) $(OBJECTS)

all: $(TARGET)
	chmod +x $(TARGET)

%.o: %.c
	$(CC) $(CFLAGS) $< -o $@
