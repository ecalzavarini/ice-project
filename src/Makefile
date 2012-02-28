CC=gcc 
CFLAGS= -Wall  
LFLAGS= -lm 
SOURCES= parameters.c fluid.c initial_conditions.c temperature.c output.c lb.c salt.c melting.c
OBJECTS=$(SOURCES:.c=.o)
HEADERS= common_object.h common_main.h typedef.h 
EXECUTABLE= ice2d

%.o: %.c $(HEADERS) 
	$(CC) -c $(CFLAGS) $< -o $@ 

$(EXECUTABLE):  main.c $(OBJECTS) $(HEADERS)
	$(CC) $(CFLAGS) -o $(EXECUTABLE) main.c $(OBJECTS) $(LFLAGS)

clean:
	rm -rf ./*.o ./*~ $(EXECUTABLE) ./*.dat ./RUN

tar:
	tar cvzf iceproject_backup.tgz ./*.c ./*.h ./*.sh Makefile 