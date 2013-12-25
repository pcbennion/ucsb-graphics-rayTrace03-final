# calls:
CC         = g++
#CC         = /usr/lib64/ccache/x86_64-redhat-linux-g++
CFLAGS     = -c -Wall -O3 -I./

#linux
#LDFLAGS    = -lfreeglut -lopengl32
LDFLAGS	   = -lGL -lGLU -lglut
#osx
#LDFLAGS    = -framework OpenGL -framework GLUT
EXECUTABLE = skeletonViewer

SOURCES    = skeletonViewer.cpp skeleton.cpp motion.cpp displaySkeleton.cpp MATRIX4.cpp Shapes.cpp
OBJECTS    = $(SOURCES:.cpp=.o)

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f *.o
