# As long as you have wx-config on your system-wide directory 
# (through something like a make install of wxwidgets)
# this should work fine!

CXX = clang++
STD = -std=c++17
PROJECT = MDsimAPP
OBJECTS = MDsim-functions.o MDsim-app.o
WXFLAGS = `wx-config --cxxflags`
WXLIBS = `wx-config --libs`

$(PROJECT): $(OBJECTS)
	$(CXX) $(STD) $(OBJECTS) -o $(PROJECT) $(WXLIBS)

MDsim-app.o:
	$(CXX) $(STD) MDsim-app.cpp -c $(WXFLAGS)

MDsim-functions.o:
	$(CXX) $(STD) MDsim-functions.cpp -c

clean:
	-rm -f *.o

