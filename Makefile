
# ------------------------------------------
# change project path here !!!
# ------------------------------------------

#PATHPROJ = ./projects/test_nd
PATHPROJ = ./tests/eigenvectors_kempe

# ------------------------------------------


PATHCORE = ./core

CXX=g++ -g
RM=rm -f
CPPFLAGS=$(shell root-config --cflags)
LDFLAGS=$(shell root-config --ldflags)
LDLIBS=$(shell root-config --libs)



include $(PATHPROJ)/Makefile



OBJS=$(subst .cc,.o,$(SRCS))
LOBJS=$(patsubst %, build/%, $(notdir $(OBJS)))

MOBJS=$(subst .cc,.o,$(MAINS))
LMOBJS=$(patsubst %, build/%, $(notdir $(MOBJS)))

all: clean sim cleanobjs

run: runsim 

sim: $(OBJS) $(MOBJS)
	$(CXX) $(LDFLAGS) -o sim $(LOBJS) ./build/main.o

%.o: %.cc
	$(CXX) -std=c++11 -c $< -o build/$(notdir $@)
	

clean:
	$(RM) build/*.o
	$(RM) sim

cleanobjs:
	$(RM) build/*.o

runsim:
	./sim

doc:
	doxygen doxygen_config

