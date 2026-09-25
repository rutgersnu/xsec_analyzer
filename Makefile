CXX = g++
CXXFLAGS = -Wall -pedantic -Wno-unused-variable -Isrc
ROOTFLAGS = `root-config --cflags --libs`
LIB_DIR = ./lib

PLOTUTILS_SRCS = src/plotutils/PlotUtils.cpp src/plotutils/HistUtils.cpp src/plotutils/GridCanvas.cpp src/plotutils/MnvColors.cpp src/plotutils/HistFolio_slim.cpp src/plotutils/UBTH2Poly.cpp

all: dirs bin/chi_square_cc0pi_christian bin/univmake

dirs:
	@mkdir -p bin lib

src/%: %.cpp
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o  $@ $<
	
%.o : %.cpp
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $*.o  -c $*.cpp 
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $* $*.o

stv_root_dict.o:
	$(RM) stv_root_dict*.*
	rootcling -f stv_root_dict.cc -c LinkDef.h
	$(CXX) $(shell root-config --cflags --libs) -O3 \
	-fPIC -o stv_root_dict.o -c stv_root_dict.cc
	$(RM) stv_root_dict.cc
	
$(LIB_DIR)/libPlotUtils.so: 
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -shared -fPIC -O3 -o $@ $^ $(PLOTUTILS_SRCS)

bin/chi_square_cc0pi_christian: src/chi_square_cc0pi_christian.cpp $(LIB_DIR)/libPlotUtils.so
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -O3 -L$(LIB_DIR) -o $@ $^ -lPlotUtils

bin/annie_stv_prep: src/annie_stv_prep.cpp
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -O3 -o $@ $^

bin/slice_plots_ccinc: ccinc/slice_plots_ccinc.cpp
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -O3 -o $@ $^

bin/univmake: src/univmake.C
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -O3 -o $@ $^
	
.PHONY: clean

.INTERMEDIATE: stv_root_dict.o

clean:
	rm -f $(wildcard *.o) $(patsubst %.cpp, %, $(wildcard *.cpp)) chi_square_cc0pi_christian univmake $(LIB_DIR)/*.so bin/*

