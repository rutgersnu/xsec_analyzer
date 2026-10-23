CXX = g++
CXXFLAGS = -Wall -pedantic -Wno-unused-variable -Isrc
ROOTFLAGS = `root-config --cflags --libs`
LIB_DIR = ./lib

PLOTUTILS_SRCS = src/plotutils/PlotUtils.cpp src/plotutils/HistUtils.cpp src/plotutils/GridCanvas.cpp src/plotutils/MnvColors.cpp src/plotutils/HistFolio_slim.cpp src/plotutils/UBTH2Poly.cpp

all: dirs bin/chi_square_cc0pi_christian bin/univmake bin/annie_stv_prep bin/slice_plots_ccinc

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

# Libraries after the sources: linkers using --as-needed (Ubuntu default) drop libraries listed first
bin/univmake: src/univmake.C src/UniverseMaker.hh src/EventCategory.hh src/SystematicsCalculator.hh src/MCC9SystematicsCalculator.hh src/FilePropertiesManager.hh
	$(CXX) $(CXXFLAGS) -O3 -o $@ $< $(ROOTFLAGS)

# nuNCQE single-bin analysis (see ncqe/README_NCQE.md). -ffp-contract=off keeps the selection
# arithmetic identical to the Python reference (no fused multiply-adds).
NCQE_FLAGS = -O2 -ffp-contract=off

ncqe: dirs bin/univmake bin/annie_ncqe_prep bin/ncqe_xsec

bin/annie_ncqe_prep: src/annie_ncqe_prep.cpp src/NCQESelection.hh src/EventCategory.hh
	$(CXX) $(CXXFLAGS) $(NCQE_FLAGS) -o $@ $< $(ROOTFLAGS)

bin/ncqe_xsec: src/ncqe_xsec.cpp src/NCQESystematicsCalculator.hh src/SystematicsCalculator.hh src/FilePropertiesManager.hh src/UniverseMaker.hh src/EventCategory.hh
	$(CXX) $(CXXFLAGS) $(NCQE_FLAGS) -o $@ $< $(ROOTFLAGS)

.PHONY: clean

.INTERMEDIATE: stv_root_dict.o

clean:
	rm -f $(wildcard *.o) $(patsubst %.cpp, %, $(wildcard *.cpp)) chi_square_cc0pi_christian univmake $(LIB_DIR)/*.so bin/*

