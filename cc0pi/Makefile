all: univmake cc0pi_analyzer 

XGBOOST=/exp/uboone/app/users/mastbaum/xgboost/build

univmake: univmake.C 
	$(CXX) $(shell root-config --cflags --libs) -O3 -o $@ $^

cc0pi_analyzer: cc0pi_analyzer.C
	$(CXX) $(shell root-config --cflags --libs) -I${XGBOOST}/include -I.. -L${XGBOOST}/lib64 -lxgboost -g -o $@ $^

.PHONY: clean

clean:
	$(RM) univmake cc0pi_analyzer
