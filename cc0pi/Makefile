all: cc0pi_analyzer

XGBOOST=/exp/uboone/app/users/mastbaum/xgboost/build

cc0pi_analyzer: analyzer.C
	$(CXX) $(shell root-config --cflags --libs) -I${XGBOOST}/include -I.. -L${XGBOOST}/lib64 -lxgboost -g -o $@ $^ ../stv_root_dict.o

.PHONY: clean

clean:
	$(RM) cc0pi_analyzer

