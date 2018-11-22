CXX ?= g++
PKG_CXXFLAGS = -O3 -march=native -std=c++11 -pthread -fPIC -fopenmp
PKG_CXXFLAGS += -I"./include" -I"./igraph-R"
LDFLAGS += -lpthread
#CXXFLAGS += -I"./include" -I"./igraph-R"

LIB = ../libleidenalg.a
lib: $(LIB)

$(LIB): CPMVertexPartition.o LinearResolutionParameterVertexPartition.o RBConfigurationVertexPartition.o SurpriseVertexPartition.o ModularityVertexPartition.o RBERVertexPartition.o GraphHelper.o MutableVertexPartition.o ResolutionParameterVertexPartition.o Optimiser.o SignificanceVertexPartition.o
	ar rvs $@ $?

#libn2.so: base.o hnsw.o hnsw_node.o distance.o heuristic.o mmap.o
#	$(CXX) $(CXXFLAGS) -shared -o $@ $(LDFLAGS) $?

clean:
	@-rm -f *.o $(LIB) 


