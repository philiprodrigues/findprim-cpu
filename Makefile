BINARIES := do_processing do_processing_store do_processing_transpose do_processing_32tap do_processing_coherent_pedsub do_processing_coherent_pedsub_store

CFLAGS := -pthread -lrt -fPIC -mavx -mavx2 -Wall -Wno-format -I$(BOOST_INC)
LDFLAGS := -lboost_program_options -L$(BOOST_LIB)

CFLAGS_DEBUG := -g
CFLAGS_OPT := -O3 -ffast-math

ifdef DEBUG
CFLAGS += $(CFLAGS_DEBUG)
else
CFLAGS += $(CFLAGS_OPT)
endif

all: $(BINARIES)
.PHONY: all

%.o: %.c
	$(CC) $(INCLUDE) $(CFLAGS) -MMD -MP -o $@ -c $<

%.o: %.cxx
	$(CXX) $(INCLUDE) $(CFLAGS) -MMD -MP -o $@ -c $<

do_processing_store.o: do_processing.cxx
	$(CXX) $(INCLUDE) $(CFLAGS) -DSTORE_INTERMEDIATE -MMD -MP -o $@ -c $<

do_processing_transpose.o: do_processing.cxx
	$(CXX) $(INCLUDE) $(CFLAGS) -DTRANSPOSE_MEMORY -MMD -MP -o $@ -c $<

do_processing_32tap.o: do_processing.cxx
	$(CXX) $(INCLUDE) $(CFLAGS) -DNTAPS=32 -MMD -MP -o $@ -c $<

do_processing_coherent_pedsub.o: do_processing.cxx
	$(CXX) $(INCLUDE) $(CFLAGS) -DCOHERENT_PEDSUB -MMD -MP -o $@ -c $<

do_processing_coherent_pedsub_store.o: do_processing.cxx
	$(CXX) $(INCLUDE) $(CFLAGS) -DCOHERENT_PEDSUB -DSTORE_INTERMEDIATE -MMD -MP -o $@ -c $<

# All of the binaries have the same format, so use a "static pattern
# rule". Each binary "foo" depends on "foo.o" and we build it with the
# recipe given ($@ will be the name of the binary)
$(BINARIES) : %: %.o
	$(CXX) $(INCLUDE) $(CFLAGS) $(LDFLAGS) -o $@ $^

clean:
	rm -f $(BINARIES) *.o *.d

.PHONY: clean

-include $(DEPS)