BINARIES := do_processing run_one_thread

CFLAGS := -std=c++11 -pthread -lrt -fPIC -march=native -mtune=native -mavx -mavx2 -Wall -Wno-format -I$(BOOST_INC)
LDFLAGS := -lboost_program_options -L$(BOOST_LIB)

CFLAGS_DEBUG := -g
CFLAGS_OPT := -O3 -ffast-math -g

ASM_FLAGS := -fverbose-asm -masm=intel -Wa,-adhln

LIBRARY := libfindprim.so

ifdef DEBUG
CFLAGS += $(CFLAGS_DEBUG)
else
CFLAGS += $(CFLAGS_OPT)
endif

all: $(BINARIES) $(LIBRARY) do_processing.S
.PHONY: all

%.o: %.c
	$(CC) $(INCLUDE) $(CFLAGS) -MMD -MP -o $@ -c $<

%.o: %.cxx
	$(CXX) $(INCLUDE) $(CFLAGS) -MMD -MP -o $@ -c $<

%.o: %.cpp
	$(CXX) $(INCLUDE) $(CFLAGS) -MMD -MP -o $@ -c $<

%.S: %.cxx
	$(CXX) $(INCLUDE) $(CFLAGS) $(ASM_FLAGS) -MMD -MP -o $@ -S $<

libfindprim.so: process_avx2.o process_naive.o design_fir.o
	$(CXX) -shared -o $@ $<

# All of the binaries have the same format, so use a "static pattern
# rule". Each binary "foo" depends on "foo.o" and we build it with the
# recipe given ($@ will be the name of the binary)
$(BINARIES) : %: %.o
	$(CXX) $(INCLUDE) $(CFLAGS) $(LDFLAGS) -o $@ $^

clean:
	rm -f $(BINARIES) $(LIBRARY) *.o *.d

.PHONY: clean

-include $(DEPS)
