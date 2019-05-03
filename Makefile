SRC  = src
LIB  = lib
BIN  = bin
INC  = $(addprefix -I,$(SRC))

vpath %.C $(SRC)

CC      = g++
CCFLAGS = -c  -Wno-deprecated -g  #-std=c++11 -Wall 

LD      = g++
#LDFLAGS = -lGui -lTreePlayer -Wall -std=c++11 -fopenmp
LDFLAGS = -lGui -lTreePlayer -Wall#-w # -Wall  -w toglie i warning
# YOU SHOULD ACTIVE -Wall  every time !!!!

ROOTINCS = $(shell root-config --cflags) # will be used as it is
ROOTLIBS = $(shell root-config --glibs) -lMinuit  # will be used as it is
LIBS     = $(LIB)                        # will get an -L for each entry


EXEC0        = ana
COMPONENTS0  = main_ana ana
BIN_EXEC0    = $(addprefix $(BIN)/,$(EXEC0) )

EXEC1        = event
COMPONENTS1  = main_event event
BIN_EXEC1    = $(addprefix $(BIN)/,$(EXEC1) )

default: all

# compile sources
$(LIB)/%.o: %.C
	@echo .
	@echo ... compiling source: $< to $@
	$(CC) $(CCFLAGS) $< $(ROOTINCS) $(INC) -o $@

# 0
$(BIN_EXEC0):
	$(LD) $(LDFLAGS) $^ $(ROOTLIBS) $(addprefix -L, $(LIBS)) -o $@
$(BIN)/$(EXEC0): $(addprefix $(LIB)/, $(addsuffix .o, $(COMPONENTS0) ) )
# 1  
$(BIN_EXEC1):
	$(LD) $(LDFLAGS) $^ $(ROOTLIBS) $(addprefix -L, $(LIBS)) -o $@
$(BIN)/$(EXEC1): $(addprefix $(LIB)/, $(addsuffix .o, $(COMPONENTS1) ) )

.PHONY : clean
clean:
	@echo ... cleaning
	rm -f $(LIB)/*.o
	rm -f $(LIB)/*.so
	rm -f $(BIN_EXEC0)
	rm -f $(BIN_EXEC1)
$(LIB):
	mkdir -p $(LIB)
$(BIN):
	mkdir -p $(BIN)
installdirs: $(LIB) $(BIN)
rec:  installdirs $(BIN_EXEC0) installdirs $(BIN_EXEC1) 
all: rec

