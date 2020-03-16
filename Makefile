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

#EXEC2	     = post_event
#COMPONENTS2  = main_post_event post_event
#BIN_EXEC2    = $(addprefix $(BIN)/,$(EXEC2) )

#EXEC3        = daq
#COMPONENTS3  = main_daq daq
#BIN_EXEC3     = $(addprefix $(BIN)/,$(EXEC3) )

EXEC4	     = ext
COMPONENTS4  = main_ext ext
BIN_EXEC4    = $(addprefix $(BIN)/,$(EXEC4) )

EXEC5        = clean
COMPONENTS5  = main_clean clean
BIN_EXEC5    =$(addprefix $(BIN)/,$(EXEC5) )

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
# 2
#$(BIN_EXEC2):
#	$(LD) $(LDFLAGS) $^ $(ROOTLIBS) $(addprefix -L, $(LIBS)) -o $@
#$(BIN)/$(EXEC2): $(addprefix $(LIB)/, $(addsuffix .o, $(COMPONENTS2) ) )
# 3
#$(BIN_EXEC3):
#	$(LD) $(LDFLAGS) $^ $(ROOTLIBS) $(addprefix -L, $(LIBS)) -o $@
#$(BIN)/$(EXEC3): $(addprefix $(LIB)/, $(addsuffix .o, $(COMPONENTS3) ) )
#4
$(BIN_EXEC4):
	$(LD) $(LDFLAGS) $^ $(ROOTLIBS) $(addprefix -L, $(LIBS)) -o $@
$(BIN)/$(EXEC4): $(addprefix $(LIB)/, $(addsuffix .o, $(COMPONENTS4) ) )
#5
$(BIN_EXEC5):
	$(LD) $(LDFLAGS) $^ $(ROOTLIBS) $(addprefix -L, $(LIBS)) -o $@
$(BIN)/$(EXEC5): $(addprefix $(LIB)/, $(addsuffix .o, $(COMPONENTS5) ) )

.PHONY : clean
clean:
	@echo ... cleaning
	rm -f $(LIB)/*.o
	rm -f $(LIB)/*.so
	rm -f $(BIN_EXEC0)
	rm -f $(BIN_EXEC1)
#rm -f $(BIN_EXEC2)
#rm -f $(BIN_EXEC3)
	rm -f $(BIN_EXEC4)
	rm -f $(BIN_EXEC5)
cleanext:
	@echo ... cleaning ext
	rm -f $(BIN_EXEC4)
$(LIB):
	mkdir -p $(LIB)
$(BIN):
	mkdir -p $(BIN)
installdirs: $(LIB) $(BIN)
rec:  installdirs $(BIN_EXEC0) installdirs $(BIN_EXEC1) installdirs $(BIN_EXEC4) installdirs $(BIN_EXEC5)
all: rec

