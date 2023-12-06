CC        = cc



#----------------------------------------------------------------------------
# Project specific path defintions.
#----------------------------------------------------------------------------
PROJDIR    = .

SRC        = $(PROJDIR)/src
INCL       = $(PROJDIR)/include
OBJ        = $(PROJDIR)/obj
BIN        = $(PROJDIR)/bin
TEST       = $(PROJDIR)/test

LDEFINES    = -DNO_RANGE_CHECK
LINCLUDES      =  -I$(INCL)

#----------------------------------------------------------------------------
# Target files
#----------------------------------------------------------------------------
RNAVIEW    = $(BIN)/rnaview

SRC_FILES = $(wildcard $(SRC)/*.c)

HEADER_FILES = $(wildcard (INCL)/*.h)

OBJ_FILES = $(patsubst $(SRC_FILES)/%.c, $(OBJ)/%.o,  $(SRC_FILES))

all: $(RNAVIEW)

CFLAGS  =  $(LINCLUDES) 

$(RNAVIEW) : $(OBJ_FILES) 
	$(CC) -g -Wall $(CFLAGS) -o $@ $(OBJ_FILES) $(LDFLAGS) -lm $(MALLOCLIB)

$(OBJ)/%.o : $(SRC)/%.c 
	mkdir -o $(@D)
	$(CC) -g -Wall $(CFLAGS) -c $^ -o $@

clean:
	@rm -f $(OBJ)/*.o
	@rm -f $(ALLTARGETS)

export:
	mkdir -p $(EXPORT_DIR)
	@cd $(EXPORT_DIR); mkdir -p $(INCL)
	$(EXPORT) $(EXPORT_LIST) $(HEADER_FILES) $(EXPORT_DIR)/$(INCL)
	@cd $(EXPORT_DIR); mkdir -p $(SRC)
	$(EXPORT) $(EXPORT_LIST) $(SRC_FILES) $(TARGETSRC) $(EXPORT_DIR)/$(SRC)
	@cd $(EXPORT_DIR); mkdir -p $(BIN)
	@cd $(EXPORT_DIR); mkdir -p $(LIB)
	@cd $(EXPORT_DIR); mkdir -p $(OBJ)
	@cd $(EXPORT_DIR); mkdir -p $(TEST)
	@cp $(TESTFILES) $(EXPORT_DIR)/$(TEST)
	@cp Makefile $(EXPORT_DIR)




