# Makefile

BINDIR := bin
LIBDIR := lib

CC := g++ -std=c++11 -pedantic -O2 -g -pipe -Wall

MP3PLAYER = mpg123

# src/ (declaracoes de funcoes, de classes + codigo)
# main/ (programas principais)
# bin/ (temporarios, .o, .exe)
# lib/ (bibliotecas) biblioteca

# making library
# - static: .a
# - shareable: .so

VPATH = main:src

ROOTLIB := $(shell root-config --libs)
ROOTINC := $(shell root-config --incdir)

SRC := $(wildcard src/*.cpp) $(wildcard src/*.C)
OBJ := $(patsubst %.cpp, $(BINDIR)/%.o, $(notdir $(SRC))) $(patsubst %.C, $(BINDIR)/%.o, $(notdir $(SRC)))
PYOBJ := $(wildcard graphics/*.png) $(wildcard graphics/*.mp4) 
INC := $(wildcard src/*.h)

lib: $(LIBDIR)/libFT.a

$(LIBDIR)/libFT.a: $(OBJ)
	@echo Creating FireTracer lib
	ar ruv $@ $^
	ranlib $@
	@echo ----------libFT done-------------

%.exe: $(BINDIR)/%.o $(LIBDIR)/libFT.a 
	@echo ------------------------"\nLinking Thou Command:\n\t$<\n"------------------------
	$(CC) -I src $< -o $(BINDIR)/$@ -L lib -l FT $(ROOTLIB)

$(BINDIR)/%.o: %.cpp | $(INC)
	@echo ------------------------"\nCompiling Thou Command:\n\t$<\n"------------------------
	$(CC) -I src -I $(ROOTINC) -c $< -o $@

$(BINDIR)/%.o: %.cc | $(INC)
	@echo ------------------------"\nCompiling Thou Command:\n\t$<\n"------------------------
	$(CC) -I src -I $(ROOTINC) -c $< -o $@

$(BINDIR)/%.o: %.C | $(INC)
	@echo ------------------------"\nCompiling Thou Command:\n\t$<\n"------------------------
	$(CC) -I src -I $(ROOTINC) -c $< -o $@


#################################################
### if you want to do things the funny way :D ###
#################################################

ussr:
	$(MP3PLAYER) -q -n 700 soundtrack/hymn.mp3

rui:
	$(MP3PLAYER) -q soundtrack/rui.mp3

lib_ussr:
	gnome-terminal -- $(MP3PLAYER) -q -n 1790 soundtrack/hymn.mp3
	make lib

exe:
	@read -p "Enter .exe name:" module; \
	executavel=./bin/$$module.exe; \
	make rui; \
	$$executavel

py:
	@echo Preparing some visualization
	python3 python/graphicmaker.py
	python3 python/gifmaker.py

################## clean #########################

tilde := $(wildcard */*~) $(wildcard *~)
exe := $(wildcard $(BINDIR)/*.exe) $(wildcard *.exe)
obj := $(wildcard */*.o) $(wildcard *.o) $(wildcard */*.so) $(wildcard */*.pcm) $(wildcard */*.d)

clean:
	@echo Cleaning
	rm -f $(exe) $(obj) $(tilde) $(generated)
	rm -f $(LIBDIR)/libFT.a
	./Cleaner

cclean:
	@echo Major Cleaning
	rm -f $(exe) $(obj) $(tilde) $(generated)
	rm -f $(LIBDIR)/libFT.a
	rm -f $(PYOBJ)
	./Cleaner

cleanpy:
	@echo Cleaning pictures
	rm -f $(PYOBJ)
	./Cleaner