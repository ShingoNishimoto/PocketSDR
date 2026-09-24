#
#  makefile of LIBFEC shared library (libfec.so)
#
#! You need to install LIBFEC source tree as follows.
#!
#! $ git clone https://github.com/quiet/libfec libfec

#! specify directory of LIBFEC source tree
SRC = ../libfec

ifeq ($(OS),Windows_NT)
    INSTALL = ../win32
    EXT = so
else ifeq ($(shell uname -sm),Darwin arm64)
    INSTALL = ../macos
    EXT = dylib
else
    INSTALL = ../linux
    EXT = so
endif

ifeq ($(shell uname -m),aarch64)
    CONF_OPT = --build=arm
endif

TARGET = libfec.$(EXT) libfec.a

all :
	DIR=`pwd`; \
	mkdir -p $(SRC)/build; \
	cd $(SRC)/build; \
	../configure $(CONF_OPT); \
	sed 's/-lc//' < makefile > makefile.p; \
	mv makefile.p makefile; \
	make; \
	cd $$DIR; \
    cp $(SRC)/build/libfec.a $(SRC)/build/libfec.$(EXT) .

clean:
	rm -rf $(SRC)/build
	rm -f $(TARGET)

install:
	cp $(TARGET) $(INSTALL)
