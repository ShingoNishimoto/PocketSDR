# When invoked under sudo, HOME becomes /root. Use the original user's home instead.
ifdef SUDO_USER
  _REAL_HOME := $(shell getent passwd $(SUDO_USER) | cut -d: -f6)
else
  _REAL_HOME := $(HOME)
endif
BINDIR   ?= /usr/local/bin
BINS     = $(wildcard bin/pocket_* bin/fftw_wisdom)
LIBDIR   = lib/build
USE_FFTW ?= 1

all:
	$(MAKE) -C $(LIBDIR) USE_FFTW=$(USE_FFTW)
	$(MAKE) -C $(LIBDIR) install USE_FFTW=$(USE_FFTW)
	$(MAKE) -C app USE_FFTW=$(USE_FFTW)

clean:
	$(MAKE) -C $(LIBDIR) clean
	$(MAKE) -C app clean

install:
	mkdir -p $(BINDIR)
	$(MAKE) -C app install BIN=$(abspath $(BINDIR)) USE_FFTW=$(USE_FFTW)

uninstall:
	cd $(BINDIR) && rm -f $(notdir $(BINS))
