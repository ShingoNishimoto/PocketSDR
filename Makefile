# When invoked under sudo, HOME becomes /root. Use the original user's home instead.
ifdef SUDO_USER
  _REAL_HOME := $(shell getent passwd $(SUDO_USER) | cut -d: -f6)
else
  _REAL_HOME := $(HOME)
endif
BINDIR   ?= $(_REAL_HOME)/bin
BINS     = $(wildcard bin/pocket_* bin/fftw_wisdom)
LIBDIR   = lib/build

all:
	$(MAKE) -C $(LIBDIR)
	$(MAKE) -C $(LIBDIR) install
	$(MAKE) -C app

clean:
	$(MAKE) -C $(LIBDIR) clean
	$(MAKE) -C app clean

install:
	mkdir -p $(BINDIR)
	$(MAKE) -C app install BIN=$(abspath $(BINDIR))
ifdef SUDO_USER
	chown -R $(SUDO_USER):$(SUDO_USER) $(abspath $(BINDIR))/pocket_* \
	    $(abspath $(BINDIR))/fftw_wisdom $(abspath $(BINDIR))/convbin 2>/dev/null || true
endif

uninstall:
	cd $(BINDIR) && rm -f $(notdir $(BINS))
