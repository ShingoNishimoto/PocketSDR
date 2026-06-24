BINDIR   = /usr/local/bin
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
	$(MAKE) -C app install
	install -m 755 $(BINS) $(BINDIR)

uninstall:
	cd $(BINDIR) && rm -f $(notdir $(BINS))
