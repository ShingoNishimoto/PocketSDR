USE_FFTW ?= 0

all:
	make -C pocket_scan USE_FFTW=$(USE_FFTW)
	make -C pocket_conf USE_FFTW=$(USE_FFTW)
	make -C pocket_dump USE_FFTW=$(USE_FFTW)
	make -C pocket_acq USE_FFTW=$(USE_FFTW)
	make -C pocket_trk USE_FFTW=$(USE_FFTW)
	make -C pocket_snap USE_FFTW=$(USE_FFTW)
	make -C pocket_calib USE_FFTW=$(USE_FFTW)
	make -C convbin USE_FFTW=$(USE_FFTW)
	make -C str2str USE_FFTW=$(USE_FFTW)
clean:
	make -C pocket_scan clean
	make -C pocket_conf clean
	make -C pocket_dump clean
	make -C pocket_acq clean
	make -C pocket_trk clean
	make -C pocket_snap clean
	make -C pocket_calib clean
	make -C convbin clean
	make -C str2str clean
install:
	make -C pocket_scan install USE_FFTW=$(USE_FFTW)
	make -C pocket_conf install USE_FFTW=$(USE_FFTW)
	make -C pocket_dump install USE_FFTW=$(USE_FFTW)
	make -C pocket_acq install USE_FFTW=$(USE_FFTW)
	make -C pocket_trk install USE_FFTW=$(USE_FFTW)
	make -C pocket_snap install USE_FFTW=$(USE_FFTW)
	make -C pocket_calib install USE_FFTW=$(USE_FFTW)
	make -C convbin install USE_FFTW=$(USE_FFTW)
	make -C str2str install USE_FFTW=$(USE_FFTW)
