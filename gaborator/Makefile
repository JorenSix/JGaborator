
all:
	# Performance can be significantly improved by using a faster FFT library. 
	# On macOS, you can use the FFT from Apple's vDSP library by defining GABORATOR_USE_VDSP and linking with the Accelerate framework
	# On Linux and NetBSD, you can use the PFFFT (Pretty Fast FFT) library which is what is done below: 
	cc -c -O3 -ffast-math -fPIC pffft/pffft.c -o pffft/pffft.o
	cc -c -O3 -ffast-math -fPIC -DFFTPACK_DOUBLE_PRECISION pffft/fftpack.c -o pffft/fftpack.o
	c++ -std=c++11 -I"gaborator-1.2" -I"pffft" -I"$(JAVA_HOME)/include" -I"$(JAVA_HOME)/include/linux" -fPIC -shared -O3 -ffast-math -DGABORATOR_USE_PFFFT  -o ../build/precompiled/libjgaborator.so jgaborator.cc pffft/pffft.o pffft/fftpack.o	

mac:
	c++ -std=c++11 -I"gaborator-1.2"  -I"$(JAVA_HOME)/include" -I"$(JAVA_HOME)/include/darwin" -O3 -ffast-math -DGABORATOR_USE_VDSP  -o ../build/precompiled/libjgaborator.dylib jgaborator.cc  -framework Accelerate

clean: 
	rm *.so *.class
