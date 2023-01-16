
JNI_INCLUDES=-I"$(JAVA_HOME)/include" -I"$(JAVA_HOME)/include/darwin" -I"$(JAVA_HOME)/include/linux" -I"$(JAVA_HOME)/include/win32" -I"$(JAVA_HOME)/include/win64"

all:
	echo $(JAVA_HOME)
	ruby precompile_for_targets.rb

plain:
	# Performance can be significantly improved by using a faster FFT library. 
	# On macOS, you can use the FFT from Apple's vDSP library by defining GABORATOR_USE_VDSP and linking with the Accelerate framework
	# On Linux and NetBSD, you can use the PFFFT (Pretty Fast FFT) library which is what is done below: 
	cc -c -O3 -ffast-math -fPIC pffft/pffft.c -o pffft/pffft.o
	cc -c -O3 -ffast-math -fPIC -DFFTPACK_DOUBLE_PRECISION pffft/fftpack.c -o pffft/fftpack.o
	c++ -std=c++11 -I"gaborator-1.7" -I"pffft" $(JNI_INCLUDES) -fPIC -shared -O3 -ffast-math -DGABORATOR_USE_PFFFT  -o ../src/main/resources/jni/libjgaborator.so jgaborator.cc pffft/pffft.o pffft/fftpack.o

mac_aarch64:
	c++ -std=c++11 -I"gaborator-1.7"  $(JNI_INCLUDES) -O3 -ffast-math -DGABORATOR_USE_VDSP  -o ../src/main/resources/jni/macos_aarch64_libjgaborator.dylib jgaborator.cc  -framework Accelerate

mac_dummy:
	c++ -std=c++11 -I"gaborator-1.7"  $(JNI_INCLUDES) -O3 -ffast-math -DGABORATOR_USE_VDSP  -o .../src/main/resources/jni/libjgaborator.dylib jgaborator_dummy.cc  -framework Accelerate

clean:
	rm *.so *.class pfft/*o

zig_windows_example:
	echo $(JNI_INCLUDES)
	zig cc -target x86_64-windows-gnu -c -O3 -ffast-math -fPIC pffft/pffft.c -o pffft/pffft.o
	zig cc -target x86_64-windows-gnu -c -O3 -ffast-math -fPIC -DFFTPACK_DOUBLE_PRECISION pffft/fftpack.c -o pffft/fftpack.o
	zig c++ -target x86_64-windows-gnu -I"pffft" -I"gaborator-1.7"  $(JNI_INCLUDES) -O3 -ffast-math -DGABORATOR_USE_PFFFT  -o libjgaborator.dll jgaborator.cc pffft/pffft.o pffft/fftpack.o

zig:
	zig cc  -target $(ZIG_TARGET_PLATFORM) -c -O3 -ffast-math -fPIC pffft/pffft.c -o pffft/pffft.o
	zig cc  -target $(ZIG_TARGET_PLATFORM) -c -O3 -ffast-math -fPIC -DFFTPACK_DOUBLE_PRECISION pffft/fftpack.c -o pffft/fftpack.o
	zig c++ -target $(ZIG_TARGET_PLATFORM) -I"pffft" -I"gaborator-1.7"  $(JNI_INCLUDES)  -shared -O3 -ffast-math -DGABORATOR_USE_PFFFT  -o $(ZIG_OUTPUT_FILE) jgaborator.cc pffft/pffft.o pffft/fftpack.o



