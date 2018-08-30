//
// Fast Fourier transform using the Apple vDSP framework
//
// Copyright (C) 2013-2018 Andreas Gustafsson.  This file is part of
// the Gaborator library source distribution.  See the file LICENSE at
// the top level of the distribution for license information.
//

#ifndef _GABORATOR_FFT_VDSP_H
#define _GABORATOR_FFT_VDSP_H

#include <assert.h>

#include <iterator>
#include <vector>

#include <memory.h>

#include <mach/mach.h>
#include <mach/task.h>
#include <mach/task_info.h>
#include <mach/vm_map.h>

#include <Accelerate/Accelerate.h>

namespace gaborator {

static inline int log2_int_exact(int n) {
    // n must be a power of two
    assert(n != 0 && ((n & (n >> 1)) == 0));
    int r = 0;
    for (;;) {
        n >>= 1;
        if (n == 0)
            break;
        r++;
    }
    return r;
}

template <>
struct fft<std::complex<float> *> {
    typedef std::complex<float> *I;
    typedef typename std::iterator_traits<I>::value_type C; // complex
    typedef typename C::value_type T; // float/double

    fft(unsigned int n_): n(n_), log2n(log2_int_exact(n)) {
        setup = vDSP_create_fftsetup(log2n, kFFTRadix2);
    }
    ~fft() {
        vDSP_destroy_fftsetup(setup);
    }

    unsigned int size() { return n; }

    // in-place
    void
    transform(I a) {
        DSPSplitComplex s;
        // XXX this result in disoptimal alignment
        s.realp = (float *) a;
        s.imagp = (float *) a + 1;
        vDSP_fft_zip(setup, &s, 2, log2n, kFFTDirection_Forward);
    }

    void
    itransform(I a) {
        DSPSplitComplex s;
        s.realp = (float *) a;
        s.imagp = (float *) a + 1;
        vDSP_fft_zip(setup, &s, 2, log2n, kFFTDirection_Inverse);
    }

    // out-of-place
    // XXX const
    void
    transform(I in, I out) {
        DSPSplitComplex si;
        si.realp = (float *) in;
        si.imagp = (float *) in + 1;
        DSPSplitComplex so;
        so.realp = (float *) out;
        so.imagp = (float *) out + 1;
        vDSP_fft_zop(setup,
                     &si, 2,
                     &so, 2,
                     log2n, kFFTDirection_Forward);
    }

    void
    itransform(I in, I out) {
        DSPSplitComplex si;
        si.realp = (float *) in;
        si.imagp = (float *) in + 1;
        DSPSplitComplex so;
        so.realp = (float *) out;
        so.imagp = (float *) out + 1;
        vDSP_fft_zop(setup,
                     &si, 2,
                     &so, 2,
                     log2n, kFFTDirection_Inverse);
    }

private:
    // Size of the transform
    unsigned int n;
    unsigned int log2n;
    FFTSetup setup;
};

// Support transforming std::vector<std::complex<float> >::iterator

template <>
struct fft<typename std::vector<std::complex<float> >::iterator>:
    public fft<std::complex<float> *>
{
    typedef fft<std::complex<float> *> base;
    typedef typename std::vector<std::complex<float> >::iterator I;
    fft(unsigned int n_): fft<std::complex<float> *>(n_) { }
    void
    transform(I a) {
        base::transform(&(*a));
    }
    void
    itransform(I a) {
        base::itransform(&(*a));
    }
    void
    transform(I in, I out) {
        base::transform(&(*in), &(*out));
    }
    void
    itransform(I in, I out) {
        base::itransform(&(*in), &(*out));
    }
};

} // Namespace

#endif
