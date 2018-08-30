//
// Rendering of spectrogram images
//
// Copyright (C) 2015-2018 Andreas Gustafsson.  This file is part of
// the Gaborator library source distribution.  See the file LICENSE at
// the top level of the distribution for license information.
//

#ifndef _GABORATOR_RENDER_H
#define _GABORATOR_RENDER_H

#include "gaborator/gaborator.h"
#include "gaborator/resample2.h"

namespace gaborator {


// Convert a floating-point linear brightness value in the range 0..1
// into an 8-bit pixel value, with clamping and (rough) gamma
// correction.  This nominally uses the sRGB gamma curve, but the
// current implementation cheats and uses a gamma of 2 because it can
// be calculated quickly using a square root.

template <class T>
unsigned int float2pixel_8bit(T val) {
    // Clamp before gamma correction so we don't take the square root
    // of a negative number; those can arise from bicubic
    // interpolation.  While we're at it, let's also skip the gamma
    // correction for small numbers that will round to zero anyway,
    // and especially denormals which could rigger GCC bug target/83240.
    static const T almost_zero = 1.0 / 65536;
    if (val < almost_zero)
        val = 0;
    if (val > 1)
        val = 1;
    return (unsigned int)(sqrtf(val) * 255.0f);
}


// Magnitude

template<class T>
struct complex_abs_fob {
    T operator()(const complex<T> &c) {
        return complex_abs(c);
    }
};


// A source object for resample2() that provides the absolute
// values of a row of spectrogram coordinates.

template <class T, class OI, class NORMF>
struct abs_row_source {
    typedef complex<T> C;

    typedef transform_output_iterator<NORMF, OI, C> abs_writer_t;
    abs_row_source(const analyzer<T> &frs_,
                   const sliced_coefs<C> &sc_,
                   int oct_, unsigned int obno_,
                   NORMF normf_):
        rs(frs_, sc_, oct_, obno_),
        normf(normf_)
    { }
    OI operator()(sample_index_t i0, sample_index_t i1, OI output) const {
        abs_writer_t abswriter(normf, output);
        abs_writer_t abswriter_end = rs(i0, i1, abswriter);
        return abswriter_end.output;
    }
    row_source<T, abs_writer_t> rs;
    NORMF normf;
};

// Helper class for abs_row_source specialization below

template <class C, class OI>
struct abs_writer_dest {
    abs_writer_dest(OI output_): output(output_) { }
    void process_existing_slice(C *bv, size_t len) {
        complex_magnitude(bv, output, len);
        output += len;
    }
    void process_missing_slice(size_t len) {
        for (size_t i = 0; i < len; i++)
            *output++ = 0;
    }
    OI output;
};

// Partial specialization of class abs_row_source for NORMF = complex_abs_fob,
// for vectorization.

template <class T, class OI>
struct abs_row_source<T, OI, struct complex_abs_fob<T> > {
    typedef complex<T> C;
    // Note unused last arg
    abs_row_source(const analyzer<T> &frs_,
                   const sliced_coefs<C> &sc_,
                   int oct_, unsigned int obno_,
                   complex_abs_fob<T>):
        slicer(frs_, sc_, oct_, obno_)
    { }
    OI operator()(coef_index_t i0, coef_index_t i1, OI output) const {
        abs_writer_dest<C, OI> dest(output);
        slicer(i0, i1, dest);
        return dest.output;
    }
    row_foreach_slice<T, abs_writer_dest<C, OI>, C> slicer;
};

// Render a single line (single frequency band), with scaling by
// powers of two in the horizontal (time) dimension, and filtering to
// avoid aliasing when minifying.

template <class OI, class T, class NORMF>
OI render_p2scale_line(const analyzer<T> &frs,
                       const coefs<T> &msc,
                       int gbno,
                       int64_t xorigin,
                       sample_index_t i0, sample_index_t i1, int e,
                       bool interpolate,
                       OI output,
                       NORMF normf)
{
    int oct;
    unsigned int obno; // Band number within octave
    bool clip = ! frs.bno_split(gbno, oct, obno, false);
    if (clip) {
        for (sample_index_t i = i0; i < i1; i++)
            *output++ = (T)0;
        return output;
    }
    abs_row_source<T, T *, NORMF>
        abs_rowsource(frs, msc.octaves[oct], oct, obno, normf);

    // Scale by the downsampling factor of the band
    int scale_exp = frs.band_scale_exp(oct, obno);
    output = resample2(abs_rowsource, xorigin,
                       i0, i1, e - scale_exp,
                       interpolate, output);
    return output;
}

// Render a two-dimensional image with scaling by powers of two in the
// horizontal direction only.  In the vertical direction, there is
// always a one-to-one correspondence between bands and pixels.
// yi0 and yi1 already have the yorigin applied, so there is no
// yorigin argument.

template <class OI, class T, class NORMF>
OI render_p2scale_noyscale(const analyzer<T> &frs,
                           const coefs<T> &msc,
                           int64_t xorigin,
                           int64_t xi0, int64_t xi1, int xe,
                           int64_t yi0, int64_t yi1,
                           bool interpolate,
                           OI output,
                           NORMF normf)
{
    assert(xi1 >= xi0);
    int w = xi1 - xi0;
    int gbno0 = yi0;
    int gbno1 = yi1;
    for (int gbno = gbno0; gbno < gbno1; gbno++) {
        int oct;
        unsigned int obno; // Band number within octave
        bool clip = ! frs.bno_split(gbno, oct, obno, false);
        if (clip) {
            for (int x = 0; x < w; x++)
                *output++ = (T)0;
        } else {
            output = render_p2scale_line(frs, msc, gbno, xorigin,
                                         xi0, xi1, xe,
                                         interpolate, output, normf);
        }
    }
    return output;
}

// Source data from a column of a row-major two-dimensional array.
// data points to the beginning of a row-major array with an x
// range of x0..x1 and an y range from y0..y1, and operator()
// returns data from column x (where x is within the range x0..x1).

template <class OI>
struct transverse_source {
    transverse_source(float *data_,
                      int64_t x0_, int64_t x1_, int64_t y0_, int64_t y1_,
                      int64_t x_):
        data(data_),
        x0(x0_), x1(x1_), y0(y0_), y1(y1_),
        x(x_),
        stride(x1 - x0)
    { }
    OI operator()(int64_t i0, int64_t i1, OI out) const {
        assert(x >= x0);
        assert(x <= x1);
        assert(i1 >= i0);
        assert(i0 >= y0);
        assert(i1 <= y1);
        float *p = data + (x - x0) + (i0 - y0) * stride;
        while (i0 != i1) {
            *out++ = *p;
            p += stride;
            ++i0;
        }
        return out;
    }
    float *data;
    int64_t x0, x1, y0, y1, x;
    size_t stride;
};

template <class I, class T>
struct stride_iterator: public std::iterator<std::forward_iterator_tag, T> {
    stride_iterator(I it_, size_t stride_): it(it_), stride(stride_) { }
    T& operator*() { return *it; }
    stride_iterator<I, T>& operator++() {
        it += stride;
        return *this;
    }
    stride_iterator operator++(int) {
        stride_iterator old = *this;
        it += stride;
        return old;
    }
    I it;
    size_t stride;
};

// Render a two-dimensional image with scaling by powers of two in
// both the horizontal (time) and vertical (frequency) directions.
// The output may be written through "output" out of order, so
// "output" must be a random access iterator.

// Note the default template argument for NORMF.  This is needed
// because the compiler won't deduce the type of NORMF from the
// default function argument "NORMF normf = complex_abs_fob<T>()"
// when the normf argument is omitted; it is considered a "non-deduced
// context", being "a template parameter used in the parameter type of
// a function parameter that has a default argument that is being used
// in the call for which argument deduction is being done".
// Unfortuantely, this work-around of providing a default template
// argument requires C++11.

template <class OI, class T, class NORMF = complex_abs_fob<T> >
void render_p2scale(const analyzer<T> &frs,
                    const coefs<T> &msc,
                    int64_t xorigin, int64_t yorigin,
                    int64_t xi0, int64_t xi1, int xe,
                    int64_t yi0, int64_t yi1, int ye,
                    OI output,
                    bool interpolate = true,
                    NORMF normf = complex_abs_fob<T>())
{
    // Construct a temporary float image of the right width,
    // but still needing scaling of the height.  Include
    // extra scanlines at the top and bottom for interpolation.

    // Find the image bounds in the spectrogram coordinate system,
    // including the interpolation margin.  The Y bounds are in
    // bands and are used both to determine what to render into the
    // temporary image and for short-circuiting; the X bounds are in
    // samples, and are only used for short-circuiting.
    int64_t ysi0, ysi1;
    resample2_support(yorigin, yi0, yi1, ye, ysi0, ysi1);
    int64_t xsi0, xsi1;
    resample2_support(xorigin, xi0, xi1, xe, xsi0, xsi1);

    // Short-circuiting: if the image to be rendered falls entirely
    // outside the data, just set it to zero instead of resampling down
    // (potentially) high-resolution zeros to the display resolution.
    // This makes a difference when zooming out by a large factor, for
    // example such that the entire spectrogram falls within a single
    // tile; that tile will necessarily be expensive to calculate, but
    // the other tiles need not be, and mustn't be if we are going to
    // keep the total amount of work bounded by O(L) with respect
    // to the signal length L regardless of zoom.
    int64_t cxi0, cxi1;
    frs.get_coef_bounds(msc, cxi0, cxi1);
    if (ysi1 < 0 || // Entirely above
        ysi0 >= frs.n_bands_total - 1 || // Entirely below
        xsi1 < cxi0 || // Entirely to the left
        xsi0 >= cxi1) // Entirely to the right
    {
        size_t n = (yi1 - yi0) * (xi1 - xi0);
        for (size_t i = 0; i < n; i++)
            output[i] = (T)0;
        return;
    }

    // Allocate buffer for temporary image resampled in the X
    // direction but not yet in the Y direction
    size_t n_pixels = (ysi1 - ysi0) * (xi1 - xi0);
    pod_vector<float> render_data(n_pixels);

    // Render data resampled in the X direction
    float *p = render_data.data();
    render_p2scale_noyscale(frs, msc, xorigin, xi0, xi1, xe,
                            ysi0, ysi1, interpolate, p, normf);

    // Resample in the Y direction
    for (int64_t xi = xi0; xi < xi1; xi++) {
        transverse_source<OI> src(render_data.data(),
                                  xi0, xi1, ysi0, ysi1,
                                  xi);
        stride_iterator<OI, float> dest(output + (xi - xi0), (xi1 - xi0));
        resample2(src, yorigin, yi0, yi1, ye, interpolate, dest);
    }
}


} // namespace

#endif
