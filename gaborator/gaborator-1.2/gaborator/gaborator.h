//
// Constant Q spectrum analysis and resynthesis
//
// Copyright (C) 2015-2018 Andreas Gustafsson.  This file is part of
// the Gaborator library source distribution.  See the file LICENSE at
// the top level of the distribution for license information.
//

#ifndef _GABORATOR_GABORATOR_H
#define _GABORATOR_GABORATOR_H

#define __STDC_LIMIT_MACROS

#include <assert.h>
#include <errno.h>
#include <inttypes.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>

#include <algorithm>
#include <complex>
#include <limits>
#include <map>
#include <typeinfo>

#include "gaborator/fft.h"
#include "gaborator/gaussian.h"
#include "gaborator/pod_vector.h"
#include "gaborator/pool.h"
#include "gaborator/ref.h"
#include "gaborator/vector_math.h"


namespace gaborator {

using std::complex;

// An integer identifying an audio sample
typedef int64_t sample_index_t;

// An integer identifying a coefficient
typedef int64_t coef_index_t;

// An integer identifying a slice
typedef int64_t slice_index_t;

// See https://tauday.com/tau-manifesto
static const double tau = 2.0 * M_PI;

// Round up to next higher or equal power of 2

inline int
next_power_of_two(int x) {
    --x;
    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
    return x + 1;
}

// Determine if x is a power of two.
// Note that this considers 0 to be a power of two.

static inline bool
is_power_of_two(unsigned int x) {
    return (x & (x - 1)) == 0;
}

// Given a power of two v, determine log2(v)
// https://graphics.stanford.edu/~seander/bithacks.html#DetermineIfPowerOf2

static inline unsigned int whichp2(unsigned int v) {
    assert(is_power_of_two(v));
    unsigned int r = (v & 0xAAAAAAAA) != 0;
    r |= ((v & 0xCCCCCCCC) != 0) << 1;
    r |= ((v & 0xF0F0F0F0) != 0) << 2;
    r |= ((v & 0xFF00FF00) != 0) << 3;
    r |= ((v & 0xFFFF0000) != 0) << 4;
    return r;
}

// Floor division: return the integer part of a / b
// rounded down (not towards zero).  For positive b only.

inline int64_t floor_div(int64_t a, int64_t b) {
    assert(b > 0);
    if (a >= 0)
        return a / b;
    else
        return (a - b + 1) / b;
}

// Floating point modulus, the remainder r of a / b
// satisfying 0 <= r < b even for negative a.
// For positive b only.

static inline double
sane_fmod(double a, double b) {
    assert(b > 0);
    double m = fmod(a, b);
    if (m < 0)
        m += b;
    return m;
}

// Do an arithmetic left shift of a 64-bit signed integer.  This is
// what a << b ought to do, but according to the C++11 draft (n3337),
// section 5.8, that invokes undefined behavior when a is negative.
// GCC is actually smart enough to optimize this into a single shlq
// instruction.
//
// No corresponding kludge is needed for right shifts, because a right
// shift of a negative signed integer is implementation-defined, not
// undefined, and we trust implementations to define it sanely.

static inline int64_t
shift_left(int64_t a, unsigned int b) {
    if (a < 0)
        return -(((uint64_t) -a) << b);
    else
        return (((uint64_t) a) << b);
}

// Convert between complex types

template<class T, class U>
complex<T> c2c(complex<U> c) { return complex<T>(c.real(), c.imag()); }

// Convert a sequence of complex values to real

template <class I, class O>
O complex2real(I b, I e, O o) {
    while (b != e) {
        *o++ = (*b++).real();
    }
    return o;
}

// Test a sequence for being all zero

template <class I>
bool is_zero(I b, I e) {
    while (b != e) {
        if (*b)
            return false;
        ++b;
    }
    return true;
}

// A vector-like object that allows arbitrary integer indices
// (positive or negative, but excluding the largest possible integer)
// and automatically resizes the storage.  Uses storage proportional
// to the difference between the smallest and largest index value (for
// example, if indices range from -102 to -100 (inclusive), memory use
// is on the order of 3 elements).
//
// T is the element type
// I is the integer index type

template <class T, class I = int>
struct range_vector {
    range_vector():
        lower(std::numeric_limits<I>::max()),
        upper(std::numeric_limits<I>::min())
    { }
private:
    T *unchecked_get(I i) {
        return &v[i & ((I)v.size() - 1)];
    }
    const T *unchecked_get(I i) const {
        return &v[i & ((I)v.size() - 1)];
    }

public:
    // Note: Pointer returned becomes invalid when range_vector
    // is changed
    T *
    get(I i, bool create) {
        if (! has_index(i)) {
            if (create)
                extend(i);
            else
                return 0;
        }
        return unchecked_get(i);
    }

    // Get a pointer to an existing element, or null if out of range
    const T *
    get(I i) const {
        if (! has_index(i))
            return 0;
        return unchecked_get(i);
    }

    // Note: Reference returned becomes invalid when range_vector
    // is changed
    T &
    get_or_create(I i) {
        return *get(i, true);
    }

    // Get a reference to the element at index i, which must be valid
    T &
    get_existing(I i) {
        assert(has_index(i));
        return *unchecked_get(i);
    }

    // Const version of the above
    const T &
    get_existing(I i) const {
        assert(has_index(i));
        return *unchecked_get(i);
    }

private:
    void extend(I i) {
        I new_lower = lower;
        I new_upper = upper;
        if (i < lower)
            new_lower = i;
        if (i + 1 > upper)
            new_upper = i + 1;
        I old_size = v.size();
        I new_need = new_upper - new_lower;
        if (new_need > old_size) {
            if (old_size == 0) {
                v.resize(1);
            } else {
                I new_size = old_size;
                while (new_size < new_need)
                    new_size *= 2;
                v.resize(new_size);
                if (old_size) {
                    for (I j = lower; j < upper; j++) {
                        I jo = j & (old_size - 1);
                        I jn = j & (new_size - 1);
                        if (jo != jn)
                            std::swap(v[jo], v[jn]);
                    }
                }
            }
        }
        lower = new_lower;
        upper = new_upper;
    }

public:
    // Erase the elements whose index is less than "limit"
    void erase_before(I limit) {
        I i = lower;
        for (;i < upper && i < limit; i++)
            *unchecked_get(i) = T();
        lower = i;
    }

    I begin_index() const { return lower; }
    I end_index() const { return upper; }
    bool empty() const { return lower >= upper; }
    bool has_index(ssize_t i) const { return i >= lower && i < upper; }

private:
    std::vector<T> v;
    I lower, upper;
};

// Calculate the size of the alias-free part (the "filet")
// of a signal slice of size "fftsize"

static inline unsigned int filet_part(unsigned int fftsize) {
    return fftsize >> 1;
}

// Calculate the size of the padding (the "fat") at each
// end of a signal slice of size "fftsize"

static inline unsigned fat_part(unsigned int fftsize) {
    return fftsize >> 2;
}

// Frequency band parameters shared between octaves

template<class T>
struct band_params: public refcounted {
    typedef complex<T> C;
    bool dc; // True iff this is the DC band
    unsigned int sftsize; // Size of "short FFT" spanning the band
    unsigned int sftsize_log2; // log2(sftsize)
    fft<C *> *sft; // Fourier transform for windows, of size sftsize
    std::vector<T> kernel; // Frequency-domain filter kernel
    std::vector<T> dual_kernel; // Dual of the above
    pod_vector<C> shift_kernel; // Complex exponential for fractional frequency compensation
    pod_vector<C> shift_kernel_conj; // Conjugate of the above
    int fq_offset_int; // Frequency offset in bins (big-fft bin of left window edge)
    double ff; // Center (bp) or corner (lp) frequency in units of the sampling frequency
    double center; // Center frequency in units of FFT bins
    int icenter; // Center frequency rounded to nearest integer FFT bin
    double ffsd; // Standard deviation of the bandpass Gaussian, as fractional frequency
    float anl_support;
    float syn_support;
};

// Downsampling parameters.  These have some similarity to band
// parameters, but only some.  For example, these may use a real
// rather than complex FFT for the "short FFT".

template <class T>
struct downsampling_params {
    typedef complex<T> C;
    unsigned int sftsize;
    std::vector<T> kernel; // Frequency-domain filter kernel
    std::vector<T> dual_kernel;
#if GABORATOR_USE_REAL_FFT
    rfft<C *> *rsft;
#else
    fft<C *> *sft;
#endif
    double time_support; // Filter time domain support, each side
};

// Forward declarations
template<class T> struct analyzer;
template<class T> struct zone;
template<class T, class C = complex<T> >  struct coefs;
template<class C> struct sliced_coefs;
template<class T, class CT = complex<T> > struct rw_row_view;

// Abstract class for tracking changes to a coefficient set.
// This may be used for updating a resolution pyramid of
// magnitude data.

template<class T>
struct shadow {
    virtual ~shadow() { }
    virtual void update(const coefs<complex<T> > &msc, sample_index_t i0, sample_index_t i1) = 0;
    // Deprecated
    virtual void update(int oct, int ze, const sliced_coefs<complex<T> > &sc, slice_index_t sli0, slice_index_t sli1) = 0;
};

// Coefficient metadata.  This contains information describing a
// "coefs" structure that is common to many instances and should not
// be duplicated in each one.  In practice, there will be two
// coefs_meta structures per zone per analyzer, one for unpadded coefs
// and one for padded coefs.

struct coefs_meta {
    typedef std::vector<unsigned int> shape_vector;
    void init(const shape_vector &shape) {
        band_offsets.resize(shape.size());
        unsigned int offset = 0;
        for (unsigned int i = 0; i < shape.size(); i++) {
            band_offsets[i] = offset;
            offset += shape[i];
        }
        total_size = offset;
    }
    // Offset of the beginning of each band in the data array
    std::vector<unsigned int> band_offsets;
    // Total size of data array (in elements, not bytes)
    unsigned int total_size;
};

// Coefficients of a single octave for a single input signal slice.
// These are used both as part of the final sliced coefficients
// and for temporary padded coefficients during analysis/synthesis.
// C is the coefficient type, typically complex<float> but can also
// be e.g. unsigned int to store cluster numbers, or float to store
// magnitudes.

template<class C>
struct oct_coefs: public refcounted {
    oct_coefs(const coefs_meta &meta_, bool clear_ = true):
        meta(meta_),
        data(meta.total_size),
        bands(*this)
    {
        if (clear_)
            clear();
    }
    uint64_t estimate_memory_usage() const {
        return meta.total_size * sizeof(C) + sizeof(*this);
    }
    void clear() {
        memset(data.data(), 0, data.size() * sizeof(C));
    }
    // Deep copy
    oct_coefs &operator=(const oct_coefs &rhs) {
        assert(data.size() == rhs.data.size());
        memcpy(data.data(), rhs.data.data(), data.size() * sizeof(C));
        return *this;
    }

    const coefs_meta &meta;

    // The data for all the bands are allocated together
    // as a single vector to reduce the number of allocations
    pod_vector<C> data;
    // Vector-like collection of pointers into "data", one for each band
    struct band_array {
        band_array(oct_coefs &outer_): outer(outer_) { }
        C *operator[](size_t i) const {
            return outer.data.data() + outer.meta.band_offsets[i];
        }
        size_t size() const { return outer.meta.band_offsets.size(); }
        oct_coefs &outer;
    } bands;
private:
    oct_coefs(const oct_coefs &);
};


// Add the oct_coefs "b" to the oct_coefs "a"

template <class T, class C>
void add(zone<T> &z, oct_coefs<C> &a, const oct_coefs<C> &b) {
    unsigned int n_bands = a.bands.size();
    assert(n_bands == b.bands.size());
    for (unsigned int obno = 0; obno < n_bands; obno++) {
        unsigned int len = filet_part(z.bandparams[obno]->sftsize);
        complex<T> *band_a = a.bands[obno];
        complex<T> *band_b = b.bands[obno];
        for (unsigned int j = 0; j < len; j++) {
            band_a[j] += band_b[j];
        }
    }
}

// Temporary coefficients used during analysis/synthesis.
// These include padding (aka fat).

template<class C>
struct padded_coefs: public oct_coefs<C> {
    padded_coefs(const coefs_meta &meta_, bool clear_ = true):
        oct_coefs<C>(meta_, clear_) { }
};

// Sliced coefficients.  These cover an arbitrary time range, but only
// a single octave.  Template argument is as for struct oct_coefs.

template<class C>
struct sliced_coefs {
    typedef range_vector<ref<oct_coefs<C> >, slice_index_t> slices_t;
    slices_t slices;
    uint64_t estimate_memory_usage() const {
        unsigned int n = 0;
        size_t size_each = 0;
        for (slice_index_t sl = slices.begin_index(); sl < slices.end_index(); sl++) {
            const ref<oct_coefs<C> > &t = slices.get_existing(sl);
            if (t) {
                if (! size_each)
                    size_each = t->estimate_memory_usage();
                n++;
            }
        }
        return n * size_each;
    }
};

// Multirate sliced coefficients.  These cover an arbitrary time
// range and the full frequency space (all octaves).
// Template arguments are as for struct coef.
// Note default for template argument C defined in forward declaration.

template<class T, class C>
struct coefs {
    coefs(const analyzer<T> &anl_, shadow<T> *shadow_ = 0):
        octaves(anl_.n_octaves), shadow0(shadow_)
    { }
    uint64_t estimate_memory_usage() const {
        uint64_t s = 0;
        for (unsigned int oct = 0; oct < octaves.size(); oct++)
            s += octaves[oct].estimate_memory_usage();
        return s;
    }
    std::vector<sliced_coefs<C> > octaves;
    shadow<T> *shadow0;
};

// Perform an fftshift of the range between iterators a and b.
// Not optimized - not for use in inner loops.

template <class I>
void fftshift(I b, I e) {
    int len = e - b;
    assert(len % 2 == 0);
    for (int i = 0; i < len / 2; i++)
        std::swap(*(b + i), *(b + len / 2 + i));
}

// Given an unsigned index i into an FFT of a power-of-two size
// "size", return the corresponding signed index, ranging from -size/2
// to size/2-1.  This is equivalent to sign extension of an integer of
// log2(size) bits; see Hacker's Delight, page 18.  This can be used
// to convert an FFT index into a frequency (positive or negative;
// note that Nyquist is considered negatitive = -fs/2).

static inline int signed_index(unsigned int i, unsigned int size) {
    unsigned int t = size >> 1;
    return (i ^ t) - t;
}

// Construct a frequency-domain lowpass filter whose response is the
// convolution of a rectangle and a gaussian.  The cutoff freqency is
// ff_cutoff (a fractional frequency), and the standard deviation of
// the gaussian is ff_sd.  The returned filter covers the full
// frequency range from 0 to fs (with negative frequencies at the end,
// the usual convention for FFT spectra).
//
// When center=true, construct a time-domain window instead,
// passing the center of the time-domain signal.
//
// The result is stored between iterators b and e, which must have a
// real value_type.

template <class I>
inline void gaussian_windowed_lowpass(double ff_cutoff, double ff_sd,
                                      I b, I e, bool center = false)
{
    size_t len = e - b;
    double inv_len = 1.0 / len;
    for (I it = b; it != e; ++it) {
        size_t i = it - b;
        double thisff;
        if (center)
            // Symmetric around center
            thisff = std::abs(i - (len * 0.5)) * inv_len;
        else
            // Symmetric around zero
            thisff = (i > len / 2 ? len - i : i) * inv_len;
        double x = thisff - ff_cutoff;
        double v = gaussian_edge(ff_sd, -x);
        *it = v;
    }
}

// A set of octaves having identical parameters form a "zone",
// and their shared parameters are stored in a "struct zone".

template <class T>
struct zone: public refcounted {
    zone(): n_bands(0) { }
    ~zone() { }
    unsigned int n_bands; // Total number of bands, including DC band if lowest octave
    // Band parameters by increasing frequency; DC band is index 0 if present
    std::vector<ref<band_params<T> > > bandparams;
    // Pseudo-bands mimicing the response of bands in the
    // neighboring octaves, used for calculating the duals only
    std::vector<ref<band_params<T> > > mock_bandparams;
    pod_vector<T> power;
    pod_vector<T> power_nodc;
    coefs_meta cmeta[2]; // Unpadded and padded
};

template <class T>
struct octave {
    zone<T> *z;
    unsigned int n_bands_above; // Total number of bands in higher octaves
};


// Helper function for pushing parameters onto the vectors in struct zone

template <class T>
void push(std::vector<ref<band_params<T> > > &v, band_params<T> *p) {
    v.push_back(ref<band_params<T> >(p));
}

// A set of spectum analysis parameters

struct parameters {
    parameters(unsigned int bands_per_octave_, double ff_min_,
               double ff_ref_ = 1.0,
               double overlap_ = 0.7,
               double max_error_ = 1e-5):
        bands_per_octave(bands_per_octave_),
        ff_min(ff_min_),
        ff_ref(ff_ref_),
        overlap(overlap_),
        max_error(max_error_)
    { }
    // Provide an operator< so that we can create a set or map of parameters
    bool operator<(const parameters &b) const {
#define GABORATOR_COMPARE_LESS(member) do { \
        if (member < b.member) \
            return true; \
        if (member > b.member) \
            return false; \
        } while(0)
        GABORATOR_COMPARE_LESS(bands_per_octave);
        GABORATOR_COMPARE_LESS(ff_min);
        GABORATOR_COMPARE_LESS(ff_ref);
        GABORATOR_COMPARE_LESS(overlap);
        GABORATOR_COMPARE_LESS(max_error);
#undef GABORATOR_COMPARE_LESS
        // Equal
        return false;
    }
    bool operator==(const parameters &b) const {
        return !((*this < b) || (b < *this));
    }
    template<class T> friend class analyzer;
    unsigned int bands_per_octave;
    double ff_min;
    double ff_ref;
    double overlap;
    double max_error;
};

// Index of first slice affected by sample at t0

// fft number i covers the sample range
// t = (i * filetsize .. i * filetsize + (fftsize - 1))
// t >= i * filetsize and t < i * filetsize + fftsize
// A sample at t affects ffts i where
//   i <= t / filetsize and
//   i > (t - fftsize) / filetsize
// the filet of fft number i covers the sample range
// (fat + (i * filetsize) .. fat + (i * filetsize) + (filetsize - 1))

static inline slice_index_t affected_slice_0(sample_index_t t0, unsigned int fftsize) {
    return floor_div(t0 - fftsize, filet_part(fftsize)) + 1;
}

// Index of first slice not affected by samples before t1

static inline slice_index_t affected_slice_1(sample_index_t t1, unsigned int fftsize) {
    return floor_div(t1 - 1, filet_part(fftsize)) + 1;
}


// Multiply a vector by a scalar, in-place.
// Used only at the setup stage, so performance is not critical.

template <class V, class S>
void scale_vector(V &v, S s) {
    for (size_t i = 0; i < v.size(); i++)
        v[i] *= s;
}

// Fill the buffer at dst, of length dstlen, with data from src where
// available, otherwise with zeroes.  The data in src covers dst indices
// from i0 (inclusive) to i1 (exclusive).

template <class T>
void copy_overlapping_zerofill(T *dst, unsigned int dstlen, const T *src,
                               int64_t src_i0, int64_t src_i1)
{
    int64_t overlap_begin = std::max((int64_t) 0, src_i0);
    int64_t overlap_end = std::min((int64_t)dstlen, src_i1);
    if (overlap_end <= overlap_begin) {
        // No overlap
        std::fill(dst, dst + dstlen, 0);
    } else {
        // Overlap
        if (overlap_begin != 0)
            std::fill(dst, dst + overlap_begin, 0);
        std::copy(src + overlap_begin - src_i0, src + overlap_end - src_i0, dst + overlap_begin);
        if (overlap_end != dstlen)
            std::fill(dst + overlap_end, dst + dstlen, 0);
    }
}

// Given a set of FFT coefficients "coefs" of a real
// sequence, where only positive-frequency coefficients
// (including DC and Nyquist) are valid, return the
// coefficient for an arbitrary frequency index "i"
// which may correspond to a negative frequency, or
// even an alias outside the range (0..fftsize-1).

template<class T>
complex<T> get_real_spectrum_coef(complex<T> *coefs, int i, unsigned int fftsize) {
    i &= fftsize - 1;
    // Note that this is >, not >=, becase fs/2 is considered nonnegative
    bool neg_fq = (i > (int)(fftsize >> 1));
    if (neg_fq) {
        i = fftsize - i;
    }
    complex<T> c = coefs[i];
    if (neg_fq) {
        c = conj(c);
    }
    return c;
}

template<class T>
struct analyzer: public refcounted {
    typedef complex<T> C;
    analyzer(const parameters &params_):
        params(params_)
    {
        // Sanity check
        assert(params.ff_min < 0.5);

        // The frequency increases by this factor from one band
        // to the next
        band_spacing_log2 = 1.0 / params.bands_per_octave;
        band_spacing = exp2(band_spacing_log2);

        // The tuning adjustment, as a log2ff.  This is a number between
        // 0 and band_spacing_log2, corresponding to a frequency at
        // or slightly above the sampling frequency where a band
        // center would fall if they actually went that high.
        // Tuning is done by increasing the center frequencies of
        // all bands by this amount relative to the untuned case
        // where one band would fall on fs exactly.
        tuning_log2ff = sane_fmod(log2(params.ff_ref), band_spacing_log2);

        // Calculate the frequency of band0, the "fs/8 band", the
        // lowest band in each octave except possibly the lowest octave.
        // Its frequency will fs/8, or a fraction of the band spacing higher
        // due to tuning. The magic -3 derives from fs/8 as log2(1/8).
        band0_log2ff = -3 + tuning_log2ff;

        // Calculate the total number of bands needed so that
        // the lowest band has a frequency <= params.ff_min.
        // end_log2ff = the log2ff of the band after the last (just past fs/2)
        double end_log2ff = tuning_log2ff - 1;
        n_bandpass_bands_total =
            (unsigned int)ceil((end_log2ff - log2(params.ff_min)) / band_spacing_log2);

        // Calculate the kernel support needed for band0,
        // and size the FFT accordingly.  This duplicates some
        // code at the beginning of make_band().
        double band0_ff = exp2(band0_log2ff);
        double band0_time_sd = time_sd(band0_ff);

        double band0_time_support = gaussian_support(band0_time_sd, params.max_error);
        double band0_time_synthesis_support = band0_time_support * synthesis_support_multiplier();

        fftsize_log2 = 3;
        fftsize = 1 << fftsize_log2;
        while (band0_time_synthesis_support > fat_part(fftsize)) {
            fftsize_log2++;
            fftsize <<= 1;
        }
        init_with_fftsize();
    }

    void init_with_fftsize() {
        // Clear everything, including things that may have been set
        // on a previous try.
        sftsize_max = 0;
        octaves.clear();
        zones.clear();

        inv_fftsize_double = 1.0 / fftsize;
        inv_fftsize_t = (T) inv_fftsize_double;

#if GABORATOR_USE_REAL_FFT
        rft = pool<rfft<C *>, int>::shared.get(fftsize);
#else
        ft = pool<fft<C *>, int>::shared.get(fftsize);
#endif

        // Band number starting at 0 close to fs/2 and increasing
        // with decreasing frequency
        int tbno = 0;
        // Band number of the fs/8 band in the current octave
        int base = 0;
        int zno = 0;
        // Loop over the octaves, from high to low frequencies,
        // creating new zones where needed
        for (;;) {
            int max_bands_this_octave = (zno == 0) ?
                params.bands_per_octave * 2 : params.bands_per_octave;
            base += max_bands_this_octave;
            int bands_remaining = n_bandpass_bands_total - tbno;
            int bands_this_octave = std::min(max_bands_this_octave, bands_remaining);
            int bands_below = bands_remaining - bands_this_octave;
            bool dc_zone = (bands_below == 0);
            bool dc_adjacent_zone = (bands_below < (int)params.bands_per_octave);
            if (zno < 2 || dc_zone || dc_adjacent_zone) {
                make_zone(zno, base - (tbno + bands_this_octave), base - tbno, dc_zone, bands_below);
                zno++;
            }
            octaves.push_back(octave<T>());
            octaves.back().z = zones[zno - 1].get();
            octaves.back().n_bands_above = tbno;
            tbno += bands_this_octave;
            if (dc_zone)
                break;
        }
        n_octaves = octaves.size();

        // Verify the total number of bands
        n_bands_total = 0;
        for (unsigned int i = 0; i < n_octaves; i++) {
            n_bands_total += octaves[i].z->n_bands;
        }
        assert(n_bandpass_bands_total + 1 == n_bands_total);

        // Set up the downsampling parameters in dsparams.

        // Downsampling is always by a factor of two.
        // dsparams.sftsize is the size of the FFT used to go back to
        // the time domain after discarding the top half of the
        // spectrum.
        dsparams.sftsize = fftsize >> 1;
        dsparams.kernel.resize(dsparams.sftsize);
        dsparams.dual_kernel.resize(dsparams.sftsize);

        // In terms of the post-downsampling sampling frequency,
        // the downsampling lowpass filter transition band
        // starts at the top edge of the highest band,
        // and ranges to fs/2.   We use the worst-case center
        // frequency for the highest band, fs/4, and add
        // the support.
        double f0 = 0.25 + ff_sd(0.25);
        double f1 = 0.5;

        // The cutoff frequency is in the center of the transition band
        double ff = (f0 + f1) * 0.5;
        double support = (f1 - f0) * 0.5;
        double ff_sd = gaussian_support_inv(support, params.max_error);
        // Fudge: make the lowpass a bit steeper (and correspondingly
        // wider in the time dimension).  This seems to improve the S/N
        // a bit.
        ff_sd *= 0.7813026596806952;

        // Calculate and save the time-domain support of the
        // downsampling lowpass filter for use in analyze_sliced().
        double time_sd = sd_f2t(ff_sd);
        dsparams.time_support = gaussian_support(time_sd, params.max_error);

        // Use the convolution of a rectagle and a gaussian.
        // A piecewise function composed from two half-gaussians
        // joined by a horizontal y=1 segment is not quite smooth
        // enough.
        gaussian_windowed_lowpass(ff, ff_sd, dsparams.kernel.begin(), dsparams.kernel.end());
        // Put the passband in the middle
        fftshift(dsparams.kernel.begin(), dsparams.kernel.end());
        // The dual_kernel field of the downsampling pseudo-band holds
        // the upsampling filter, identical to the downsampling filter
        // except for amplitude scaling.
        std::copy(dsparams.kernel.begin(), dsparams.kernel.end(),
                  dsparams.dual_kernel.begin());
        // Prescale the downsampling filter
        scale_vector(dsparams.kernel, inv_fftsize_double);
        // Prescale the upsampling filter
        scale_vector(dsparams.dual_kernel, 1.0 / dsparams.sftsize);
#if GABORATOR_USE_REAL_FFT
        dsparams.rsft = pool<rfft<C *>, int>::shared.get(dsparams.sftsize);
#else
        dsparams.sft = pool<fft<C *>, int>::shared.get(dsparams.sftsize);
#endif
        top_band_log2ff = band0_log2ff + band_spacing_log2 * (2 * params.bands_per_octave - 1);

        ffref_gbno = (int)rint((top_band_log2ff - log2(params.ff_ref)) / band_spacing_log2);
    }

    void make_zone(unsigned int zno, int bb, int be, bool dc_zone, int bands_below) {
        assert(zones.size() == zno);
        zone<T> *z = new zone<T>();
        zones.push_back(ref<zone<T> >(z));

        // The maximum number of mock bands to insert between the DC
        // band and the real bands.
        int n_mock_bands = 6;

        if (dc_zone) {
            // This zone goes all the way to DC
            band_params<T> *dc_band = make_band(bb - 1, true);
            push(z->bandparams, dc_band);
        } else {
            // There are other zones below this; add mock bands
            // to simulate them for purposes of calculating the dual
            if (n_mock_bands > bands_below)
                n_mock_bands = bands_below;
            // Mock DC band
            push(z->mock_bandparams,
                 make_band(bb - 1 - n_mock_bands, true));
            // Mock bandpass bands
            for (int i = bb - 1; i > bb - 1 - n_mock_bands; i--)
                push(z->mock_bandparams, make_band(i, false));
        }

        // The actual bandpass bands of this zone
        for (int i = bb; i < be; i++)
            push(z->bandparams, make_band(i, false));

        // If there are other zones above this, add mock bands
        // to simulate them for purposes of calculating the dual
        for (int i = be; i < (int)params.bands_per_octave * 2; i++)
            push(z->mock_bandparams, make_band(i, false));

        z->n_bands = z->bandparams.size();

        // Accumulate window power for calculating dual

        z->power.resize(fftsize);
        std::fill(z->power.data(), z->power.data() + fftsize, 0);
        for (unsigned int obno = 0; obno < z->bandparams.size(); obno++)
            accumulate_power(z->bandparams[obno].get(), z->power.data());
        for (unsigned int obno = 0; obno < z->mock_bandparams.size(); obno++)
            accumulate_power(z->mock_bandparams[obno].get(), z->power.data());

        // Calculate complex exponentials for non-integer center frequency adjustment
        for (unsigned int obno = 0; obno < z->bandparams.size(); obno++) {
            band_params<T> *bp = z->bandparams[obno].get();
            for (unsigned int i = 0; i < bp->sftsize; i++) {
                unsigned int ii = i;
                double arg = tau * ((double)ii / bp->sftsize) * -(bp->center - bp->icenter);
                C t(cos(arg), sin(arg));
                bp->shift_kernel[i] = t;
                bp->shift_kernel_conj[i] = conj(t);
            }
        }

        // Calculate duals
        for (unsigned int obno = 0; obno < z->bandparams.size(); obno++) {
            band_params<T> *bp = z->bandparams[obno].get();
            for (unsigned int i = 0; i < bp->sftsize; i++) {
                // ii = large-FFT bin number
                int ii = i + bp->fq_offset_int;
                bp->dual_kernel[i] = bp->kernel[i] / z->power[ii & (fftsize - 1)];
            }
        }

        // Initialize coefficient metadata for unpadded and padded coefficients
        for (int padded = 0; padded < 2; padded++) {
            typename coefs_meta::shape_vector shape(z->bandparams.size());
            for (unsigned int i = 0; i < z->bandparams.size(); i++) {
                unsigned int len = z->bandparams[i]->sftsize;
                if (! padded)
                    len >>= 1;
                shape[i] = len;
            }
            z->cmeta[padded].init(shape);
        }

        // Free the mock band parameters
        z->mock_bandparams.clear();
        z->power.clear();
    }

    // Add the power of the kernel in "*bp" to "power"
    void
    accumulate_power(band_params<T> *bp, T *power) {
        for (unsigned int i = 0; i < bp->sftsize; i++) {
            // ii = large-FFT bin number
            unsigned int ii = (i + bp->fq_offset_int) & (fftsize - 1);
            assert(ii >= 0 && ii < fftsize);
            T y = bp->kernel[i];
            T p = y * y;
            power[ii] += p;
            if (! bp->dc) {
                unsigned int ni = fftsize - ii;
                ni &= fftsize - 1; // XXX is this correct?
                assert(ni < fftsize);
                power[ni] += p;
            }
        }
    }

    // Given a fractional frequency, return the standard deviation
    // of the frequency-domain window as a fractional frequency
    double ff_sd(double ff) const { return params.overlap * (band_spacing - 1) * ff; }

    // Given a fractional frequency, return the standard deviation
    // of the time-domain window in samples.
    //
    // ff_sd = 1.0 / (tau * t_sd)
    // per http://users.ece.gatech.edu/mrichard/Gaussian%20FT%20and%20random%20process.pdf
    // and python test program gaussian-overlap.py
    // => (tau * t_sd) * ff_sd = 1.0
    // => t_sd = 1.0 / (tau * f_sd)
    double time_sd(double ff) const { return 1.0 / (tau * ff_sd(ff)); }

    // Defining Q as the frequency divided by the half-power bandwidth,
    // we get
    //
    // norm_gaussian(sd, hbw) = sqrt(2)
    //
    // (%i1) e1: exp(-(hbw * hbw) / (2 * sd * sd)) = 1 / sqrt(2);
    // (%i2) solve(e1, hbw);
    // (%o2) [hbw = - sqrt(log(2)) sd, hbw = sqrt(log(2)) sd]
    //
    // Q = ff / (2 * sqrt(log(2)) * ff_sd(ff))
    //   = 1.0 / ((2 * sqrt(log(2)) * (params.overlap * band_spacing - 1)))
    double q() const {
        return 1.0 / ((2 * sqrt(log(2)) * params.overlap * (band_spacing - 1)));
    }

    // Find the worst-case time support of the analysis filters, i.e.,
    // the largest distance in time between a signal sample and a
    // coefficient affected by that sample.
    double analysis_support() const {
        int gbno = n_bands_total - 2; // Last before DC
        int oct;
        unsigned int obno;
        bool valid = bno_split(gbno, oct, obno, false);
        assert(valid);
        double ff = band_ff(oct, obno);
        double timesd = time_sd(ff);
        double time_support = gaussian_support(timesd, params.max_error);
        return time_support;
    }

    // Ditto for the resynthesis filters.
    double synthesis_support() const {
        return analysis_support() * synthesis_support_multiplier();
    }

    // Empirical formula for synthesis support multiplier
    double synthesis_support_multiplier() const {
        return 1.8 / params.overlap;
    }

    // Return the fractional frequency of relative band number
    // "rbno" (relative to the sampling frequency of its octave).
    double rbno_ff(double rbno) const {
        double log2ff = band0_log2ff + rbno * band_spacing_log2;
        return exp2(log2ff);
    }

    // Return the fractional frequency of band number "obno"
    // in octave "oct", scaling according to the octave
    double band_ff(int oct, int obno) const {
        int gbno = bno_merge(oct, obno);
        return exp2(tuning_log2ff - 1 - (gbno + 1) * band_spacing_log2);
    }

    // Return the coefficient index (the time in terms of coefficient
    // subsamples) of the first cofficient of slice "sli" of band
    // "obno" in octave "oct"
    coef_index_t coef_time(slice_index_t sli, int oct, int obno) const {
        int sftsize = octaves[oct].z->bandparams[obno]->sftsize;
        return fat_part(sftsize) + sli * filet_part(sftsize);
    }

    // Return the sample index (the time in terms of samples) time of
    // coefficient "i" in slice "sli" of band "obno" in octave "oct"
    sample_index_t sample_time(slice_index_t sli, int i, int oct, int obno) const {
        coef_index_t sst = coef_time(sli, oct, obno) + i;
        return shift_left(sst, band_scale_exp(oct, obno));
    }

    // Calculate band parameters for a single band.
    //
    // rbno is a "relative band number" indicating a frequency
    // within its octave: it is 0 for the fs/8 band and increases
    // with frequency.  It is not always the same as the index
    // into its octaves' bandparams[] or the coefficient bands
    // (those would be called obno, not rbno).
    //
    // If dc is true, this is the DC band, and rbno indicates
    // the cutoff frequency; it is one less than the rbno of
    // the lowest-frequency bandpass band.

    band_params<T> *
    make_band(double rbno, bool dc) {
        if (dc)
            // Make the actual DC band cutoff frequency a bit higher,
            // by an empirically chosen fraction of a band, to reduce
            // power fluctuations.
            rbno += 0.8750526596806952;

        // For bandpass bands, the center frequency, or
        // for the DC band, the lowpass cutoff frequency,
        // as a fractional frequency.
        double ff = rbno_ff(rbno);

        // Standard deviation of the bandpass Gaussian,
        // as a fractional frequency
        double ffsd = ff_sd(ff);
        // The support of the Gaussian, i.e., the smallest standard
        // deviation at which it can be truncated on each side
        // without the error exceeding our part of the error budget,
        // which is some fraction of params.max.  Note
        // that this is one-sided; the full width of the support
        // is 2 * ff_support.
        double ff_support = gaussian_support(ffsd, params.max_error * 0.5);
        // Additional support for the flat portion of the DC band lowpass
        double dc_support = dc ? ff : 0;
        // The support as the number of FFT frequency bands needed,
        // allowing for both sides of the Gaussian.
        int fq_2support = int(ceil((ff_support + dc_support) * 2 * fftsize));

        band_params<T> *bp = new band_params<T>;
        bp->dc = dc;
        bp->sftsize = next_power_of_two(fq_2support);
        bp->sftsize_log2 = whichp2(bp->sftsize);
        sftsize_max = std::max(sftsize_max, bp->sftsize);
        bp->sft = pool<fft<C *>, int>::shared.get(bp->sftsize);

        bp->kernel.resize(bp->sftsize);
        bp->dual_kernel.resize(bp->sftsize);
        bp->shift_kernel.resize(bp->sftsize);
        bp->shift_kernel_conj.resize(bp->sftsize);

        if (dc)
            bp->center = 0;
        else
            bp->center = ff * fftsize;
        bp->icenter = (int)rint(bp->center);
        bp->ff = ff;
        bp->fq_offset_int = bp->icenter - (bp->sftsize >> 1);
        bp->ffsd = ffsd;

        // Calculate frequency-domain window kernel

        if (dc) {
            // The cutoff frequency is a fraction of the
            // fftsize, but gaussian_windowed_lowpass()
            // designs the filter in terms of the the
            // sftsize, so we need to scale the frequencies
            // accordingly.
            double scale = fftsize / bp->sftsize;
            gaussian_windowed_lowpass(
                ff * scale,
                ffsd * scale,
                bp->kernel.begin(), bp->kernel.end());
            fftshift(bp->kernel.begin(), bp->kernel.end());
        } else {
            for (unsigned int i = 0; i < bp->sftsize; i++) {
                // ii = large-FFT band number
                unsigned int ii = i + bp->fq_offset_int;
                // this_ff = fractional frequency of this kernel sample
                double this_ff = ii * inv_fftsize_double;
                T y = norm_gaussian(ffsd, this_ff - ff);
                bp->kernel[i] = y;
            }
        }

        return bp;
    }

    // Index of first slice affected by sample at t0
    slice_index_t affected_slice_b(sample_index_t t0) const {
        return affected_slice_0(t0, fftsize);
    }

    // Index of first slice not affected by sample at t1
    slice_index_t affected_slice_e(sample_index_t t1) const {
        return affected_slice_1(t1, fftsize);
    }

    // Sample index of the first sample in the filet of slice si
    sample_index_t slice_filet_begin(slice_index_t si) const {
        // XXX optimize using shift
        return si * (sample_index_t)filet_part(fftsize) + fat_part(fftsize);
    }

    // Analyze a single slice of signal.
    // The signal in "real_signal", which is fftsize samples long.
    // The sample time of the first sample is "t0".
    // Returns a set of spectrogram coefficients through "c1",
    // and the decimated-by-2 signal through *downsampled_dst.
    // Uses temporary buffers passed through buf0...buf3.

    void
    analyze_one_slice(int oct, T *real_signal, sample_index_t t0,
                      oct_coefs<C> &c1,
                      T *downsampled_dst,
                      pod_vector<C> &buf0, // fftsize
                      pod_vector<C> &buf1, // fftsize
                      pod_vector<C> &buf2, // largest sftsize
                      pod_vector<C> &buf3  // largest sftsize
                      ) const
    {
        zone<T> &z = *octaves[oct].z;
        assert(c1.bands.size() == z.n_bands);

        pod_vector<C> &spectrum(buf1);
#if GABORATOR_USE_REAL_FFT
        rft->transform(real_signal, spectrum.data());
#else
        // Real to complex
        pod_vector<C> &signal(buf0);
        std::copy(real_signal, real_signal + fftsize, signal.begin());
        ft->transform(signal.data(), spectrum.data());
#endif
        pod_vector<C> &tmp(buf2);

        T scale_factor = inv_fftsize_t;

        for (unsigned int obno = 0; obno < z.bandparams.size(); obno++) {
            band_params<T> *bp = z.bandparams[obno].get();
            C *sdata = tmp.data();

            // Multiply a slice of the spectrum by the frequency-
            // domain window and store in sdata.  The band center
            // frequency is at the center of the spectrum slice and
            // the center of the window, but in sdata, it needs to be
            // at the beginning to be ready for the inverse FFT.
            // Therefore, we need to perform an ifftshift.  Also,
            // we need to take care not to overrun the beginning or
            // end of the spectrum - for the dc band, we always
            // need to wrap around to negative frequencies, and
            // potentially it could happen with other bands, too
            // if they are really wide.  To avoid the overhead of
            // checking in the inner loop, use a separate slow path
            // for the rare cases where wrapping happens.

            size_t half_size = bp->sftsize >> 1;
            int start_index = bp->fq_offset_int;
            int end_index = bp->fq_offset_int + bp->sftsize;
            if (start_index >= 0 && end_index < (int)((fftsize >> 1) + 1)) {
                // Fast path: the slice lies entirely within the
                // positive-frequency half of the spectrum (including
                // DC and Nyquist).  We still need to handle the
                // positive and negative frequency halves of the slice
                // (as opposed to the whole spectrum) separately.
                // Positive frequencies
                elementwise_product(sdata, spectrum.data() + start_index + half_size,
                                    bp->kernel.data() + half_size, half_size);
                // Negative frequencies
                elementwise_product(sdata + half_size, spectrum.data() + start_index,
                                    bp->kernel.data(), half_size);
            } else {
                // Slow path
                // Positive frequencies
                for (size_t i = 0; i < half_size; i++)
                    sdata[i] = get_real_spectrum_coef(spectrum.data(),
                        start_index + half_size + i, fftsize) * bp->kernel[half_size + i];
                // Negative frequencies
                for (size_t i = 0; i < half_size; i++)
                    sdata[half_size + i] = get_real_spectrum_coef(spectrum.data(),
                        start_index + i, fftsize) * bp->kernel[i];
            }

            // Switch to time domain
            C *band = buf3.data();
            bp->sft->itransform(sdata, band);

            // Extract filet, adjust for non-integer center frequency,
            // correct phase, scale amplitude, and write to the output
            // coefficients.
            double ff = bp->center * inv_fftsize_double;
            double arg = -tau * t0 * ff;
            C phase_times_scale = C(cos(arg), sin(arg)) * scale_factor;
            elementwise_product_times_scalar(c1.bands[obno], band + fat_part(bp->sftsize),
                                             bp->shift_kernel.data() + fat_part(bp->sftsize),
                                             phase_times_scale, filet_part(bp->sftsize));
        }

        // Downsample
        if (oct + 1 < (int) n_octaves) {
            pod_vector<C> &sdata(buf2);
            // This is using a larger buffer than we actually need
            pod_vector<C> &ddata(buf0);
            assert(ddata.size() >= dsparams.sftsize);
            // Extract the low-frequency part of "spectrum" into "sdata"
            // and multiply it by the lowpass filter frequency response.
            // This means both positive and negative low frequencies.
            size_t half_size = dsparams.sftsize >> 1;
            assert(fftsize - half_size == 3 * half_size);
#if GABORATOR_USE_REAL_FFT
            // Positive frequencies
            elementwise_product(sdata.data(), spectrum.data(),
                                dsparams.kernel.data() + half_size, half_size);
            // Nyquist
            sdata[half_size] = 0;
            // Use the same buffer as the complex FFT, but as floats
            T *real_ddata = reinterpret_cast<T *>(ddata.data());
            dsparams.rsft->itransform(sdata.data(), real_ddata);
            // Beginning and end of filet part
            T *b = real_ddata + fat_part(dsparams.sftsize);
            T *e = b + filet_part(dsparams.sftsize);
            std::copy(b, e, downsampled_dst);
#else
            // Positive frequencies
            elementwise_product(sdata.data(), spectrum.data(),
                dsparams.kernel.data() + half_size, half_size);
            // Negative requencies
            elementwise_product(sdata.data() + half_size, spectrum.data() + fftsize - half_size,
                dsparams.kernel.data(), half_size);
            dsparams.sft->itransform(sdata.data(), ddata.data());
            // Beginning and end of filet part
            C *b = ddata.data() + fat_part(dsparams.sftsize);
            C *e = b + filet_part(dsparams.sftsize);
            complex2real(b, e, downsampled_dst);
#endif
        }
    }

    void
    synthesize_one_slice(int oct, const padded_coefs<C> &c,
                         const pod_vector<T> &downsampled,
                         sample_index_t t0,
                         T *signal_out,
                         pod_vector<C> &buf0, // fftsize
                         pod_vector<C> &buf2  // largest sftsize
                         ) const
    {
        zone<T> &z = *octaves[oct].z;
        pod_vector<C> &signal(buf0);
        std::fill(signal.begin(), signal.end(), 0);

        for (unsigned int obno = 0; obno < z.bandparams.size(); obno++) {
            band_params<T> *bp = z.bandparams[obno].get();

            C *indata = c.bands[obno];
            pod_vector<C> &sdata(buf2);

            T scale_factor = (T)1 / bp->sftsize;

            // Apply phase correction, adjust for non-integer center frequency,
            // and apply scale factor.  Note that phase must be calculated in double
            // precision.
            double ff = bp->center * inv_fftsize_double;
            double arg = tau * t0 * ff;
            C phase_times_scale = C(cos(arg), sin(arg)) * scale_factor;
            elementwise_product_times_scalar(sdata.data(), indata, bp->shift_kernel_conj.data(),
                                             phase_times_scale, bp->sftsize);

            // Switch to frequency domain
            bp->sft->transform(sdata.data());

            // Multiply signal spectrum by frequency-domain dual window,
            // accumulating result in signal.

            for (unsigned int i = 0; i < bp->sftsize; i++) {
                int iii = (bp->fq_offset_int + i) & (fftsize - 1);
                // Note the ifftshift of the input index, as f=0 appears in the middle
                // of the window
                C v = sdata[i ^ (bp->sftsize >> 1)] * bp->dual_kernel[i];
                // Frequency symmetry
                signal[iii] += v;
                if (! bp->dc)
                    signal[(fftsize - iii) & (fftsize - 1)] += conj(v);
            }
        }

        if (oct + 1 < (int) n_octaves) {
            // Upsample the downsampled data from the lower octaves
            pod_vector<C> &sdata(buf2);
            assert(downsampled.size() == dsparams.sftsize);
            assert(sdata.size() >= dsparams.sftsize);
#if GABORATOR_USE_REAL_FFT
            dsparams.rsft->transform(downsampled.data(), sdata.begin());
#else
            // Real to complex
            std::copy(downsampled.begin(), downsampled.end(), sdata.begin());
            dsparams.sft->transform(sdata.data());
#endif
            for (unsigned int i = 0; i < dsparams.sftsize; i++) {
                sdata[i] *= dsparams.dual_kernel[i ^ (dsparams.sftsize >> 1)];
            }

            // This implicitly zero pads the spectrum, by
            // not adding anything to the middle part
            // The splitting of the nyquist band is per
            // http://dsp.stackexchange.com/questions/14919/upsample-data-using-ffts-how-is-this-exactly-done
            // but should not really matter because
            // there should be no energy there to speak of
            // thanks to the windowing above.
            assert(dsparams.sftsize == fftsize / 2);
            unsigned int i;
            for (i = 0; i < dsparams.sftsize / 2; i++)
                signal[i] += sdata[i];
            //C nyquist = sdata[i] * (T)0.5;
            C nyquist = sdata[i] * (T)0.5;
            signal[i] += nyquist;
            signal[i + fftsize / 2] += nyquist;
            i++;
            for (;i < dsparams.sftsize; i++)
                signal[i + fftsize / 2] += sdata[i];
        }

        // Switch to time domain
#if GABORATOR_USE_REAL_FFT
        rft->itransform(signal.data(), signal_out);
#else
        ft->itransform(signal.data());
        // Copy real part to output
        complex2real(signal.begin(), signal.end(), signal_out);
#endif
    }

private:

    // Analyze a signal segment consisting of any number of samples.
    // oct is the octave; this is 0 except in recursive calls
    // real_signal points to the first sample
    // t0 is the sample time of the first sample
    // t1 is the sample time of the sample after the last sample
    // coefficients are added to msc

    void
    analyze_sliced(int oct, const T *real_signal, sample_index_t t0, sample_index_t t1, coefs<T> &msc) const {
        sliced_coefs<C> &sc = msc.octaves[oct];

        // Length of alias-free section; also the overlap period
        unsigned int filetsize = filet_part(fftsize);
        unsigned int fatsize = fat_part(fftsize);

        // Find the range of slices affected by the sample range
        slice_index_t si0 = affected_slice_b(t0);
        slice_index_t si1 = affected_slice_e(t1);

        // Length of each downsampled slice (including padding)
        unsigned int dslen = fftsize >> 1;
        // Ditto without padding
        unsigned int dsfiletlen = filet_part(dslen);
        // Total length of downsampled data
        unsigned int dstotlen = (si1 - si0) * dsfiletlen;

        // The range of sample times covered by the "downsampled" array
        sample_index_t tmp = (int64_t)si0 * (int)filetsize + (int)fatsize;
        assert((tmp & 1) == 0);
        sample_index_t dst0 = tmp >> 1;
        sample_index_t dst1 = dst0 + (int)dstotlen;

        // Not all of the "downsampled" array actually contains
        // nonzero data.  Calculate adjusted bounds to use in the
        // recursive analysis so that we don't needlessly analyze
        // zeroes.
        int ds_support = dsparams.time_support;
        sample_index_t dst0a = std::max(dst0, (t0 >> 1) - ds_support);
        sample_index_t dst1a = std::min(dst1, (t1 >> 1) + 1 + ds_support);

        pod_vector<T> downsampled(dstotlen);
        pod_vector<T> slice(fftsize);

        // Allocate buffers for analyze_one_slice(), to be shared between
        // successive calls to avoid repeated allocation
        pod_vector<C> buf0(fftsize);
        pod_vector<C> buf1(fftsize);
        pod_vector<C> buf2(sftsize_max);
        pod_vector<C> buf3(sftsize_max);

        zone<T> &z = *octaves[oct].z;

        // Temporary coefficients used in slow path below; allocated only
        // if needed, and only once per call
        ref<oct_coefs<C> > temp_coefs;

        // For each slice
        for (slice_index_t si = si0; si < si1; si++) {
            sample_index_t slice_t0 = si * (sample_index_t)filetsize;
            sample_index_t slice_t1 = slice_t0 + (sample_index_t)fftsize;

            // Intersection of signal and slice
            sample_index_t ss_t0 = std::max(t0, slice_t0);
            sample_index_t ss_t1 = std::min(t1, slice_t1);
            // If the intersection is empty, the affected_*
            // calculations above must be wrong
            assert(ss_t1 >= ss_t0);

            // If this part of the signal is all zero, we don't need
            // the coefficients, and we know the downsampling will
            // produce zero, too.
            if (is_zero(real_signal + (ss_t0 - t0),
                        real_signal + (ss_t1 - t0)))
            {
                // Gather zeroes in lieu of downsampled data
                if (oct + 1 < (int)n_octaves) {
                    unsigned int dsi0 = (si - si0) * dslen / 2;
                    unsigned int dsi1 = dsi0 + dslen / 2;
                    std::fill(downsampled.data() + dsi0, downsampled.data() + dsi1, 0);
                }
            } else {
                // For each sample in the slice
                // XXX optimize away copy when no zero fill is needed?

                copy_overlapping_zerofill(slice.data(), fftsize, real_signal, t0 - slice_t0, t1 - slice_t0);

                T *downsampled_dst = downsampled.data() + (si - si0) * filet_part(dslen);

                bool created;
                oct_coefs<C> &ssc(get_or_create_coefs_uninit(sc, si, oct, created));

                if (created) {
                    // Fast path: just store new coefficients directly
                    analyze_one_slice(oct, slice.data(), slice_t0, ssc,
                                      downsampled_dst, buf0, buf1, buf2, buf3);
                } else {
                    // Slow path: add to existing coefficients
                    if (! temp_coefs)
                        temp_coefs.reset(new oct_coefs<C>(coef_meta(oct), false));
                    analyze_one_slice(oct, slice.data(), slice_t0, *temp_coefs.get(),
                                      downsampled_dst, buf0, buf1, buf2, buf3);
                    add(z, ssc, *temp_coefs.get());
                }
            }
        }

        if (msc.shadow0) {
            // Note the "0" argument to band_scale_exp; we assume
            // band 0 in the octave has the lowest exponent and
            // therefore corresponds to the deepest level of resolution
            // pyramid needed.  If there are bands with a higher exponent
            // in the octave, their resolution pyramid will go deeper
            // than strictly necessary, but that's harmless.
            msc.shadow0->update(oct, band_scale_exp(oct, 0), sc, si0, si1);
        }

        // Recurse
        if (oct + 1 < (int)n_octaves)
            analyze_sliced(oct + 1, downsampled.data() + (dst0a - dst0), dst0a, dst1a, msc);
    }

    // Resynthesize audio from the coefficients in "msc".  The audio will
    // cover samples from t0 (inclusive) to t1 (exclusive), and is stored
    // starting at *real_signal, which must have room for (t1 - t0)
    // samples.  The octave "oct" is 0 except in recursive calls.

    void
    synthesize_sliced(int oct, const coefs<T> &msc, sample_index_t t0, sample_index_t t1, T *real_signal) const {
        const sliced_coefs<C> &sc = msc.octaves[oct];

        int filetsize = filet_part(fftsize);
        int fatsize = fat_part(fftsize);

        slice_index_t si0 = affected_slice_b(t0);
        slice_index_t si1 = affected_slice_e(t1);

        // sub_signal holds the reconstructed subsampled signal from the lower octaves,
        // for the entire time interval covered by the slices
        int sub_signal_len = ((si1 - si0) * filetsize + 2 * fatsize) / 2;
        pod_vector<T> sub_signal(sub_signal_len);
        memset(sub_signal.data(), 0, sub_signal_len * sizeof(T));
        if (oct + 1 < (int)n_octaves) {
            ssize_t sub_t0 = si0 * (filetsize / 2);
            ssize_t sub_t1 = sub_t0 + sub_signal_len;
            // Recurse
            assert(sub_t1 - sub_t0 == (ssize_t)sub_signal.size());
            synthesize_sliced(oct + 1, msc, sub_t0, sub_t1, sub_signal.data());
        }

        // Allocate buffers for synthesize_one_slice(), to be shared
        // between successive calls to avoid repeated allocation
        pod_vector<C> buf0(fftsize);
        //pod_vector<C> buf1(fftsize);
        pod_vector<C> buf2(sftsize_max);
        pod_vector<T> downsampled(dsparams.sftsize);

        // For each slice
        for (slice_index_t si = si0; si < si1; si++) {
            sample_index_t slice_t0 = si * filetsize;
            if (! sc.slices.has_index(si)) {
                // Zero fill.  Some code duplication with the copying
                // at the end of the function.
                sample_index_t b = std::max(slice_t0 + fatsize, t0);
                sample_index_t e = std::min(slice_t0 + fftsize - fatsize, t1);
                for (sample_index_t i = b; i < e; i++)
                    real_signal[i - t0] = 0;
                continue;
            }
            padded_coefs<C> c(coef_meta(oct, true), false);

            for (unsigned int obno = 0; obno < c.bands.size(); obno++) {
                unsigned int padded_len = octaves[oct].z->bandparams[obno]->sftsize;
                unsigned int len = filet_part(padded_len);
                for (int neighbor = -1; neighbor <= 1; neighbor++) {
                    slice_index_t ni = si + neighbor;
                    // Location in destination where sample 0 of source is copied,
                    // or would be if it was inside the buffer range
                    int copy_dest_offset = (len >> 1) + len * neighbor;
                    int copy_source_offset = 0;
                    int copy_len = len;
                    if (copy_dest_offset < 0) {
                        int adj = -copy_dest_offset;
                        copy_dest_offset += adj;
                        copy_source_offset += adj;
                        copy_len -= adj;
                    } else if (copy_dest_offset + copy_len >= (int)padded_len) {
                        int adj = copy_dest_offset + copy_len - padded_len;
                        copy_len -= adj;
                    }

                    C *destp = &c.bands[obno][copy_dest_offset];
                    if (ni < sc.slices.begin_index() || ni >= sc.slices.end_index()) {
                        // Zero pad
                        for (int j = 0; j < copy_len; j++)
                            destp[j] = 0;
                    } else {
                        const ref<oct_coefs<C> > &t = sc.slices.get_existing(ni);
                        if (t) {
                            const C *srcp = &t->bands[obno][copy_source_offset];
                            for (int j = 0; j < copy_len; j++)
                                destp[j] = srcp[j];
                        }
                    }
                }
            }

            // Copy downsampled signal to "downsampled" for upsampling
            if (oct + 1 < (int) n_octaves) {
                int bi = (si - si0) * filet_part(dsparams.sftsize);
                int ei = bi + dsparams.sftsize;
                assert(bi >= 0);
                assert(ei <= (int)sub_signal.size());
                std::copy(sub_signal.begin() + bi,
                          sub_signal.begin() + ei,
                          downsampled.begin());
            }

            T signal_slice[fftsize];
            synthesize_one_slice(oct, c, downsampled, slice_t0, signal_slice, buf0, buf2);

            // Copy overlapping part
            sample_index_t b = std::max(slice_t0 + fatsize, t0);
            sample_index_t e = std::min(slice_t0 + fftsize - fatsize, t1);
            for (sample_index_t i = b; i < e; i++)
                real_signal[i - t0] = signal_slice[i - slice_t0];
        }
    }

public:
    // The main analysis entry point.
    // The resulting coefficients are added to any existing coefficients in "msc".

    void analyze(const T *real_signal, sample_index_t t0, sample_index_t t1,
                 coefs<T> &msc, int n_threads = 1) const
    {
        analyze1(real_signal, t0, t1, msc, n_threads, 1);
    }

    void analyze1(const T *real_signal, sample_index_t t0, sample_index_t t1,
                  coefs<T> &msc, int n_threads, int level) const
    {
        assert(msc.octaves.size() == n_octaves);
        (void)n_threads;
        analyze_sliced(0, real_signal, t0, t1, msc);
    }

    // The main synthesis entry point

    void
    synthesize(const coefs<T> &msc, sample_index_t t0, sample_index_t t1,
               T *real_signal, int n_threads = 1) const {
        (void)n_threads;
        synthesize_sliced(0, msc, t0, t1, real_signal);
    }

    // Get an existing coefficient slice, or create a new one.
    // Note that this hides the distinction between two types
    // of nonexistence: that of slices outside the range
    // of the range_vector, and that of missing slices within
    // the range (having a null ref).  CT is the coefficient
    // type, which is typically C aka complex<T>, but can
    // be different, for example float to represent magnitudes.
    template <class CT>
    oct_coefs<CT> &get_or_create_coefs(sliced_coefs<CT> &sc,
        slice_index_t i, unsigned int oct) const
    {
        ref<oct_coefs<CT> > &p(sc.slices.get_or_create(i));
        if (! p)
            p.reset(new oct_coefs<CT>(coef_meta(oct)));
        return *p;
    }

    // As above, but return uninitialized coefficients and set created
    // to true if the coefficients did not already exist.  This is
    // just an optimization, turning a memset and add-to-memory
    // into a store.
    template <class CT>
    oct_coefs<CT> &get_or_create_coefs_uninit(sliced_coefs<CT> &sc,
        slice_index_t i, unsigned int oct, bool &created) const
    {
        ref<oct_coefs<CT> > &p(sc.slices.get_or_create(i));
        if (! p) {
            p.reset(new oct_coefs<CT>(coef_meta(oct), false));
            created = true;
        } else {
            created = false;
        }
        return *p;
    }

    // Get a pointer to an existing existing coefficient slice,
    // or null if one does not exist.  Like get_or_create_coefs(),
    // this hides the distinction between the two types of nonexistence.
    template <class CT>
    oct_coefs<CT> *get_existing_coefs(const sliced_coefs<CT> &sc,
                                      slice_index_t i) const
    {
        // XXX optimize to not lookup twice
        if (! sc.slices.has_index(i))
            return 0;
        const ref<oct_coefs<CT> > &p(sc.slices.get_existing(i));
        return p.get();
    }


    // Split a "global band number" gbno into an octave and band
    // number within octave ("obno").
    //
    // Global band numbers start at 0 for the band at or close to
    // fs/2, and increase towards lower frequencies.
    //
    // Include the DC band if "dc" is true.
    // Returns true iff valid.

    bool bno_split(int gbno, int &oct, unsigned int &obno, bool dc) const {
        if (gbno < 0) {
            // Above top octave
            return false;
        } else if (gbno < 2 * (int)params.bands_per_octave) {
            // Within top octave
            oct = 0;
            obno = 2 * params.bands_per_octave - 1 - gbno;
            return true;
        } else if (gbno < (int)n_bands_total - 1) {
            // Within a middle octave, or within non-DC part of bottom octave
            int t = gbno - 2 * params.bands_per_octave;
            assert(t >= 0);
            oct = 1 + t / params.bands_per_octave;
            obno = params.bands_per_octave - 1 - (t % params.bands_per_octave);
            if (oct == (int)n_octaves - 1) {
                // This octave may be shorter and has DC band at the
                // beginning, count from the end instead
                obno -= params.bands_per_octave;
                obno += octaves[oct].z->n_bands;
                //obno++; // Skip the DC band
            }
            assert(obno >= 0 && obno < octaves[oct].z->n_bands);
            return true;
        } else if (gbno == (int)n_bands_total - 1 && dc) {
            // DC
            oct = n_octaves - 1;
            obno = 0;
            return true;
        } else {
            return false;
        }
    }

    // The inverse of the above.  Returns a gbno.  The arguments must
    // be valid.

    int bno_merge(int oct, unsigned int obno) const {
        unsigned int n_bands = octaves[oct].z->n_bands;
        assert(obno < n_bands);
        int bno_from_end = n_bands - 1 - obno;
        return bno_from_end + octaves[oct].n_bands_above;
    }


    // Get the range of sample indices that a set of coefficients
    // pertain to.
    void get_coef_bounds(const coefs<T> &msc, sample_index_t &si0, sample_index_t &si1) const {
        // Look at the lowest-frequency band, since it has the greatest support
        // XXX what about DC?
        int gbno = n_bands_total - 2;

        int oct;
        unsigned int obno; // Band number within octave
        bool r = bno_split(gbno, oct, obno, false);
        assert(r);

        const typename sliced_coefs<C>::slices_t &slices = msc.octaves[oct].slices;
        // signal samples per band sample
        int exp = band_scale_exp(oct, obno);
        // times number of samples in band
        exp += octaves[oct].z->bandparams[obno]->sftsize_log2 - 1;
        si0 = shift_left((sample_index_t)slices.begin_index(), exp);
        si1 = shift_left((sample_index_t)slices.end_index(), exp);
    }

    // Return the time step (aka downsampling factor) of band "gbno".
    // If gbno is out of range, zero is returned.
    unsigned int band_step_log2(int gbno) const {
        int oct;
        unsigned int obno;
        bool valid = bno_split(gbno, oct, obno, true);
        if (! valid)
            return 0;
        return band_scale_exp(oct, obno);
    }

    int bandpass_bands_begin() const { return 0; }
    int bandpass_bands_end() const { return n_bands_total - 1; }

    int bands_begin() const { return 0; }
    int bands_end() const { return n_bands_total; }

    // Get the band number of the lowpass band
    int band_lowpass() const { return n_bands_total - 1; }

    // XXX simplify!
    double band_ff(int gbno) {
        if (gbno == band_lowpass())
            return 0;
        int oct;
        unsigned int obno;
        bool valid = bno_split(gbno, oct, obno, true);
        assert(valid);
        return band_ff(oct, obno);
    }


    // Convenience function to get the coefficient metadata for a given octave,
    // with or without padding
    coefs_meta &coef_meta(unsigned int oct, bool padded = false) const {
        return octaves[oct].z->cmeta[padded];
    }

    ~analyzer() {
    }

    // Get the base 2 logarithm of the downsampling factor of
    // band "obno" in octave "oct"
    int band_scale_exp(int oct, unsigned int obno) const {
        return fftsize_log2 - octaves[oct].z->bandparams[obno]->sftsize_log2 + oct;
    }
    // Members initialized in the constructor, and listed in
    // order of initialization
    parameters params;
    double band_spacing_log2;
    double band_spacing;
    double tuning_log2ff;
    double band0_log2ff;
    unsigned int n_bandpass_bands_total;
    unsigned int fftsize_log2; // log2(fftsize)
    unsigned int fftsize; // The size of the main FFT, a power of two.

    // The following members may get assigned more than once, if we
    // need to try more than one FFT size.

    double inv_fftsize_double; // 1.0 / fftsize
    T inv_fftsize_t; // 1.0f / fftsize (if using floats)
    unsigned int sftsize_max; // The size of the largest band FFT, a power of two
    unsigned int n_octaves;

    std::vector<ref<zone<T> > > zones;
    downsampling_params<T> dsparams;

    // Fourier transform object for transforming a full slice
#if GABORATOR_USE_REAL_FFT
    rfft<C *> *rft;
#else
    fft<C *> *ft;
#endif
    std::vector<octave<T> > octaves; // Per-octave parameters
    unsigned int n_bands_total; // Total number of frequency bands, including DC
    double top_band_log2ff; // log2 of fractional frequency of the highest-frequency band
    int ffref_gbno; // Band number of the reference frequency
};



// Iterate over the slices holding coefficients for a row (band)
// in the spectrogram with indices ranging from i0 to i1, and call
// the "process_existing_slice" method of the given "dest" object
// for each full or partial slice of coefficients, and/or the
// "process_missing_slice" method for each nonexistent slice.
//
// Template parameters:
//   T is the spectrogram value type
//   D is the dest object type
//   C is the coefficient type

template <class T, class D, class C = complex<T> >
struct row_foreach_slice {
    typedef C value_type;
    // With sc arg
    row_foreach_slice(const analyzer<T> &anl_,
                      const sliced_coefs<C> &sc_,
                      int oct_, unsigned int obno_):
        anl(anl_), oct(oct_), obno(obno_), sc(sc_)
    {
        init();
    }
    // With msc arg
    row_foreach_slice(const analyzer<T> &anl_,
                      const coefs<T, C> &msc,
                      int oct_, unsigned int obno_):
        anl(anl_), oct(oct_), obno(obno_), sc(msc.octaves[oct])
    {
        init();
        assert(oct < (int)msc.octaves.size());
    }
private:
    void init() {
        // This works for power-of-two-sized filets only.  To extend
        // it to other lengths, we will need a divmod function that
        // works correctly for negative arguments.
        unsigned int slice_len =
            filet_part(anl.octaves[oct].z->bandparams[obno]->sftsize);
        sh = whichp2(slice_len);
    }
public:
    void operator()(coef_index_t i0, coef_index_t i1, D &dest) const {
        assert(i0 <= i1);
        // Band size (power of two)
        int bsize = 1 << sh;
        // Adjust for t=0 being outside the filet
        int fatsize = bsize >> 1;
        i0 -= fatsize;
        i1 -= fatsize;
        coef_index_t i = i0;
        while (i < i1) {
            // Slice index
            slice_index_t sli = i >> sh;
            // Band vector index
            int bvi = i & (bsize - 1);
            int len = bsize - bvi;
            coef_index_t remain = i1 - i;
            if (remain < len)
                len = remain;
            oct_coefs<C> *c = anl.get_existing_coefs(sc, sli);
            if (c) {
                dest.process_existing_slice(c->bands[obno] + bvi, len);
            } else {
                dest.process_missing_slice(len);
            }
            i += len;
        }
    }
    const analyzer<T> &anl;
    int oct;
    unsigned int obno;
    unsigned int sh;
    const sliced_coefs<C> &sc;
};

// Helper class for row_source

template <class C, class OI>
struct writer_dest {
    writer_dest(OI output_): output(output_) { }
    void process_existing_slice(C *bv, size_t len) {
        // Can't use std::copy here because it takes the output
        // iterator by value, and using the return value does not
        // work, either.
        for (size_t i = 0; i < len; i++)
            *output++ = bv[i];
    }
    void process_missing_slice(size_t len) {
        for (size_t i = 0; i < len; i++)
            *output++ = 0;
    }
    OI output;
};

// Retrieve a sequence of coefficients from a row (band) in the
// spectrogram, with indices ranging from i0 to i1.  The indices can
// be negative, and can extend outside the available data, in which
// case zero is returned.  The coefficients are written through the
// output iterator "output".
// Template arguments:
//   T is the spectrogram value type
//   OI is the output iterator type
//   C is the coefficient value type

template <class T, class OI, class C = complex<T> >
struct row_source {
    // With sc arg
    row_source(const analyzer<T> &anl_,
               const sliced_coefs<C> &sc_,
               int oct_, unsigned int obno_):
        slicer(anl_, sc_, oct_, obno_)
    { }
    // With msc arg
    row_source(const analyzer<T> &anl_,
               const coefs<T, C> &msc_,
               int oct_, unsigned int obno_):
        slicer(anl_, msc_, oct_, obno_)
    { }
    OI operator()(coef_index_t i0, coef_index_t i1, OI output) const {
        writer_dest<C, OI> dest(output);
        slicer(i0, i1, dest);
        return dest.output;
    }
    row_foreach_slice<T, writer_dest<C, OI>, C> slicer;
};


// T -> f() -> OI::value_type

template <class F, class OI, class T>
struct transform_output_iterator: public std::iterator<std::output_iterator_tag, T> {
    typedef T value_type;
    transform_output_iterator(F f_, OI output_): f(f_), output(output_) { }
    transform_output_iterator<F, OI, T>& operator=(T v) {
        *output++ = f(v);
        return *this;
    }
    transform_output_iterator<F, OI, T>& operator*() { return *this; }
    transform_output_iterator<F, OI, T>& operator++() { return *this; }
    transform_output_iterator<F, OI, T>& operator++(int) { return *this; }
    F f;
    OI output;
};

// Apply the function f to each existing coefficient in the
// coefficient set msc.

template <class T, class F>
void apply(const analyzer<T> &anl, const coefs<T> &msc, F f) {
    typedef complex<T> C;
    unsigned int n_oct = msc.octaves.size();
    for (unsigned int oct = 0; oct < n_oct; oct++) {
        const sliced_coefs<C> &sc = msc.octaves[oct];
        slice_index_t bi = sc.slices.begin_index();
        slice_index_t ei = sc.slices.end_index();
        for (slice_index_t s = bi; s < ei; s++) {
            const ref<oct_coefs<C> > &t = sc.slices.get_existing(s);
            if (! t)
                continue;
            const oct_coefs<C> &c = *t;
            unsigned int n_bands = c.bands.size();
            for (unsigned int obno = 0; obno < n_bands; obno++) {
                C *band = c.bands[obno];
                unsigned int len = filet_part(anl.octaves[oct].z->bandparams[obno]->sftsize);
                int bno = anl.bno_merge(oct, obno);
                sample_index_t st = anl.sample_time(s, 0, oct, obno);
                int time_step = 1 << anl.band_scale_exp(oct, obno);
                for (unsigned int i = 0; i < len; i++) {
                    f(band[i], bno, st);
                    st += time_step;
                }
            }
        }
    }
}

// Apply the function f to each existing coefficient in the
// coefficient set msc within the time range st0 to st1.

template <class T, class F>
void apply(const analyzer<T> &anl, const coefs<T> &msc, F f,
           sample_index_t st0,
           sample_index_t st1)
{
    typedef complex<T> C;
    unsigned int n_oct = msc.octaves.size();
    for (unsigned int oct = 0; oct < n_oct; oct++) {
        const sliced_coefs<C> &sc = msc.octaves[oct];
        unsigned int n_bands = anl.octaves[oct].z->n_bands;
        for (unsigned int obno = 0; obno < n_bands; obno++) {
            int exp = anl.band_scale_exp(oct, obno);
            int time_step = 1 << exp;
            // Find the range of valid coefficient indices within the
            // given range of sample times.
            coef_index_t ci0 = (st0 + time_step - 1) >> exp;
            coef_index_t ci1 = (st1 + time_step - 1) >> exp;
            int bno = anl.bno_merge(oct, obno);
            // Helper class for row_foreach_slice()
            struct apply_dest {
                // st is the sample time of the first coefficient
                // sample in the range
                apply_dest(int bno_, sample_index_t st_, int step_, F f_):
                    bno(bno_), st(st_), step(step_), f(f_)
                { }
                void process_existing_slice(C *bv, size_t len) {
                    for (size_t i = 0; i < len; i++) {
                        f(bv[i], bno, st);
                        st += step;
                    }
                }
                void process_missing_slice(size_t len) {
                    st += len * step;
                }
                int bno;
                sample_index_t st;
                int step;
                F f;
            } dest(bno, shift_left(ci0, exp), 1 << exp, f);
            row_foreach_slice<T, apply_dest>(anl, sc, oct,  obno)(ci0, ci1, dest);
        }
    }
}

template <class T>
void forget_before(const analyzer<T> &anl, coefs<T> &msc,
                   sample_index_t limit)
{
    typedef complex<T> C;
    unsigned int n_oct = msc.octaves.size();
    for (unsigned int oct = 0; oct < n_oct; oct++) {
        sliced_coefs<C> &sc = msc.octaves[oct];
        // Convert limit from samples to slices, rounding down.
        // The "- 1" is because we only use the filet part.
        sample_index_t fat = fat_part(anl.fftsize) << oct;
        slice_index_t sli = (limit - fat) >> (oct + anl.fftsize_log2 - 1);
        sc.slices.erase_before(sli);
    }
}


} // namespace

#endif
