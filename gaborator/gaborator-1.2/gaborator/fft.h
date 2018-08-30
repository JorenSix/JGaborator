//
// Fast Fourier transform
//
// Copyright (C) 2016-2018 Andreas Gustafsson.  This file is part of
// the Gaborator library source distribution.  See the file LICENSE at
// the top level of the distribution for license information.
//

#ifndef _GABORATOR_FFT_H
#define _GABORATOR_FFT_H

#include "gaborator/fft_naive.h"

#if GABORATOR_USE_VDSP
#include "gaborator/fft_vdsp.h"
#elif GABORATOR_USE_PFFFT
#include "gaborator/fft_pffft.h"
#endif

#endif
