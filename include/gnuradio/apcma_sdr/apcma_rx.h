/* -*- c++ -*- */
/*
 * Copyright 2023 gr-apcma_sdr author.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef INCLUDED_APCMA_SDR_APCMA_RX_H
#define INCLUDED_APCMA_SDR_APCMA_RX_H

#include <gnuradio/apcma_sdr/api.h>
#include <gnuradio/block.h>

namespace gr {
namespace apcma_sdr {

/*!
 * \brief <+description of block+>
 * \ingroup apcma_sdr
 *
 */
class APCMA_SDR_API apcma_rx: virtual public gr::block {
  public:
    typedef std::shared_ptr<apcma_rx> sptr;

    /*!
     * \brief Return a shared_ptr to a new instance of apcma_sdr::apcma_rx.
     *
     * To avoid accidental use of raw pointers, apcma_sdr::apcma_rx's
     * constructor is in a private implementation
     * class. apcma_sdr::apcma_rx::make is the public interface for
     * creating new instances.
     */
    static sptr make( int   sf,
                      int   samp_rate,
                      int   os_factor,
                      int   code_def,
                      int   N_bits,
                      int   subslot_width,
                      int   sliding_width,
                      float threshold );
};

}    // namespace apcma_sdr
}    // namespace gr

#endif /* INCLUDED_APCMA_SDR_APCMA_RX_H */
