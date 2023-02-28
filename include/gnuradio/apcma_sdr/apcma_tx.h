/* -*- c++ -*- */
/*
 * Copyright 2023 gr-apcma_sdr author.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef INCLUDED_APCMA_SDR_APCMA_TX_H
#define INCLUDED_APCMA_SDR_APCMA_TX_H

#include <gnuradio/apcma_sdr/api.h>
#include <gnuradio/block.h>

namespace gr {
namespace apcma_sdr {

/*!
     * \brief <+description of block+>
     * \ingroup apcma_sdr
     *
     */
class APCMA_SDR_API apcma_tx: virtual public gr::block {
  public:
    typedef std::shared_ptr<apcma_tx> sptr;

    /*!
       * \brief Return a shared_ptr to a new instance of apcma_sdr::apcma_tx.
       *
       * To avoid accidental use of raw pointers, apcma_sdr::apcma_tx's
       * constructor is in a private implementation
       * class. apcma_sdr::apcma_tx::make is the public interface for
       * creating new instances.
       */
    static sptr make( uint8_t sf, uint32_t samp_rate, uint32_t os_factor, uint8_t code_def, uint32_t N_bits, uint32_t subslot_width, uint32_t start_var, uint32_t end_var, uint32_t number_of_sending );
};

}    // namespace apcma_sdr
}    // namespace gr

#endif /* INCLUDED_APCMA_SDR_APCMA_TX_H */
