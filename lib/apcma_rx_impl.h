/* -*- c++ -*- */
/*
 * Copyright 2023 gr-apcma_sdr author.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef INCLUDED_APCMA_SDR_APCMA_RX_IMPL_H
#define INCLUDED_APCMA_SDR_APCMA_RX_IMPL_H

#include "mid_pulse_list.h"
#include "utilities.h"

#include <algorithm>
#include <boost/dynamic_bitset.hpp>
#include <fftw3.h>
#include <gnuradio/apcma_sdr/apcma_rx.h>
#include <inttypes.h>
#include <iostream>
#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <vector>

namespace gr {
namespace apcma_sdr {

template <typename _Ty>
std::ostream&
    operator<<( std::ostream& ostr, const std::vector<_Ty>& v ) {
    if ( v.empty() ) {
        ostr << "{ }";
        return ostr;
    }
    ostr << "{" << v.front();
    for ( auto itr = ++v.begin(); itr != v.end(); itr++ ) {
        ostr << ", " << *itr;
    }
    ostr << "}";
    return ostr;
}

class apcma_rx_impl: public apcma_rx {
  private:
    uint32_t m_bw;           ///< Bandwidth
    uint32_t m_samp_rate;    ///< Sampling rate
    uint8_t  m_sf;           ///< Spreading factor
    uint8_t  m_os_factor;    ///< oversampling factor

    uint8_t  m_N_bits;           ///< payload length
    uint32_t m_band_width;       ///< Band width
    uint8_t  m_code_def;         ///< Code definition
    uint32_t m_subslot_width;    ///< Subslot width

    uint32_t m_sliding_width;         ///< Sliding window width
    float    m_threshold;             ///< Power threshold for subslot detecting
    uint32_t m_number_of_bins;        ///< Number of bins in each lora Symbol
    uint32_t m_samples_per_symbol;    ///< Number of samples received per slots

    std::vector<gr_complex> in_downed;      ///< downsampled input
    std::vector<gr_complex> m_downchirp;    ///< Reference downchirp
    std::vector<gr_complex> m_upchirp;      ///< Reference upchirp

    uint32_t                             m_length_c;              ///< Code length C of gengeral APCMA
    uint8_t                              m_N_pulses;              ///< Number of pulses per single frame message
    uint32_t                             m_pulse_train_length;    ///< Pulse train length for high-throughput APCMA
    boost::dynamic_bitset<>              pulse_train;             ///< Pulse train which shift
    std::vector<boost::dynamic_bitset<>> codeword_table;          ///< Codeword table

    fftwf_complex* fftw_in;
    fftwf_complex* fftw_out;
    fftwf_plan     fftw_p;

    uint64_t subslot_cnt;    ///< Number of slots already received
    uint64_t sliding_cnt;
    uint32_t previous_decoded_val;

    void make_codeword_table( bool do_print );
    void get_freq_power_peak( const gr_complex* samples,
                              const gr_complex* ref_chirp,
                              bool*             is_peak_bin );
    bool pulse_detection( std::vector<int>& virtual_freq_offset );

  public:
    apcma_rx_impl( int   sf,
                   int   samp_rate,
                   int   os_factor,
                   int   code_def,
                   int   N_bits,
                   int   subslot_width,
                   int   sliding_width,
                   float threshold );
    ~apcma_rx_impl();

    // Where all the action really happens
    void forecast( int noutput_items, gr_vector_int& ninput_items_required );

    int general_work( int                        noutput_items,
                      gr_vector_int&             ninput_items,
                      gr_vector_const_void_star& input_items,
                      gr_vector_void_star&       output_items );
};

}    // namespace apcma_sdr
}    // namespace gr

#endif /* INCLUDED_APCMA_SDR_APCMA_RX_IMPL_H */
