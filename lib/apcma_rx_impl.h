/* -*- c++ -*- */
/*
 * Copyright 2023 gr-apcma_sdr author.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef INCLUDED_APCMA_SDR_APCMA_RX_IMPL_H
#define INCLUDED_APCMA_SDR_APCMA_RX_IMPL_H

#include "utilities.h"
#include <gnuradio/apcma_sdr/apcma_rx.h>
#include <fftw3.h>
#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <boost/dynamic_bitset.hpp>
#include <algorithm>
#include <iostream>
#include <vector>
extern "C" {
#include "kiss_fft.h"
}

namespace gr {
namespace apcma_sdr {

class apcma_rx_impl : public apcma_rx
{
private:
    uint32_t m_bw;        ///< Bandwidth
    uint32_t m_samp_rate; ///< Sampling rate
    uint8_t m_sf;         ///< Spreading factor
    uint8_t m_os_factor;  ///< oversampling factor

    uint8_t m_N_bits;         ///< payload length
    uint32_t m_band_width;    ///< Band width
    uint8_t m_code_def;       ///< Code definition
    uint32_t m_subslot_width; ///< Subslot width

    uint32_t m_sliding_width;      ///< Sliding window width
    float m_threshold;             ///< Power threshold for subslot detecting
    uint32_t m_number_of_bins;     ///< Number of bins in each lora Symbol
    uint32_t m_samples_per_symbol; ///< Number of samples received per lora symbols

    std::vector<gr_complex> in_downed;   ///< downsampled input
    std::vector<gr_complex> m_downchirp; ///< Reference downchirp
    std::vector<gr_complex> m_upchirp;   ///< Reference upchirp

    int64_t slot_cnt;                    ///< Number of slots already received
    uint32_t m_length_c;                 ///< Code length C of gengeral APCMA
    uint8_t m_N_pulses;                  ///< Number of pulses per single frame message
    uint32_t pulse_train_length;         ///< Pulse train length for high-throughput APCMA
    boost::dynamic_bitset<> pulse_train; ///< Pulse train which shift
    std::vector<boost::dynamic_bitset<>> codeword_tabel; ///< Codeword table

    kiss_fft_cpx* cx_in;  ///< input of the FFT
    kiss_fft_cpx* cx_out; ///< output of the FFT

    fftwf_complex* fftw_in;
    fftwf_complex* fftw_out;
    fftwf_plan fftw_p;

    uint64_t m_num_consumed_samples;


    std::vector<int> get_css_pulse_detect_kiss(const gr_complex* samples,
                                               gr_complex* ref_chirp);
    std::vector<int> get_css_pulse_detect_fftw(const gr_complex* samples,
                                               gr_complex* ref_chirp);

public:
    apcma_rx_impl(int sf,
                  int samp_rate,
                  int os_factor,
                  int code_def,
                  int N_bits,
                  int subslot_width,
                  int sliding_width,
                  float threshold);
    ~apcma_rx_impl();

    // Where all the action really happens
    void forecast(int noutput_items, gr_vector_int& ninput_items_required);

    int general_work(int noutput_items,
                     gr_vector_int& ninput_items,
                     gr_vector_const_void_star& input_items,
                     gr_vector_void_star& output_items);
};

} // namespace apcma_sdr
} // namespace gr

#endif /* INCLUDED_APCMA_SDR_APCMA_RX_IMPL_H */
