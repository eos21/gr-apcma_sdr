/* -*- c++ -*- */
/*
 * Copyright 2023 gr-apcma_sdr author.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef INCLUDED_APCMA_SDR_APCMA_TX_IMPL_H
#define INCLUDED_APCMA_SDR_APCMA_TX_IMPL_H

#include "mid_pulse_list.h"
#include "utilities.h"

#include <boost/dynamic_bitset.hpp>
#include <fstream>
#include <gnuradio/apcma_sdr/apcma_tx.h>
#include <gnuradio/block.h>
#include <gnuradio/io_signature.h>
#include <inttypes.h>
#include <iostream>
#include <random>
#include <sstream>

namespace gr {
namespace apcma_sdr {
template <typename _Ty>
std::ostream &
    operator<<( std::ostream &ostr, const std::vector<_Ty> &v ) {
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
class apcma_tx_impl: public apcma_tx {
  private:
    uint8_t  m_sf;            ///< Transmission spreading factor
    uint32_t m_samp_rate;     ///< Transmission sampling rate
    uint32_t m_band_width;    ///< Transmission bandwidth
    int      m_os_factor;     ///< oversampling factor based on sampling rate and bandwidth

    uint32_t m_number_of_bins;        ///< number of bin per one chirp signal
    uint32_t m_samples_per_symbol;    ///< samples per slot

    uint8_t               m_code_def;
    uint32_t              m_send_index = 0;
    uint32_t              m_subslot_width;
    uint32_t              m_N_bits;
    uint32_t              m_N_pulses;
    uint32_t              m_length_c;
    std::vector<uint32_t> m_l_pulse_interval;
    std::vector<uint32_t> m_l_mid_pulse;
    uint32_t              message_length;

    uint32_t              m_start_var;
    uint32_t              m_end_var;
    std::vector<uint32_t> m_list_send_data;
    uint32_t              m_send_data;
    uint64_t              m_nunber_of_transmission;

    uint32_t           m_interval_slots;
    std::mt19937       mt;     // メルセンヌ・ツイスタの32ビット版
    std::random_device rnd;    // 非決定的な乱数生成器
    uint32_t           m_max_interval_slots;

    std::vector<gr_complex> m_upchirp;      ///< reference upchirp
    std::vector<gr_complex> m_downchirp;    ///< reference downchirp

    uint32_t                             m_pulse_train_length;    ///< Pulse train length for high-throughput APCMA
    boost::dynamic_bitset<>              pulse_train;             ///< Pulse train which shift
    std::vector<boost::dynamic_bitset<>> codeword_table;          ///< Codeword table

    uint64_t         m_frame_cnt;
    uint64_t         m_output_buffer_size;
    std::vector<int> init_transmission();
    void             make_codeword_table( bool do_print );

  public:
    apcma_tx_impl( uint8_t sf, uint32_t samp_rate, uint32_t os_factor, uint8_t code_def, uint32_t N_bits, uint32_t subslot_width, uint32_t start_var, uint32_t end_var, uint32_t number_of_sending );
    ~apcma_tx_impl();

    // Where all the action really happens
    void forecast( int noutput_items, gr_vector_int &ninput_items_required );

    void set_min_output_buffer( int port, long min_output_buffer );

    int general_work( int                        noutput_items,
                      gr_vector_int             &ninput_items,
                      gr_vector_const_void_star &input_items,
                      gr_vector_void_star       &output_items );
};

}    // namespace apcma_sdr
}    // namespace gr

#endif /* INCLUDED_APCMA_SDR_APCMA_TX_IMPL_H */
