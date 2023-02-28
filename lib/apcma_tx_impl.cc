/* -*- c++ -*- */
/*
 * Copyright 2023 gr-apcma_sdr author.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "apcma_tx_impl.h"


namespace gr {
namespace apcma_sdr {

using input_type  = float;
using output_type = gr_complex;
apcma_tx::sptr
    apcma_tx::make( uint8_t sf, uint32_t samp_rate, uint32_t os_factor, uint8_t code_def, uint32_t N_bits, uint32_t subslot_width, uint32_t start_var, uint32_t end_var, uint32_t number_of_sending ) {
    return gnuradio::make_block_sptr<apcma_tx_impl>(
        sf, samp_rate, os_factor, code_def, N_bits, subslot_width, start_var, end_var, number_of_sending );
}


/*
     * The private constructor
     */
apcma_tx_impl::apcma_tx_impl( uint8_t sf, uint32_t samp_rate, uint32_t os_factor, uint8_t code_def, uint32_t N_bits, uint32_t subslot_width, uint32_t start_var, uint32_t end_var, uint32_t number_of_sending ):
    gr::block( "apcma_tx",
               gr::io_signature::make( 0, 0, sizeof( input_type ) ),
               gr::io_signature::make( 1, 1, sizeof( output_type ) ) ) {
    m_sf         = sf;
    m_samp_rate  = samp_rate;
    m_os_factor  = os_factor;
    m_band_width = samp_rate / os_factor;

    m_code_def      = code_def;
    m_N_bits        = N_bits;
    m_subslot_width = subslot_width;

    m_number_of_bins     = ( 1u << m_sf );
    m_samples_per_symbol = m_number_of_bins * m_os_factor;
    m_upchirp.resize( m_number_of_bins );
    m_downchirp.resize( m_number_of_bins );

    m_start_var = start_var;
    m_end_var   = end_var;
    m_list_send_data.resize( m_end_var - m_start_var + 1 );
    for ( int i = 0; i < ( m_end_var - m_start_var + 1 ); i++ ) {
        m_list_send_data[i] = m_start_var + i;
    }
    m_nunber_of_transmission = number_of_sending;


    m_l_pulse_interval.resize( m_N_pulses - 1 );
    m_l_mid_pulse = get_mid_pulse_vec( m_code_def, m_N_bits );

    m_upchirp.resize( m_samples_per_symbol );
    m_downchirp.resize( m_samples_per_symbol );
    build_ref_chirps( &m_upchirp[0], &m_downchirp[0], m_sf );

    m_output_buffer_size = 512;
    m_frame_cnt          = 0;
}

/*
     * Our virtual destructor.
     */
apcma_tx_impl::~apcma_tx_impl() {
}

void
    apcma_tx_impl::forecast( int noutput_items, gr_vector_int &ninput_items_required ) {
    ninput_items_required[0] = 0;
}

void
    apcma_tx_impl::set_min_output_buffer( int port, long min_output_buffer ) {
    port              = 0;
    min_output_buffer = m_output_buffer_size;
}

void
    apcma_tx_impl::make_codeword_table( bool do_print ) {
    // APCMA パラメータの設定
    switch ( m_code_def ) {
        case 0:    // Andi-4pulse
            m_length_c = ( 1u << ( m_N_bits ) ) * 2 + 5;
            m_N_pulses = 4;    // パルスの数を定義
            break;
        case 1:    // PaDi-4pulse(1)
            m_length_c = ( 1u << ( m_N_bits ) ) * 2 + 5;
            m_N_pulses = 4;    // パルスの数を定義
            break;
        case 2:    // PaDi-4pulse(2)
            m_length_c = ( 1u << ( m_N_bits ) ) * 2 + 5;
            m_N_pulses = 4;    // パルスの数を定義
            break;
        case 3:    // PaDi-4pulse(3)
            m_length_c = ( 1u << ( m_N_bits ) ) * 2 + 6;
            m_N_pulses = 4;    // パルスの数を定義
            break;
        case 4:    // Andi-5pulse
            m_length_c = ( 1u << ( m_N_bits ) ) * 3 + 5;
            m_N_pulses = 5;    // パルスの数を定義
            break;
        case 5:    // PaDi-5pulse
            m_length_c = ( 1u << ( m_N_bits ) ) * 3 + 6;
            m_N_pulses = 5;    // パルスの数を定義
        default:
            std::cerr << "Error occured. Code definition doesn't match any of them."
                      << std::endl;
            break;
    }

    /*
    この式の意味 :   第1項 -> 最終パルス以外のguard slotを含めた部分
                    第2項 -> guard slot以外のパルス間隔
                    第3項 -> 最終パルス
    */
    m_pulse_train_length = m_number_of_bins / m_subslot_width * ( 2 * ( m_N_pulses - 1 ) ) + ( m_length_c - ( m_N_pulses * 2 - 1 ) ) + 1;
    pulse_train.resize( m_pulse_train_length, false );

    auto v_mid_pulse = get_mid_pulse_vec( m_code_def, m_N_bits );
    codeword_table.resize( 1u << m_N_bits );
    for ( int i = 0; i < ( 1u << m_N_bits ); i++ ) {
        codeword_table[i].resize( m_pulse_train_length, false );
    }
    for ( int i = 0; i < ( 1u << m_N_bits ); i++ ) {
        switch ( m_code_def ) {
            case 0:
                // 最初のパルス
                // 1stパルス&guard slot(=m_samples_per_symbol / m_subslot_width * 2) + send data(=i)
                // 最後のパルスのインデックス(=m_pulse_train_length - 1) - (1stパルス&guard slot(=m_samples_per_symbol / m_subslot_width * 2) + send data(=i))
                // 最後のパルス(=m_pulse_train_length - 1)
                codeword_table[i][0]                                              = true;
                codeword_table[i][m_samples_per_symbol / m_subslot_width * 2 + i] = true;

                codeword_table[i][m_pulse_train_length - 1 - ( m_samples_per_symbol / m_subslot_width * 2 + i )] = true;
                codeword_table[i][m_pulse_train_length - 1]                                                      = true;
            case 4:
                // 最初のパルス
                // 1stパルス&guard slot(=m_samples_per_symbol / m_subslot_width * 2) + send data(=i)
                // 最後のパルスのインデックス(=m_pulse_train_length - 1) - (1stパルス&guard slot(=m_samples_per_symbol / m_subslot_width * 2) + send data(=i))
                // 最後のパルス(=m_pulse_train_length - 1)
                codeword_table[i][0]                                                                         = true;
                codeword_table[i][m_samples_per_symbol / m_subslot_width * 2 + i]                            = true;
                codeword_table[i][m_samples_per_symbol / m_subslot_width * 4 + i + v_mid_pulse[i] - 1]       = true;
                codeword_table[i][m_pulse_train_length - 1 - m_samples_per_symbol / m_subslot_width * 2 - i] = true;
                codeword_table[i][m_pulse_train_length - 1]                                                  = true;
        }
        if ( do_print )
            std::cout << i << " : " << codeword_table[i] << std::endl;
    }
}

std::vector<int>
    apcma_tx_impl::init_transmission() {
    //////////////送信データを決めるブロック/////////////////////////////////////
    m_send_data = m_list_send_data[m_send_index];
    printf( "Frame_cnt: %d,\tSend data: %d\n", m_frame_cnt, m_send_data );
    m_send_index = ( m_send_index + 1 ) % ( m_list_send_data.size() );
    ///////////////////////////////////////////////////////////////////////////
}

int
    apcma_tx_impl::general_work( int                        noutput_items,
                                 gr_vector_int             &ninput_items,
                                 gr_vector_const_void_star &input_items,
                                 gr_vector_void_star       &output_items ) {
    auto in  = static_cast<const input_type *>( input_items[0] );
    auto out = static_cast<output_type *>( output_items[0] );
    if ( m_frame_cnt == m_nunber_of_transmission ) {
        printf( "Finished sending!\n" );
        std::exit( 1 );
    }


    m_frame_cnt++;
    consume_each( noutput_items );
    return noutput_items;
}

} /* namespace apcma_sdr */
} /* namespace gr */
