/* -*- c++ -*- */
/*
 * Copyright 2023 gr-apcma_sdr author.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "apcma_rx_impl.h"

#include <gnuradio/io_signature.h>

namespace gr {
namespace apcma_sdr {

using input_type  = gr_complex;
using output_type = float;
apcma_rx::sptr
    apcma_rx::make( int sf, int samp_rate, int os_factor, int code_def, int N_bits, int subslot_width, int sliding_width, float threshold ) {
    return gnuradio::make_block_sptr<apcma_rx_impl>( sf, samp_rate, os_factor, code_def, N_bits, subslot_width, sliding_width, threshold );
}


/*
 * The private constructor
 */
apcma_rx_impl::apcma_rx_impl( int sf, int samp_rate, int os_factor, int code_def, int N_bits, int subslot_width, int sliding_width, float threshold ):
    gr::block( "apcma_rx",
               gr::io_signature::make( 1, 1, sizeof( input_type ) ),
               gr::io_signature::make( 0, 0, sizeof( output_type ) ) ) {
    m_sf            = sf;
    m_samp_rate     = samp_rate;
    m_band_width    = samp_rate / os_factor;
    m_N_bits        = N_bits;
    m_os_factor     = os_factor;
    m_code_def      = code_def;
    m_sliding_width = sliding_width;
    m_subslot_width = subslot_width;
    m_threshold     = threshold;

    if ( m_sliding_width > m_subslot_width ) {
        std::cerr << "Sliding width must be smaller than subslot width!" << std::endl;
        std::exit( 1 );
    }

    subslot_cnt = 0;
    sliding_cnt = 0;

    m_number_of_bins     = ( 1u << m_sf );
    m_samples_per_symbol = m_number_of_bins * m_os_factor;
    m_upchirp.resize( m_number_of_bins );
    m_downchirp.resize( m_number_of_bins );
    in_downed.resize( m_number_of_bins );

    fftw_in  = (fftwf_complex*)fftwf_malloc( sizeof( fftwf_complex ) * m_number_of_bins );
    fftw_out = (fftwf_complex*)fftwf_malloc( sizeof( fftwf_complex ) * m_number_of_bins );
    fftw_p   = fftwf_plan_dft_1d(
        m_number_of_bins, fftw_in, fftw_out, FFTW_BACKWARD, FFTW_ESTIMATE );
    build_ref_chirps( &m_upchirp[0], &m_downchirp[0], m_sf );

    // prepare for kiss_fft
    cx_in  = new kiss_fft_cpx[m_number_of_bins];
    cx_out = new kiss_fft_cpx[m_number_of_bins];

    switch ( m_code_def ) {
        case 0:    // Andi-4pulse
            m_length_c = std::pow( 2, m_N_bits + 1 ) + 5;
            m_N_pulses = 4;    // パルスの数を定義
            break;
        case 1:    // PaDi-4pulse(1)
            m_length_c = std::pow( 2, m_N_bits + 1 ) + 5;
            m_N_pulses = 4;    // パルスの数を定義
            break;
        case 2:    // PaDi-4pulse(2)
            m_length_c = std::pow( 2, m_N_bits + 1 ) + 5;
            m_N_pulses = 4;    // パルスの数を定義
            break;
        case 3:    // PaDi-4pulse(3)
            m_length_c = std::pow( 2, m_N_bits + 1 ) + 6;
            m_N_pulses = 4;    // パルスの数を定義
            break;
        case 4:    // Andi-5pulse
            m_length_c = std::pow( 2, m_N_bits ) * 3 + 5;
            m_N_pulses = 5;    // パルスの数を定義
            break;
        case 5:    // PaDi-5pulse
            m_length_c = std::pow( 2, m_N_bits ) * 3 + 6;
            m_N_pulses = 5;    // パルスの数を定義
        default:
            std::cerr << "Error occured. Code definition doesn't match any of them."
                      << std::endl;
            break;
    }

    make_codeword_table( false );
}

/*
 * Our virtual destructor.
 */
apcma_rx_impl::~apcma_rx_impl() {
    fftwf_free( fftw_in );
    fftwf_free( fftw_out );
    fftwf_destroy_plan( fftw_p );
}

void
    apcma_rx_impl::forecast( int noutput_items, gr_vector_int& ninput_items_required ) {
    ninput_items_required[0] = m_os_factor * 1 << ( m_sf + 2 );
}

void
    apcma_rx_impl::make_codeword_table( bool do_print ) {
    /*
    この式の意味 :   第1項 -> 最終パルス以外のguard slotを含めた部分
                    第2項 -> guard slot以外のパルス間隔
                    第3項 -> 最終パルス
    */
    m_pulse_train_length = m_number_of_bins / m_subslot_width * ( 2 * ( m_N_pulses - 1 ) ) + ( m_length_c - ( m_N_pulses * 2 - 1 ) ) + 1;
    pulse_train.resize( m_pulse_train_length, false );

    auto v_mid_pulse = get_mid_pulse_vec( m_code_def, m_N_bits );
    codeword_table.resize( std::pow( 2, m_N_bits ) );
    for ( int i = 0; i < std::pow( 2, m_N_bits ); i++ ) {
        codeword_table[i].resize( m_pulse_train_length, false );
    }
    for ( int i = 0; i < std::pow( 2, m_N_bits ); i++ ) {
        switch ( m_code_def ) {
            case 0:
                codeword_table[i][0]                                                                         = true;
                codeword_table[i][m_samples_per_symbol / m_subslot_width * 2 + i]                            = true;
                codeword_table[i][m_pulse_train_length - 1 - m_samples_per_symbol / m_subslot_width * 2 - i] = true;
                codeword_table[i][m_pulse_train_length - 1]                                                  = true;
            case 4:
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

/// @brief 入力されたIQ sampleをDechirp +
/// FFTしてthresholdを超えた周波数binのvectorをreturn
/// @param samples
/// @param ref_chirp
/// @return bin_above_threshold
bool
    apcma_rx_impl::get_css_pulse_detect_kiss( const gr_complex* samples,
                                              gr_complex*       ref_chirp ) {
    float                   fft_mag[m_number_of_bins];
    std::vector<gr_complex> dechirped( m_number_of_bins );

    // Prepare for kiss FFT
    kiss_fft_cfg cfg = kiss_fft_alloc( m_number_of_bins, 0, 0, 0 );

    // Multiply with ideal downchirp
    volk_32fc_x2_multiply_32fc( &dechirped[0], samples, ref_chirp, m_number_of_bins );
    for ( uint32_t i = 0; i < m_number_of_bins; i++ ) {
        cx_in[i].r = dechirped[i].real();
        cx_in[i].i = dechirped[i].imag();
    }
    // do the FFT
    kiss_fft( cfg, cx_in, cx_out );
    bool             has_detected_pulse;
    std::vector<int> bin_above_threshold;

    // Get freq. offset above threshold
    for ( uint32_t i = m_number_of_bins - m_sliding_width; i < m_number_of_bins; i++ ) {
        uint32_t index = ( m_number_of_bins - ( m_sliding_width - 1 ) + i ) % m_number_of_bins;
        fft_mag[index] = cx_out[index].r * cx_out[index].r + cx_out[index].i * cx_out[index].i;
        if ( fft_mag[index] > m_threshold ) {
            has_detected_pulse = true;
        }
    }
    free( cfg );

    return has_detected_pulse;
}

/// @brief 入力されたIQ sampleをDechirp +
/// FFTしてthresholdを超えた周波数binのvectorをreturn
/// @param samples
/// @param ref_chirp
/// @return bin_above_threshold
bool
    apcma_rx_impl::get_css_pulse_detect_fftw( const gr_complex* samples,
                                              gr_complex*       ref_chirp ) {
    std::vector<gr_complex> dechirped( m_number_of_bins );

    // Multiply with ideal downchirp
    volk_32fc_x2_multiply_32fc( &dechirped[0], samples, ref_chirp, m_number_of_bins );
    // FFT using FFTW3 library
    for ( uint32_t i = 0u; i < m_number_of_bins; i++ ) {
        fftw_in[i][0] = dechirped[i].imag();
        fftw_in[i][1] = dechirped[i].real();
    }

    fftwf_execute( fftw_p );
    // Get pulse detection if power is above threshold
    std::vector<int> bin_above_threshold;
    float            fft_mag[m_number_of_bins];
    bool             has_detected_pulse = false;

    for ( uint32_t i = 0; i < m_sliding_width; i++ ) {
        uint32_t index = ( m_number_of_bins - ( m_sliding_width - 1 ) + i ) % m_number_of_bins;
        fft_mag[index] = fftw_out[index][0] * fftw_out[index][0] + fftw_out[index][1] * fftw_out[index][1];
        if ( fft_mag[index] > m_threshold ) {
            has_detected_pulse = true;
        }
    }
    return has_detected_pulse;
}


std::vector<int>
    apcma_rx_impl::shift_register( std::vector<int> fft_result ) {
    std::vector<int> decoded_val;
    // 検出した周波数binから必要なbinがある場合
    // 必要なbinは2^sf-(m_sliding_width-1), ... , 2^sf-1, 2^sf(=0)
    bool             exist_needed_bins = false;
    for ( uint32_t i = 0; i < m_sliding_width; i++ ) {
        if ( std::binary_search( fft_result.begin(),
                                 fft_result.end(),
                                 ( m_number_of_bins - ( m_sliding_width - 1 ) + i ) % m_number_of_bins ) ) {
            exist_needed_bins = true;
            break;
        }
    }
    if ( sliding_cnt == 0 ) {
        pulse_train >>= 1;
    }
    if ( exist_needed_bins ) {
        pulse_train.set( pulse_train.size() - 1 );

        if ( pulse_train[0] ) {
#pragma omp parallel
            {
#pragma omp for
                for ( int i = 0; i < ( 1 << m_N_bits ); i++ ) {
#pragma omp critical( crit_cout )
                    {
                        if ( codeword_table[i].is_subset_of( pulse_train ) ) {
                            decoded_val.push_back( i );
                        }
                    }
                }
#pragma omp barrier
            }
        }
    }

    return decoded_val;
}    // namespace apcma_sdr


int
    apcma_rx_impl::general_work( int                        noutput_items,
                                 gr_vector_int&             ninput_items,
                                 gr_vector_const_void_star& input_items,
                                 gr_vector_void_star&       output_items ) {
    auto in = static_cast<const input_type*>( input_items[0] );

    std::vector<int> decoded_val;
    // Down sampling
    for ( uint32_t i = 0; i < m_number_of_bins; i++ )
        in_downed[i] = in[int( i * m_os_factor )];
    // Pulse detecting
    bool has_detected_pulse =
        get_css_pulse_detect_fftw( &in_downed[0], &m_downchirp[0] );
    if ( has_detected_pulse ) {
        pulse_train.set( m_pulse_train_length - 1 );
    }


    sliding_cnt++;
    if ( sliding_cnt == ( m_subslot_width / m_sliding_width ) ) {
        if ( pulse_train[0] && pulse_train[m_pulse_train_length - 1] ) {
            for ( uint32_t i = 0; i < ( 1 << m_N_bits ); i++ ) {
                if ( codeword_table[i].is_subset_of( pulse_train ) ) {
                    decoded_val.push_back( i );
                }
            }
        }

        for ( int i = 0; i < decoded_val.size(); i++ ) {
            printf( "Slot: %d,\tData: %d\n", subslot_cnt, decoded_val[i] );
        }
        sliding_cnt = 0;
        pulse_train >>= 1;
        subslot_cnt++;
    }


    consume_each( m_sliding_width );
    return m_sliding_width;
}

}    // namespace apcma_sdr
} /* namespace gr */
