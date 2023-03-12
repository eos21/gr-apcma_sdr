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
    apcma_rx::make( int sf, int samp_rate, int os_factor, int code_def, int N_bits, int subslot_width, int sliding_width, float threshold, int extremum_weight ) {
    return gnuradio::make_block_sptr<apcma_rx_impl>( sf, samp_rate, os_factor, code_def, N_bits, subslot_width, sliding_width, threshold, extremum_weight );
}


/*
 * The private constructor
 */
apcma_rx_impl::apcma_rx_impl( int sf, int samp_rate, int os_factor, int code_def, int N_bits, int subslot_width, int sliding_width, float threshold, int extremum_weight ):
    gr::block( "apcma_rx",
               gr::io_signature::make( 1, 1, sizeof( input_type ) ),
               gr::io_signature::make( 0, 0, sizeof( output_type ) ) ) {
    m_sf              = sf;
    m_samp_rate       = samp_rate;
    m_band_width      = samp_rate / os_factor;
    m_N_bits          = N_bits;
    m_os_factor       = os_factor;
    m_code_def        = code_def;
    m_sliding_width   = sliding_width;
    m_subslot_width   = subslot_width;
    m_threshold       = threshold;
    m_extremum_weight = extremum_weight;

    if ( m_sliding_width > m_subslot_width ) {
        std::cerr << "Sliding width must be smaller than subslot width!" << std::endl;
        std::exit( 1 );
    }

    subslot_cnt = 0;
    sliding_cnt = 0;

    m_number_of_bins     = ( 1u << m_sf );
    m_samples_per_symbol = m_subslot_width * m_os_factor;
    m_upchirp.resize( m_number_of_bins );
    m_downchirp.resize( m_number_of_bins );

    fftw_in  = (fftwf_complex*)fftwf_malloc( sizeof( fftwf_complex ) * m_number_of_bins );
    fftw_out = (fftwf_complex*)fftwf_malloc( sizeof( fftwf_complex ) * m_number_of_bins );
    fftw_p   = fftwf_plan_dft_1d(
        m_number_of_bins, fftw_in, fftw_out, FFTW_BACKWARD, FFTW_ESTIMATE );
    build_ref_chirps( &m_upchirp[0], &m_downchirp[0], m_sf );


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
    ninput_items_required[0] = m_os_factor * 1 << ( m_sf + 1 );
}


/**
 * @brief APCMAのコードワードを作る
 * 
 * @param do_print 標準出力にコードワードの一覧を表示させるか
 */
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
                codeword_table[i][0]                                                                     = true;
                codeword_table[i][m_number_of_bins / m_subslot_width * 2 + i]                            = true;
                codeword_table[i][m_pulse_train_length - 1 - m_number_of_bins / m_subslot_width * 2 - i] = true;
                codeword_table[i][m_pulse_train_length - 1]                                              = true;
            case 4:
                codeword_table[i][0]                                                                     = true;
                codeword_table[i][m_number_of_bins / m_subslot_width * 2 + i]                            = true;
                codeword_table[i][m_number_of_bins / m_subslot_width * 4 + i + v_mid_pulse[i] - 1]       = true;
                codeword_table[i][m_pulse_train_length - 1 - m_number_of_bins / m_subslot_width * 2 - i] = true;
                codeword_table[i][m_pulse_train_length - 1]                                              = true;
        }
        if ( do_print )
            std::cout << i << " : " << codeword_table[i] << std::endl;
    }
}


/**
 * @brief 
 * 入力信号に対しdechirp->FFT->power取得->フィルター->ピーク検出を行う
 * @param samples 入力信号
 * @param ref_chirp ダウンチャープ
 * @return boost::dynamic_bitset<>  is_peak_bin ピークが立っているbinを示すbitset
 */
boost::dynamic_bitset<>
    apcma_rx_impl::get_freq_power_peak( const gr_complex* samples,
                                        const gr_complex* ref_chirp ) {
    std::vector<gr_complex> dechirped( m_number_of_bins );
    // Multiply with ideal downchirp
    volk_32fc_x2_multiply_32fc( &dechirped[0], samples, ref_chirp, m_number_of_bins );

    // FFT using FFTW3 library
    for ( uint32_t i = 0u; i < m_number_of_bins; i++ ) {
        fftw_in[i][0] = dechirped[i].imag();
        fftw_in[i][1] = dechirped[i].real();
    }
    fftwf_execute( fftw_p );

    // Magnitudeを計算
    float fft_mag[m_number_of_bins];
    for ( uint32_t i = 0; i < m_number_of_bins; i++ ) {
        fft_mag[i] = fftw_out[i][0] * fftw_out[i][0] + fftw_out[i][1] * fftw_out[i][1];
    }

    // 長さnの平均値filter
    float    fft_mag_filtered[m_number_of_bins];
    uint32_t window_size = 1;    // 平均値フィルタの窓サイズ
    for ( uint32_t i = 0; i < m_number_of_bins; i++ ) {
        if ( i < ( window_size - 1 ) / 2 ) {    // 先頭インデックスの例外処理
            fft_mag_filtered[i] = std::accumulate( fft_mag, fft_mag + 2 * i + 1, 0.0 )
                / (float)( 2 * i + 1 );
        } else if ( m_number_of_bins - ( window_size + 1 ) / 2 < i ) {    // 末尾インデックスの例外処理
            fft_mag_filtered[i] = std::accumulate( fft_mag + m_number_of_bins - ( 2 * ( m_number_of_bins - i - 1 ) + 1 ), fft_mag + m_number_of_bins, 0.0 )
                / (float)( 2 * ( m_number_of_bins - i - 1 ) + 1 );
        } else {
            fft_mag_filtered[i] = std::accumulate( fft_mag - ( window_size - 1 ) / 2 + i, fft_mag + ( window_size + 1 ) / 2 + i, 0.0 )
                / (float)window_size;
        }
    }


    // ピークを検出
    boost::dynamic_bitset<> is_peak_bin( m_number_of_bins, false );
    for ( int i = 0; i < m_number_of_bins; i++ ) {
        bool is_extremum = true;
        for ( int j = 0; j < ( m_extremum_weight - 1 ) / 2; j++ ) {
            if ( fft_mag_filtered[( i - j < 0 ) ? 0 : ( i - j )] < fft_mag_filtered[( ( i - j - 1 ) < 0 ) ? 0 : ( i - j - 1 )] ) {
                is_extremum = false;
                break;
            }
            if ( fft_mag_filtered[( ( i + j ) > m_number_of_bins - 1 ) ? ( m_number_of_bins - 1 ) : ( i + j )] < fft_mag_filtered[( ( i + j + 1 ) > m_number_of_bins - 1 ) ? ( m_number_of_bins - 1 ) : ( i + j + 1 )] ) {
                is_extremum = false;
                break;
            }
        }

        bool is_above_threshold = ( fft_mag_filtered[i] > m_threshold );    // threholdを超えているか判断
        if ( is_extremum && is_above_threshold ) {
            is_peak_bin[i] = true;
        }
    }
    return is_peak_bin;
}


/**
 * @brief 
 * 
 * @return true subslot is on
 * @return false subslot is off
 */

/**
 * @brief 
 * 該当のsubslotにパルスの立ち上がりが存在するかを返す関数
 * @param input 入力信号のポインタ
 * @param virtual_freq_offset 周波数オフセットの参照渡し
 * @return true subslot is on 
 * @return false subslot is off
 */
bool
    apcma_rx_impl::pulse_detection( const gr_complex* input, std::vector<int>& virtual_freq_offset ) {
    bool has_detected_pulse = false;
    for ( uint32_t sliding_cnt = 0; sliding_cnt < m_subslot_width; sliding_cnt += m_sliding_width ) {
        // ピークが表れる周波数binを取得
        auto is_peak_bin = get_freq_power_peak( &input[sliding_cnt], &m_downchirp[0] );
        // subslotに対応した周波数binのいずれかがtrueか否か
        for ( uint32_t i = 0; i < m_sliding_width; i++ ) {
            uint32_t index = ( m_number_of_bins - ( m_sliding_width - 1 ) + i ) % m_number_of_bins;
            if ( is_peak_bin[index] ) {
                has_detected_pulse = true;
                virtual_freq_offset.push_back( sliding_cnt + ( m_number_of_bins - index ) % m_number_of_bins );
            }
        }
    }
    return has_detected_pulse;
}

int
    apcma_rx_impl::general_work( int                        noutput_items,
                                 gr_vector_int&             ninput_items,
                                 gr_vector_const_void_star& input_items,
                                 gr_vector_void_star&       output_items ) {
    auto in = static_cast<const input_type*>( input_items[0] );


    gr_complex       in_downed[m_number_of_bins * 2];
    std::vector<int> decoded_val;
    // Down sampling
    for ( uint32_t i = 0; i < m_number_of_bins * 2; i++ )
        in_downed[i] = in[int( i * m_os_factor )];

    // Pulse detecting
    std::vector<int> virtual_freq_offset;
    bool             is_subslot_on = pulse_detection( &in_downed[0], virtual_freq_offset );
    if ( is_subslot_on ) {
        pulse_train.set( m_pulse_train_length - 1 );
    }

    //shift register
    if ( pulse_train[0] && pulse_train[m_pulse_train_length - 1] ) {
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
    pulse_train >>= 1;

    // print decoded value
    for ( uint32_t i = 0; i < decoded_val.size(); i++ ) {
        printf( "%d\n", decoded_val[i] );
    }

    subslot_cnt++;
    consume_each( m_samples_per_symbol );
    return m_samples_per_symbol;
}
}    // namespace apcma_sdr
} /* namespace gr */
