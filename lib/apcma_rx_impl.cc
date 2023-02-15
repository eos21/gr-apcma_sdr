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

using input_type = gr_complex;
using output_type = float;
apcma_rx::sptr apcma_rx::make(int sf,
                              int samp_rate,
                              int os_factor,
                              int code_def,
                              int N_bits,
                              int subslot_width,
                              int sliding_width,
                              float threshold)
{
    return gnuradio::make_block_sptr<apcma_rx_impl>(sf,
                                                    samp_rate,
                                                    os_factor,
                                                    code_def,
                                                    N_bits,
                                                    subslot_width,
                                                    sliding_width,
                                                    threshold);
}


/*
 * The private constructor
 */
apcma_rx_impl::apcma_rx_impl(int sf,
                             int samp_rate,
                             int os_factor,
                             int code_def,
                             int N_bits,
                             int subslot_width,
                             int sliding_width,
                             float threshold)
    : gr::block("apcma_rx",
                gr::io_signature::make(1, 1, sizeof(input_type)),
                gr::io_signature::make(0, 0, sizeof(output_type)))
{
    m_sf = sf;
    m_samp_rate = samp_rate;
    m_band_width = samp_rate / os_factor;
    m_N_bits = N_bits;
    m_os_factor = os_factor;
    m_code_def = code_def;
    m_sliding_width = sliding_width;
    m_subslot_width = subslot_width;
    m_threshold = threshold;
    slot_cnt = 0;
    m_num_consumed_samples = 0;

    m_number_of_bins = (1u << m_sf);
    m_samples_per_symbol = m_number_of_bins * m_os_factor;
    m_upchirp.resize(m_number_of_bins);
    m_downchirp.resize(m_number_of_bins);
    in_downed.resize(m_number_of_bins);

    fftw_in = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * m_number_of_bins);
    fftw_out = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * m_number_of_bins);
    fftw_p = fftwf_plan_dft_1d(
        m_number_of_bins, fftw_in, fftw_out, FFTW_BACKWARD, FFTW_MEASURE);
    build_ref_chirps(&m_upchirp[0], &m_downchirp[0], m_sf);

    // prepare for kiss_fft
    cx_in = new kiss_fft_cpx[m_number_of_bins];
    cx_out = new kiss_fft_cpx[m_number_of_bins];

    switch (m_code_def) {
    case 0: // Andi-4pulse
        m_length_c = std::pow(2, m_N_bits + 1) + 5;
        m_N_pulses = 4; // パルスの数を定義
        break;
    case 1: // PaDi-4pulse(1)
        m_length_c = std::pow(2, m_N_bits + 1) + 5;
        m_N_pulses = 4; // パルスの数を定義
        break;
    case 2: // PaDi-4pulse(2)
        m_length_c = std::pow(2, m_N_bits + 1) + 5;
        m_N_pulses = 4; // パルスの数を定義
        break;
    case 3: // PaDi-4pulse(3)
        m_length_c = std::pow(2, m_N_bits + 1) + 6;
        m_N_pulses = 4; // パルスの数を定義
        break;
    case 4: // Andi-5pulse
        m_length_c = std::pow(2, m_N_bits) * 3 + 5;
        m_N_pulses = 5; // パルスの数を定義
        break;
    case 5: // PaDi-5pulse
        m_length_c = std::pow(2, m_N_bits) * 3 + 6;
        m_N_pulses = 5; // パルスの数を定義
    default:
        std::cerr << "Error occured. Code definition doesn't match any of them."
                  << std::endl;
        break;
    }
}

/*
 * Our virtual destructor.
 */
apcma_rx_impl::~apcma_rx_impl()
{
    fftwf_free(fftw_in);
    fftwf_free(fftw_out);
    fftwf_destroy_plan(fftw_p);
}

void apcma_rx_impl::forecast(int noutput_items, gr_vector_int& ninput_items_required)
{
    ninput_items_required[0] = m_os_factor * 1 << (m_sf + 2);
}

/// @brief 入力されたIQ sampleをDechirp +
/// FFTしてthresholdを超えた周波数binのvectorをreturn
/// @param samples
/// @param ref_chirp
/// @return bin_above_threshold
std::vector<int> apcma_rx_impl::get_css_pulse_detect_kiss(const gr_complex* samples,
                                                          gr_complex* ref_chirp)
{
    float fft_mag[m_number_of_bins];
    std::vector<gr_complex> dechirped(m_number_of_bins);

    // Prepare for kiss FFT
    kiss_fft_cfg cfg = kiss_fft_alloc(m_number_of_bins, 0, 0, 0);

    // Multiply with ideal downchirp
    volk_32fc_x2_multiply_32fc(&dechirped[0], samples, ref_chirp, m_number_of_bins);
    for (uint32_t i = 0; i < m_number_of_bins; i++) {
        cx_in[i].r = dechirped[i].real();
        cx_in[i].i = dechirped[i].imag();
    }
    // do the FFT
    kiss_fft(cfg, cx_in, cx_out);
    std::vector<int> bin_above_threshold;
    // Get freq. offset above threshold
    for (uint32_t i = 0u; i < m_number_of_bins; i++) {
        fft_mag[i] = cx_out[i].r * cx_out[i].r + cx_out[i].i * cx_out[i].i;
        if (fft_mag[i] > m_threshold) {
            bin_above_threshold.push_back(i);
            // printf("%3d,\t", i);
        }
    }
    free(cfg);

    return bin_above_threshold;
}

/// @brief 入力されたIQ sampleをDechirp +
/// FFTしてthresholdを超えた周波数binのvectorをreturn
/// @param samples
/// @param ref_chirp
/// @return bin_above_threshold
std::vector<int> apcma_rx_impl::get_css_pulse_detect_fftw(const gr_complex* samples,
                                                          gr_complex* ref_chirp)
{
    std::vector<gr_complex> dechirped(m_number_of_bins);

    // Multiply with ideal downchirp
    volk_32fc_x2_multiply_32fc(&dechirped[0], samples, ref_chirp, m_number_of_bins);
    // FFT using FFTW3 library
    for (uint32_t i = 0u; i < m_number_of_bins; i++) {
        fftw_in[i][0] = dechirped[i].imag();
        fftw_in[i][1] = dechirped[i].real();
    }

    fftwf_execute(fftw_p);
    // Get freq. offset above threshold
    std::vector<int> bin_above_threshold;
    float fft_mag[m_number_of_bins];

    for (uint32_t i = 0u; i < m_number_of_bins; i++) {
        fft_mag[i] = fftw_out[i][0] * fftw_out[i][0] + fftw_out[i][1] * fftw_out[i][1];
        if (fft_mag[i] > m_threshold) {
            bin_above_threshold.push_back(i);
        }
    }
    return bin_above_threshold;
}

int apcma_rx_impl::general_work(int noutput_items,
                                gr_vector_int& ninput_items,
                                gr_vector_const_void_star& input_items,
                                gr_vector_void_star& output_items)
{
    auto in = static_cast<const input_type*>(input_items[0]);

    // Down sampling
    for (uint32_t i = 0; i < m_number_of_bins; i++)
        in_downed[i] = in[int(i * m_os_factor)];
    // Pulse detecting
    std::vector<int> detected_result =
        get_css_pulse_detect_fftw(&in_downed[0], &m_downchirp[0]);


    // 検出した周波数binから必要なbinだけfilterする
    // 必要なbinは2^sf-(m_sliding_width-1), ... , 2^sf-1, 2^sf(=0)
    std::vector<int> filterd_detected_result;
    for (int i = 0; i < m_sliding_width; i++) {
        std::binary_search()
    }


    // if (detected_result.size() != 0) {
    //     printf("%d:::", slot_cnt);
    //     for (uint32_t i = 0; i < detected_result.size(); i++) {
    //         printf("%3d,\t", detected_result[i]);
    //     }
    //     std::cout << std::endl;
    // }


    m_num_consumed_samples += m_sliding_width;
    slot_cnt++;
    consume_each(m_sliding_width);
    return m_sliding_width;
}

} /* namespace apcma_sdr */
} /* namespace gr */
