
// Compile command
// g++ -I /home/haselab/miniconda3/envs/gr310/include/ -lfftw3f ./fft_time_measure.cc ../lib/kiss_fft.c
// Description
// FFT ライブラリの処理時間計測のためのプログラム
extern "C" {
#include </home/haselab/gr-apcma_sdr/lib/kiss_fft.h>
}
#include </home/haselab/miniconda3/envs/gr310/include/fftw3.h>
#include <chrono>
#include <complex>
#include <iostream>
#include <vector>
int
    main( int argc, char const *argv[] ) {
    int                              N = 128;
    std::vector<std::complex<float>> upchirp( N );
    std::vector<std::complex<float>> downchirp( N );
    fftwf_complex                   *fftw_in;
    fftwf_complex                   *fftw_out;
    fftwf_plan                       fftw_p;


    fftw_in  = (fftwf_complex *)fftwf_malloc( sizeof( fftwf_complex ) * N );
    fftw_out = (fftwf_complex *)fftwf_malloc( sizeof( fftwf_complex ) * N );

    auto start = std::chrono::system_clock::now();    // 計測スタート時刻を保存
    fftw_p     = fftwf_plan_dft_1d( N, fftw_in, fftw_out, FFTW_FORWARD, FFTW_MEASURE );
    auto end   = std::chrono::system_clock::now();    // 計測終了時刻を保存

    kiss_fft_cpx *cx_in;     ///< input of the FFT
    kiss_fft_cpx *cx_out;    ///< output of the FFT
    cx_in            = new kiss_fft_cpx[N];
    cx_out           = new kiss_fft_cpx[N];
    kiss_fft_cfg cfg = kiss_fft_alloc( N, 0, 0, 0 );

    float phase;
    for ( int n = 0; n < N; n++ ) {
        //the scaling factor of 0.9 is here to avoid to saturate the USRP_SINK
        phase        = 2.0 * M_PI * ( n * n / ( 2 * N ) - 0.5 * n );
        upchirp[n]   = std::complex<float>( float( std::cos( phase ) ), float( std::sin( phase ) ) );
        downchirp[n] = std::complex<float>( float( std::cos( -phase ) ), float( std::sin( -phase ) ) );
    }

    std::vector<std::complex<float>> dechirped( N );


    const int times = 1000 * 1000;
    // auto      start = std::chrono::system_clock::now();    // 計測スタート時刻を保存
    for ( int i = 0; i < times; i++ ) {
        for ( int ii = 0; ii < N; ii++ ) {
            dechirped[ii]  = upchirp[ii] * downchirp[ii];
            fftw_in[ii][0] = dechirped[ii].real();
            fftw_in[ii][1] = dechirped[ii].imag();
        }
        fftwf_execute( fftw_p );
    }

    // for ( int i = 0; i < times; i++ ) {
    //     for ( int ii = 0; ii < N; ii++ ) {
    //         dechirped[ii] = upchirp[ii] * downchirp[ii];
    //         cx_in[ii].r   = dechirped[ii].real();
    //         cx_in[ii].i   = dechirped[ii].imag();
    //     }
    //     kiss_fft( cfg, cx_in, cx_out );
    // }

    // auto end  = std::chrono::system_clock::now();    // 計測終了時刻を保存
    auto dur  = end - start;    // 要した時間を計算
    auto msec = std::chrono::duration_cast<std::chrono::milliseconds>( dur ).count();
    // 要した時間をミリ秒（1/1000秒）に変換して表示
    std::cout << msec << " milli sec \n";


    return 0;
}