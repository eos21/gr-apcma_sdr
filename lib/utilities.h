#ifndef INCLUDED_UTILITIES_H
#define INCLUDED_UTILITIES_H

#include <algorithm>
#include <cstdint>
#include <gnuradio/expj.h>
#include <iomanip>
#include <numeric>
#include <string.h>
#include <sys/resource.h>
#include <sys/syscall.h>
#include <volk/volk.h>

namespace gr {
namespace apcma_sdr {

/**
 *  \brief  Return an modulated upchirp using s_f=bw
 *
 *  \param  chirp
 *          The pointer to the modulated upchirp
 *  \param  id
 *          The number used to modulate the chirp
 * \param   sf
 *          The spreading factor to use
 * \param os_factor
 *          The oversampling factor used to generate the upchirp
 */
inline void
    build_upchirp( gr_complex* chirp, uint32_t id, uint8_t sf, uint8_t os_factor = 1 ) {
    double N      = ( 1 << sf );
    int    n_fold = N * os_factor - id * os_factor;
    for ( int n = 0; n < N * os_factor; n++ ) {
        if ( n < n_fold )
            chirp[n] =
                gr_complex( 1.0, 0.0 ) * gr_expj( 2.0 * M_PI * ( n * n / ( 2 * N ) / pow( os_factor, 2 ) + ( id / N - 0.5 ) * n / os_factor ) );
        else
            chirp[n] =
                gr_complex( 1.0, 0.0 ) * gr_expj( 2.0 * M_PI * ( n * n / ( 2 * N ) / pow( os_factor, 2 ) + ( id / N - 1.5 ) * n / os_factor ) );
    }
}

/**
 *  \brief  Return the reference chirps using s_f=bw
 *
 *  \param  upchirp
 *          The pointer to the reference upchirp
 *  \param  downchirp
 *          The pointer to the reference downchirp
 * \param   sf
 *          The spreading factor to use
 */
inline void
    build_ref_chirps( gr_complex* upchirp,
                      gr_complex* downchirp,
                      uint8_t     sf,
                      uint8_t     os_factor = 1 ) {
    double N = ( 1 << sf );
    build_upchirp( upchirp, 0, sf, os_factor );
    volk_32fc_conjugate_32fc( &downchirp[0], &upchirp[0], N * os_factor );

    // for(uint n = 0; n < N ;n++){
    //     //the scaling factor of 0.9 is here to avoid to saturate the USRP_SINK
    //     upchirp[n] =  gr_complex(0.9f, 0.0f)*gr_expj(2.0 * M_PI * (n*n/(2*N)-0.5*n));
    //     downchirp[n] = gr_complex(0.9f, 0.0f)*gr_expj(-2.0 * M_PI * (n*n/(2*N)-0.5*n));
    // }
}
}    // namespace apcma_sdr
}    // namespace gr
#endif /* INCLUDED_UTILITIES_H */