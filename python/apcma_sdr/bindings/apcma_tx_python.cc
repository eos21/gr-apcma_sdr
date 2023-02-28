/*
 * Copyright 2023 Free Software Foundation, Inc.
 *
 * This file is part of GNU Radio
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 */

/***********************************************************************************/
/* This file is automatically generated using bindtool and can be manually edited  */
/* The following lines can be configured to regenerate this file during cmake      */
/* If manual edits are made, the following tags should be modified accordingly.    */
/* BINDTOOL_GEN_AUTOMATIC(0)                                                       */
/* BINDTOOL_USE_PYGCCXML(0)                                                        */
/* BINDTOOL_HEADER_FILE(apcma_tx.h)                                        */
/* BINDTOOL_HEADER_FILE_HASH(362218f2f47ff80e594f77572da525e5)                     */
/***********************************************************************************/

#include <pybind11/complex.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

#include <gnuradio/apcma_sdr/apcma_tx.h>
// pydoc.h is automatically generated in the build directory
#include <apcma_tx_pydoc.h>

void
    bind_apcma_tx( py::module& m ) {
    using apcma_tx = ::gr::apcma_sdr::apcma_tx;


    py::class_<apcma_tx, gr::block, gr::basic_block, std::shared_ptr<apcma_tx>>( m, "apcma_tx", D( apcma_tx ) )

        .def( py::init( &apcma_tx::make ),
              py::arg( "sf" ),
              py::arg( "samp_rate" ),
              py::arg( "os_factor" ),
              py::arg( "code_def" ),
              py::arg( "N_bits" ),
              py::arg( "subslot_width" ),
              py::arg( "start_var" ),
              py::arg( "end_var" ),
              py::arg( "number_of_sending" ),
              D( apcma_tx, make ) )


        ;
}
