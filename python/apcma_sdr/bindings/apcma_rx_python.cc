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
/* BINDTOOL_HEADER_FILE(apcma_rx.h)                                        */
/* BINDTOOL_HEADER_FILE_HASH(c6a33d1e3287cb9fd3ff255fcfc5d322)                     */
/***********************************************************************************/

#include <pybind11/complex.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

#include <gnuradio/apcma_sdr/apcma_rx.h>
// pydoc.h is automatically generated in the build directory
#include <apcma_rx_pydoc.h>

void bind_apcma_rx(py::module& m)
{

    using apcma_rx    = gr::apcma_sdr::apcma_rx;


    py::class_<apcma_rx, gr::block, gr::basic_block,
        std::shared_ptr<apcma_rx>>(m, "apcma_rx", D(apcma_rx))

        .def(py::init(&apcma_rx::make),
           D(apcma_rx,make)
        )
        



        ;




}








