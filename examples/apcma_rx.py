#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#
# SPDX-License-Identifier: GPL-3.0
#
# GNU Radio Python Flow Graph
# Title: APCMA receiving program
# Author: Atsushi Nakamura
# GNU Radio version: 3.10.3.0

from gnuradio import apcma_sdr
from gnuradio import gr
from gnuradio.filter import firdes
from gnuradio.fft import window
import sys
import signal
from argparse import ArgumentParser
from gnuradio.eng_arg import eng_float, intx
from gnuradio import eng_notation
from gnuradio import uhd
import time


class apcma_rx(gr.top_block):

    def __init__(self):
        gr.top_block.__init__(self, "APCMA receiving program", catch_exceptions=True)

        ##################################################
        # Variables
        ##################################################
        self.threshold = threshold = 0.1
        self.subslot_width = subslot_width = 2**10
        self.sliding_width = sliding_width = 128
        self.sf = sf = 10
        self.samp_rate = samp_rate = 250000
        self.os_factor = os_factor = 1
        self.number_of_bits = number_of_bits = 12
        self.code_definition = code_definition = 4

        ##################################################
        # Blocks
        ##################################################
        self.uhd_usrp_source_0 = uhd.usrp_source(
            ",".join(("", '')),
            uhd.stream_args(
                cpu_format="fc32",
                args='',
                channels=list(range(0,1)),
            ),
        )
        self.uhd_usrp_source_0.set_samp_rate(samp_rate)
        # No synchronization enforced.

        self.uhd_usrp_source_0.set_center_freq(922.8e6, 0)
        self.uhd_usrp_source_0.set_antenna("RX2", 0)
        self.uhd_usrp_source_0.set_bandwidth((samp_rate / os_factor), 0)
        self.uhd_usrp_source_0.set_gain(65, 0)
        self.uhd_usrp_source_0.set_min_output_buffer((2**(sf+2)) * os_factor)
        self.apcma_sdr_apcma_rx_0 = apcma_sdr.apcma_rx(sf, samp_rate, os_factor, code_definition, number_of_bits, subslot_width, sliding_width, threshold)


        ##################################################
        # Connections
        ##################################################
        self.connect((self.uhd_usrp_source_0, 0), (self.apcma_sdr_apcma_rx_0, 0))


    def get_threshold(self):
        return self.threshold

    def set_threshold(self, threshold):
        self.threshold = threshold

    def get_subslot_width(self):
        return self.subslot_width

    def set_subslot_width(self, subslot_width):
        self.subslot_width = subslot_width

    def get_sliding_width(self):
        return self.sliding_width

    def set_sliding_width(self, sliding_width):
        self.sliding_width = sliding_width

    def get_sf(self):
        return self.sf

    def set_sf(self, sf):
        self.sf = sf

    def get_samp_rate(self):
        return self.samp_rate

    def set_samp_rate(self, samp_rate):
        self.samp_rate = samp_rate
        self.uhd_usrp_source_0.set_samp_rate(self.samp_rate)
        self.uhd_usrp_source_0.set_bandwidth((self.samp_rate / self.os_factor), 0)

    def get_os_factor(self):
        return self.os_factor

    def set_os_factor(self, os_factor):
        self.os_factor = os_factor
        self.uhd_usrp_source_0.set_bandwidth((self.samp_rate / self.os_factor), 0)

    def get_number_of_bits(self):
        return self.number_of_bits

    def set_number_of_bits(self, number_of_bits):
        self.number_of_bits = number_of_bits

    def get_code_definition(self):
        return self.code_definition

    def set_code_definition(self, code_definition):
        self.code_definition = code_definition




def main(top_block_cls=apcma_rx, options=None):
    tb = top_block_cls()

    def sig_handler(sig=None, frame=None):
        tb.stop()
        tb.wait()

        sys.exit(0)

    signal.signal(signal.SIGINT, sig_handler)
    signal.signal(signal.SIGTERM, sig_handler)

    tb.start()

    try:
        input('')
    except EOFError:
        pass
    tb.stop()
    tb.wait()


if __name__ == '__main__':
    main()
