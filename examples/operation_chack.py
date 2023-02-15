#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#
# SPDX-License-Identifier: GPL-3.0
#
# GNU Radio Python Flow Graph
# Title: Operation check
# Author: Atsushi Nakamura
# GNU Radio version: 3.10.3.0

from gnuradio import analog
from gnuradio import apcma_sdr
from gnuradio import blocks
from gnuradio import gr
from gnuradio.filter import firdes
from gnuradio.fft import window
import sys
import signal
from argparse import ArgumentParser
from gnuradio.eng_arg import eng_float, intx
from gnuradio import eng_notation




class operation_chack(gr.top_block):

    def __init__(self):
        gr.top_block.__init__(self, "Operation check", catch_exceptions=True)

        ##################################################
        # Variables
        ##################################################
        self.threshold = threshold = 0.1
        self.subslot_width = subslot_width = 128
        self.sliding_width = sliding_width = 32
        self.sf = sf = 7
        self.samp_rate = samp_rate = 250000
        self.os_factor = os_factor = 1
        self.number_of_bits = number_of_bits = 8
        self.code_definition = code_definition = 4

        ##################################################
        # Blocks
        ##################################################
        self.blocks_throttle_0 = blocks.throttle(gr.sizeof_gr_complex*1, samp_rate,True)
        self.apcma_sdr_apcma_rx_0 = apcma_sdr.apcma_rx(sf, samp_rate, os_factor, code_definition, number_of_bits, subslot_width, sliding_width, threshold)
        self.analog_noise_source_x_0 = analog.noise_source_c(analog.GR_GAUSSIAN, 1, 0)


        ##################################################
        # Connections
        ##################################################
        self.connect((self.analog_noise_source_x_0, 0), (self.blocks_throttle_0, 0))
        self.connect((self.blocks_throttle_0, 0), (self.apcma_sdr_apcma_rx_0, 0))


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
        self.blocks_throttle_0.set_sample_rate(self.samp_rate)

    def get_os_factor(self):
        return self.os_factor

    def set_os_factor(self, os_factor):
        self.os_factor = os_factor

    def get_number_of_bits(self):
        return self.number_of_bits

    def set_number_of_bits(self, number_of_bits):
        self.number_of_bits = number_of_bits

    def get_code_definition(self):
        return self.code_definition

    def set_code_definition(self, code_definition):
        self.code_definition = code_definition




def main(top_block_cls=operation_chack, options=None):
    if gr.enable_realtime_scheduling() != gr.RT_OK:
        print("Error: failed to enable real-time scheduling.")
    tb = top_block_cls()

    def sig_handler(sig=None, frame=None):
        tb.stop()
        tb.wait()

        sys.exit(0)

    signal.signal(signal.SIGINT, sig_handler)
    signal.signal(signal.SIGTERM, sig_handler)

    tb.start()

    try:
        input('Press Enter to quit: ')
    except EOFError:
        pass
    tb.stop()
    tb.wait()


if __name__ == '__main__':
    main()
