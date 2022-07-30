#!/usr/bin/env python3
# -*- coding: utf-8 -*-

##########################################################################################
##########################################################################################
##
##                                Library
##
##########################################################################################
##########################################################################################

import sys, os
import pandas as pd
import logging
from logging.handlers import QueueHandler, QueueListener
import multiprocessing


##########################################################################################
##########################################################################################
##
##                                Functions
##
##########################################################################################
##########################################################################################


def logger_init(level, args):

    mpQueue = multiprocessing.Queue()

    # this is the handler for all log records
    handler = logging.FileHandler(
        filename=os.path.join(args.output, "ictv_downloading.log"), mode="w"
    )
    LOG_FORMAT_HANDLER = "%(levelname)s: %(asctime)s - %(process)s - %(message)s"
    handler.setFormatter(logging.Formatter(LOG_FORMAT_HANDLER))

    # queueListerner gets records from the queue and sends them to the handler
    queueListerner = QueueListener(mpQueue, handler)
    queueListerner.start()

    logger = logging.getLogger()
    logger.setLevel(level)
    # add the handler to the logger so records from this process are handled
    logger.addHandler(handler)

    return queueListerner, mpQueue


##########################################################################################


def init_process(mpQueue, level):

    # all records from worker processes go to queueHandler and then into mpQueue
    queueHandler = QueueHandler(mpQueue)
    logger = logging.getLogger()
    logger.setLevel(level)
    # logger.addHandler(queueHandler)


##########################################################################################


def create_folder(mypath):

    """
    Created the folder that I need to store my result if it doesn't exist
    :param mypath: path where I want the folder (write at the end of the path)
    :type: string
    :return: Nothing
    """

    try:
        os.makedirs(mypath)
    except OSError:
        pass

    return


##########################################################################################

def _make_gen(reader):
    b = reader(1024 * 1024)
    while b:
        yield b
        b = reader(1024*1024)

##########################################################################################

def rawgencount(filename):
    f = open(filename, 'rb')
    f_gen = _make_gen(f.raw.read)
    return sum( buf.count(b'\n') for buf in f_gen )

##########################################################################################
