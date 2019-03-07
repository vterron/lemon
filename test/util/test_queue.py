
from __future__ import division

import StringIO
import operator
import os
import Queue
import random
import time
import warnings
from Queue import Empty, Full


# LEMON modules
from test import unittest
from util import Queue


class QueueTest(unittest.TestCase):

    def shortDescription(self):
        """Don't use the first line of the test method's docstring."""
        pass

    def test_put_and_get(self):
        """Test that this is a FIFO queue."""

        q = Queue()
        q.put(5)
        q.put(7)
        q.put(9)

        self.assertEqual(5, q.get())
        self.assertEqual(7, q.get())
        self.assertEqual(9, q.get())

    def test_put_when_queue_is_full(self):
        """Test that if Full is raised, qsize() doesn't break."""

        q = Queue(maxsize=1)
        q.put(5)
        self.assertEqual(1, q.qsize())
        with self.assertRaises(Full):
            q.put(7, block=False)

        # Internal counter doesn't change, see:
        # https://github.com/vterron/lemon/pull/103
        self.assertEqual(1, q.qsize())

        # Pause to allow queue to complete operation before terminating
        # Otherwise sometimes results in 'IOError: [Errno 32] Broken pipe'
        time.sleep(0.1)

    def test_get_when_queue_is_empty(self):
        """Test that if Empty is raised, qsize() doesn't break."""

        q = Queue()
        self.assertEqual(0, q.qsize())
        with self.assertRaises(Empty):
            q.get(block=False)

        # Internal counter doesn't change, see:
        # https://github.com/vterron/lemon/pull/103
        self.assertEqual(0, q.qsize())
        time.sleep(0.1)

    def test_qsize(self):
        """Test that qsize() gets updated correctly as we put() and get()."""

        q = Queue()
        self.assertEqual(0, q.qsize())
        q.put(5)
        self.assertEqual(1, q.qsize())
        q.put(7)
        self.assertEqual(2, q.qsize())
        q.put(9)
        self.assertEqual(3, q.qsize())

        q.get()
        self.assertEqual(2, q.qsize())
        q.get()
        self.assertEqual(1, q.qsize())
        q.get()
        self.assertEqual(0, q.qsize())

    def test_empty(self):
        """Test empty()."""

        q = Queue()
        self.assertTrue(q.empty())
        q.put(1)
        self.assertFalse(q.empty())
        time.sleep(0.1)

    def test_clear(self):
        """Test that clear() empties the queue."""

        q = Queue()
        q.put(5)
        q.put(7)
        q.put(9)

        q.clear()
        self.assertEqual(0, q.qsize())
