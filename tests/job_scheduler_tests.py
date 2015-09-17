#!/usr/bin/env python

import unittest
import time
from ngspyeasy import job_scheduler


class JobSchedulerTests(unittest.TestCase):

    def setUp(self):
        job_scheduler.JobScheduler(timeout=30).start()

    def tearDown(self):
        job_scheduler.stop()
        pass

    def test_submit(self):
        job_scheduler.submit("id0","ls -la")
        job_scheduler.submit("id1","pwd",["id0"])
        job_scheduler.submit("id2","date",["id0"])
        job_scheduler.submit("id3","id",["id2"])
        job_scheduler.submit("id4","python",["id3"])
        time.sleep(5.0)