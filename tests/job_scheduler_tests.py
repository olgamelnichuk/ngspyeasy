#!/usr/bin/env python

import unittest
import time
from ngspyeasy import job_scheduler


class JobSchedulerTests(unittest.TestCase):

    def setUp(self):
        job_scheduler.JobScheduler().start()

    def tearDown(self):
        job_scheduler.stop()

    def test_submit(self):
        job_scheduler.submit("id0","ls -la")
        job_scheduler.submit("id1","pwd",["id0"])
        time.sleep(5.0)