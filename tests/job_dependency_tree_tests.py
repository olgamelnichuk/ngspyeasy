#!/usr/bin/env python

import unittest
from ngspyeasy.job_dependency_tree import JobDependencyTree


class JobDependencyTreeTests(unittest.TestCase):
    def test_empty_tree(self):
        tree = JobDependencyTree()
        self.assertTrue(tree.get()[0] is None)

    def test_one_job_tree(self):
        tree = JobDependencyTree()
        tree.append("id0")
        (id, details) = tree.get()
        self.assertEquals("id0", id)
        self.assertEquals(None, details)

        (id, details) = tree.get()
        self.assertEquals(None, id)
        self.assertEquals(None, details)

    def test_duplicated_id(self):
        tree = JobDependencyTree()
        tree.append("id0")
        try:
            tree.append("id0")
            self.fail("A duplicate id exception should be raised")
        except ValueError:
            pass

    def test_non_existed_dependency(self):
        tree = JobDependencyTree()
        try:
            tree.append("id0",["non_existen_id"])
            self.fail("A non-existent dependency id exception should be raised")
        except ValueError:
            pass

    def test_job_dependency(self):
        tree = JobDependencyTree()
        tree.append("id0",)
        tree.append("id1", ["id0"])

        (id, details) = tree.get()
        self.assertEquals("id0", id)

        (id, details) = tree.get()
        self.assertEqual(None, id)

        tree.done("id0", 0)

        (id, details) = tree.get()
        self.assertEqual("id1", id)

        (id, details) = tree.get()
        self.assertEqual(None, id)

        tree.done("id0", 0)

        (id, details) = tree.get()
        self.assertEqual(None, id)
