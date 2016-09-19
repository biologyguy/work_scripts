#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on: Sep 19 2016
from .. import practice


def test_protein_parse(capsys):
    practice.protein_parse('/Users/jonesalm/practice_data/tps1.fa')
    out, error = capsys.readouterr()
    assert "Number of Sequences: 15" in out
