#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = 'karl'
import re
import urllib.request


class IfThan():
    def __init__(self):
        book_url = input("Please input a book url:")
        book_text = urllib.request.urlopen(book_url).read()
        book_text = book_text.decode("utf-8").split(".")
        output = self.parse_book_text(book_text)
        length = len(output[0])
        for x in range(0, length):
            print(output[0][x] + " " + output[1][x])

    def parse_book_text(self, strings):
        patterns = ["if", "then"]
        output = []
        output.append([])
        output.append([])
        for sentence in strings:
            match_if = re.search(patterns[0], sentence)
            match_then = re.search(patterns[1], sentence)
            if not (match_if == None or match_then == None):
                start_if = match_if.end() + 1
                end_if = match_then.start() - 1
                output[0].append(sentence[start_if:end_if])
                start_then = match_then.end() + 1
                end_then = len(sentence)
                output[1].append(sentence[start_then:end_then])
        return output


if __name__ == "__main__":
    test_obj = IfThan()