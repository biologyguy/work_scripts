#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on: Apr 22 2015 

"""
Generate a similarity metric between two homolog groups files

- Split files up into list of groups, and each group is itself a list of ids
- Iterate over the subject list, comparing each group against the groups in the query list
    - Count the number of matches between groups
    - If one or more sequences in subject group are found in query group, divide the number of matching ids by the total
    size of both the query and subject
    - Sum the scores for all matching query groups
    - Stop iterating through query once all subject ids have been identified
    - Final tally for each subject group is standardized against the size of the total size of all groups in subject
    - Final score is the sum of all tallies from subject/query search, between max value of 1 and min of 0
- The metric is not symmetric between subject and query, so for a true comparison run compare() in both directions and
take the average (not currently implemented).
"""

import MyFuncs


class Clusters():
    def __init__(self, path):
        with open(path, "r") as ifile:
            self.input = ifile.read()

        self.clusters = self.input.strip().split("group")[1:]
        self.clusters = [[y for y in x.strip().split(" ")[1:]] for x in self.clusters]
        self.size = 0.
        for group in self.clusters:
            self.size += len(group)
        self.printer = MyFuncs.DynamicPrint()

    def compare(self, query_clusters):
        score = 0.
        counter = 1
        for subj in self.clusters:
            printer.write("Cluster %s of %s" % (counter, len(self.clusters)))
            counter += 1
            tally = 0.
            len_subj = len(subj)

            # This is quite inefficient, bc it iterates from the top of query list every time... Need a way to better
            # manage search.
            for query in query_clusters.clusters:
                matches = self.num_matches(subj, query)
                if not matches:
                    continue
                else:
                    tally += (matches * 2.) / (len(subj) + len(query))
                    len_subj -= matches
                    if len_subj == 0:
                        score += tally * (len(subj) / self.size)
                        break
        print("")
        return score

    @staticmethod
    def num_matches(_subj, _query):
        count = 0.
        for next_item in _subj:
            if next_item in _query:
                count += 1
        return count if count > 0 else None


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(prog="compare_homolog_groups", description="", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("subject", help="Input file 1", action="store")
    parser.add_argument("query", help="Input file 2", action="store")

    in_args = parser.parse_args()

    timer = MyFuncs.Timer()
    printer = MyFuncs.DynamicPrint()

    groups1 = Clusters(in_args.subject)
    groups2 = Clusters(in_args.query)

    print("Score: %s\n%s" % (groups1.compare(groups2), timer.end()))