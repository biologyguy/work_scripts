#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on: Mar 4 2015 

"""
DESCRIPTION OF PROGRAM
"""

import sys
import os
import re
import shutil
import requests
from bs4 import BeautifulSoup
import time
import datetime


def timestamp(date_string):
    if date_string == "n/a":
        return None
    new_datetime = datetime.datetime.strptime(date_string, "%m/%d/%Y").timetuple()
    return round(time.mktime(new_datetime))

if __name__ == '__main__':
    #url = "http://www.nasdaq.com/earnings/earnings-calendar.aspx"
    #content = requests.get(url).text
    with open("test_files/nasdaq.html", "r") as ifile:
        content = ifile.read()

    soup = BeautifulSoup(content, 'html5lib')
    #paginate = soup.find('div', {"class": 'paginate'}).find_all('a')

    week = soup.find(id="two_column_main_content_lreportweek").find_all('a')
    for day in week:
        print(day.getText())
    sys.exit()

    tablerows = soup.find(id='ECCompaniesTable').find_all('tr')

    comp_dict = {}
    print("#Ticker --\tMarket_cap --\tReport_date --\tLast_year's_report --\tExpected_EPS --\tLast_year's_EPS")
    for row in tablerows[1:]:
        output = ""
        columns = row.find_all("td")

        company_info = columns[1].a.getText()
        ticker = re.search("\((.*)\)", company_info).group(1)
        market_cap = re.search("\$[0-9]*\.*[0-9]*[A-Z]*|n/a", company_info).group(0)
        report_date = columns[2].getText().strip()
        last_years_report_date = columns[6].getText().strip()
        expected_eps = columns[4].getText().strip()
        last_years_eps = columns[7].getText().strip()

        if timestamp(last_years_report_date) and timestamp(report_date) - timestamp(last_years_report_date) - 31536000 + 1209600 < 0:
            output += "%s\t\t\t" % ticker
            output += "%s\t\t\t" % market_cap
            output += "%s\t\t" % report_date
            output += "%s\t\t\t\t" % last_years_report_date
            output += "%s\t\t\t" % expected_eps
            output += "%s" % last_years_eps

            print(output)
