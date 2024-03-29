#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on: Mar 4 2015 

import sys
import re
import requests
from bs4 import BeautifulSoup
import time
import datetime
import MyFuncs

MIN_MARKET_CAP = 250000000.  # $250M
MIN_DAYS_EARLY = 15
MIN_DAYS_LATE = 17


def timestamp(date_string):
    if date_string == "n/a":
        return None
    new_datetime = datetime.datetime.strptime(date_string, "%m/%d/%Y").timetuple()
    return round(time.mktime(new_datetime))


def convert_value(value):  # Market Cap: input as string --> eg $12.43M
    magnitude = value[-1]
    value = float(value[1:-1])
    if magnitude == "M":
        value *= 1000000
    elif magnitude == "B":
        value *= 1000000000
    elif magnitude == "T":
        value *= 1000000000000
    else:
        sys.exit(magnitude)
    return value


def company_info(tablerows):
    for row in tablerows[1:]:
        output = ""
        columns = row.find_all("td")

        _company_info = columns[1].a.getText()
        ticker = re.search("\((.*)\)", _company_info).group(1)
        market_cap = re.search("\$[0-9]*\.*[0-9]*[A-Z]*|n/a", _company_info).group(0)
        if market_cap == "n/a" or convert_value(market_cap) < MIN_MARKET_CAP:
            continue
        if columns[0].a:
            release_time = columns[0].a.get('href')
            release_time = release_time.split("/")[-1]
        else:
            release_time = "unknown"

        report_date = columns[2].getText().strip()
        prev_report_date = columns[6].getText().strip()
        expected_eps = columns[4].getText().strip()
        last_years_eps = columns[7].getText().strip()

        if timestamp(prev_report_date) \
                and (timestamp(report_date) - timestamp(prev_report_date) - 31536000 + (80640 * MIN_DAYS_EARLY) < 0
                     or timestamp(report_date) - timestamp(prev_report_date) - 31536000 - (80640 * MIN_DAYS_LATE) > 0):
            output += "%s\t" % ticker
            output += "%s\t" % market_cap
            output += "%s\t" % report_date
            output += "%s\t" % prev_report_date
            change = ((timestamp(report_date) - timestamp(prev_report_date) - 31536000) / -86400)
            output += "%s\t" % int(change)
            output += "%s\t" % expected_eps
            output += "%s\t" % last_years_eps
            output += "%s" % release_time

            dyn_print.write("%s\n" % output)


if __name__ == '__main__':
    dyn_print = MyFuncs.DynamicPrint()
    date = datetime.date.today()
    today = date.strftime("%d")
    month = date.strftime("%m")
    valid_days = []
    for _ in range(32):
        if date.strftime("%a") not in ["Sat", "Sun"]:
            valid_days.append(date.strftime("%Y-%b-%d"))

        if date.strftime("%d") == today and date.strftime("%m") != month:
            break

        date += datetime.timedelta(days=1)

    dyn_print.write("#Ticker\tMarket cap\tReport date\tLast year's report"
                    "\tChange (days)\tExpected EPS\tLast year's EPS\tTime\n")

    for next_day in valid_days:
        dyn_print.write(next_day)
        url = "http://www.nasdaq.com/earnings/earnings-calendar.aspx?date=%s" % next_day
        content = requests.get(url).text
        soup = BeautifulSoup(content, 'html5lib')
        if soup.find(id='ECCompaniesTable'):
            company_info(soup.find(id='ECCompaniesTable').find_all('tr'))

    dyn_print.write("Done: %s - %s" % (valid_days[0], valid_days[-1]))