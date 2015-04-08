__author__ = 'karl'
import re

class test:
    patterns = ["if","then"]
    strings = ["if it is cold then I will wear a coat", "if I am hungry then I will eat a snack",
               "if it is dark then I will bring a flashlight"]
    ifs = []
    thens = []

    for sentence in strings:
        match_if = re.search(patterns[0],sentence)
        match_then = re.search(patterns[1],sentence)
        start_if = match_if.end()+1
        end_if = match_then.start()-1
        ifs.append(sentence[start_if:end_if])
        start_then = match_then.end()+1
        end_then = len(sentence)
        thens.append(sentence[start_then:end_then])
    print("Ifs\t\t\t|Thens")
    for x in range(0,len(ifs)):
        print(ifs[x]+"\t|"+thens[x])