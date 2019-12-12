'''
Final Project

Author: Y. Teruyama
Date: 11/20/2019
Version: 01

'''
import math
import json

n_file = input('Input file => ')
print(in_file)
in_file = in_file.strip() #to delete any extra spaces before and after the name of input file

data = json.loads(open(in_file).read())

individuals_list = []
#Create the list of the Person class form the data file