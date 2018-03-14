import pandas
import xlrd


def main():

    data = pandas.read_excel('Data/NIHMS413369-supplement-1_si_001.xls',
                             sheet_name='Supplementary Table 5',
                             usecols=[4, 6, 8, 9])
