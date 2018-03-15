import pandas


def main():
    data = pandas.read_excel('Data/NIHMS413369-supplement-1_si_001.xls',
                             sheet_name='Supplementary Table 1',
                             usecols=[4, 9])
    names = list(data.get(
                'SwissProt Annotation/Homology with Highest Percentage'))
    spectra = list(data.get('Number of Identified Spectra'))
    spectra_dict = {}
    for index, string in enumerate(names):
        if str(string) == 'nan':
            continue
        ID = string.split(' ')[0]
        spectra_dict[ID] = int(spectra[index])

<<<<<<< HEAD

=======
>>>>>>> 99a8dadb20aa8e52a2932c4f51ca615c38fabc67
if __name__ == '__main__':
    main()
