'''
    Script to scour metabolite exchange graphs for differences between
    directories.

    Test from command line with:
        python ImageDiff.py -1  Unit\ Tests/Image1 -2 Unit\ Tests/Image2
'''

import os
import argparse
import numpy as np
from PIL import Image
from PIL import ImageChops


def compare_images(path_one, path_two):
    """
    Prints a list of metabolites whose behavior changes between conditions

    @param: path_one: The path to the first folder
    @param: path_two: The path to the second folder
    """
    different = [['Metabolite', 'Directories']]
    same = []
    print()
    directory = os.fsencode(path_one)
    for filename in os.listdir(directory):
        filename = str(filename).split('\'')[1]
        image_one = Image.open(path_one + '/' + filename)
        if os.path.isfile(path_two + '/' + filename):
            image_two = Image.open(path_two + '/' + filename)
        else:
            different.append([filename, 'dir1 only'])
            continue
        diff = ImageChops.difference(image_one, image_two)
        h = diff.histogram()
        blah = np.array(h, dtype=int)
        okay = (np.arange(1024) % 256)**2
        rms = np.sqrt(np.sum(blah*okay)/float(image_one.size[0] *
                                              image_two.size[1]))
        if rms > 0:  # no cutoff tolerance defined yet.
            different.append([filename, 'both'])
            if not os.path.exists('diff'):
                os.mkdir('diff')
            diff.save('diff/' + filename + '.png')
        else:
            same.append(filename)

    directory = os.fsencode(path_two)  # check for metabolites not in dir1
    directory = str(os.fsencode(path_two)).split('\'')[1]
    for filename in os.listdir(directory):
        # filename = str(filename).split('\'')[1]
        if not os.path.isfile(path_one + '/' + filename):
            different.append([filename, 'dir2 only'])
    '''

    '''
    for met in different:
        print('{:<25} {:>20}'.format(* met))
    print()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='This script checks if images have changed between'
        'conditions.')
parser.add_argument('-1', '--folder1', help='Path to folder 1', required=True)
parser.add_argument('-2', '--folder2', help='Path to folder 2', required=True)
args = parser.parse_args()
compare_images(args.folder1, args.folder2)
