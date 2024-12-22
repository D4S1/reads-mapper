import utils
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('mapper_output', type=str)
    parser.add_argument('reference', type=str)
    args = parser.parse_args()
    print(utils.accuracy(args.mapper_output, args.reference))