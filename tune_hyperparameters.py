import mapper
import argparse

def run_mapper(params: dict, reads_filename: str, genome_filename: str, out_file: str) -> None:
    """ Helper function for hyperparameter tuning. """
    mapper.main(reads_filename, 
                genome_filename, 
                params['wind_size'], 
                params['k'], 
                params['hash'], 
                params['err_max'],
                params['delta'],
                out_file
                )
    
if __name__ == '__main__':

    params = {'wind_size': 100,
              'k': 100,
              'hash': 1234567,
              'err_max': 0.1,
              'delta': 0.12
              }
    
    parser = argparse.ArgumentParser()
    parser.add_argument('reference', type=str)
    parser.add_argument('reads', type=str)
    parser.add_argument('output', type=str)
    args = parser.parse_args()

    run_mapper(params, args.reads, args.reference, args.output)
