from ._extract_errors_impl import HashesCounter

import argparse
import sourmash
from sourmash.plugins import CommandLinePlugin
import tempfile
import os
import json

def cli():
    
    argparser = argparse.ArgumentParser()
    argparser.add_argument("-k", "--ksize", type=int, required=True, help="k-mer size")
    argparser.add_argument("-t", "--num-threads", type=int, required=True, help="number of threads")
    argparser.add_argument("sig_paths", nargs="+", help="signatures to process", required=True)
    args = argparser.parse_args()
    
    hashes_counter = HashesCounter(
        kSize=args.ksize,
        sig_paths = args.sig_paths,
        num_threads = args.num_threads
    )
    
    hashes_counter.process()

    error_kmers = hashes_counter.get_error_hashes()

    print(len(error_kmers))






class Command_ExtractErrors(CommandLinePlugin):
    
    command = "extract_errors"
    description = "Extract error k-mers from a set of signatures."
    formatter_class = argparse.RawTextHelpFormatter

    
    def __init__(self, subparser):
        super().__init__(subparser)

        subparser.add_argument("sig_paths", nargs="+", help="signatures to process", required=True)
        subparser.add_argument('-k', "--ksize", type=int, required=True, help="k-mer size")
        subparser.add_argument('-t', "--num-threads", type=int, required=True, help="number of threads")
        subparser.add_argument('-o', "--out", type=str, required=True, help="output signature path")
    
    
    def main(self, args):
        
        
        hashes_counter = HashesCounter(
            kSize=args.ksize,
            sig_paths = args.sig_paths,
            num_threads = args.num_threads
        )
    
        hashes_counter.process()
        print("Done processing signatures.")

        error_kmers = hashes_counter.get_error_hashes()
        print("Number of error k-mers:", len(error_kmers))

        random_sig = sourmash.load_one_signature(args.sig_paths[0], ksize=args.ksize)
        final_mh = random_sig.minhash.copy_and_clear().flatten()
        finalSig = sourmash.SourmashSignature(final_mh, name="error_kmers", filename=args.out) 
        
        if len(error_kmers) < 1_000_000:
            with sourmash.sourmash_args.FileOutput(args.out, 'wt') as fp:
                sourmash.save_signatures([finalSig], fp=fp)    
        else:
            with tempfile.NamedTemporaryFile(delete=False) as temp_file:
                with sourmash.sourmash_args.FileOutput(temp_file.name, 'wt') as fp:
                    sourmash.save_signatures([finalSig], fp=fp)
                with open(temp_file.name, 'r') as file:
                    final_output_str = file.read()

            os.remove(temp_file.name)
            json_obj = json.loads(final_output_str)
            json_obj[0]['signatures'][0]['mins'] = error_kmers
            with open(args.out, 'w') as f:
                json.dump(json_obj, f, separators=(',', ':'))

        print("Final signature saved to", args.out)
        