from ._extract_errors_impl import HashesCounter

import argparse


def cli():
    
    argparser = argparse.ArgumentParser()
    argparser.add_argument("-k", "--ksize", type=int, required=True, help="k-mer size")
    argparser.add_argument("-t", "--num-threads", type=int, required=True, help="number of threads")
    argparser.add_argument("sig_paths", nargs="+")
    args = argparser.parse_args()
    
    hashes_counter = HashesCounter(
        kSize=args.ksize,
        sig_paths = args.sig_paths,
        num_threads = args.num_threads
    )
    
    hashes_counter.process()

    error_kmers = hashes_counter.get_error_hashes()

    print(len(error_kmers))



import sourmash
from sourmash.plugins import CommandLinePlugin


class Command_ExtractErrors(CommandLinePlugin):
    
    command = "extract_errors"
    description = "Extract error k-mers from a set of signatures."
    formatter_class = argparse.RawTextHelpFormatter

    
    def __init__(self, subparser):
        super().__init__(subparser)

        subparser.add_argument("sig_paths", nargs="+")
        subparser.add_argument("--ksize", type=int, default=31)
        subparser.add_argument("--num-threads", type=int, default=1)
        subparser.add_argument("--out", type=str, default=None)
    
    
    def main(self, args):
        
        hashes_counter = HashesCounter(
            kSize=args.ksize,
            sig_paths = args.sig_paths,
            num_threads = args.num_threads
        )
    
        
        hashes_counter.process()

        error_kmers = hashes_counter.get_error_hashes()
        print("Number of error k-mers:", len(error_kmers))

        random_sig = sourmash.load_one_signature(args.sig_paths[0], ksize=args.ksize)
        final_mh = random_sig.minhash.copy_and_clear().flatten()
        final_mh.add_many(error_kmers)
        finalSig = sourmash.SourmashSignature(final_mh, name="error_kmers", filename=args.out)
        with sourmash.sourmash_args.FileOutput(args.out, 'wt') as fp:
            sourmash.save_signatures([finalSig], fp=fp)

        print("Final signature saved to", args.out)
        