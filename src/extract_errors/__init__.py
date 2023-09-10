from ._extract_errors_impl import HashesCounter

import argparse




def cli():
    
    argparser = argparse.ArgumentParser()
    argparser.add_argument("--ksize", type=int, default=31)
    argparser.add_argument("--num-threads", type=int, default=1)
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