from ._extract_errors_impl import HashesCounter

import argparse
import sourmash
from sourmash.plugins import CommandLinePlugin
import tempfile
import os
import json
from tqdm import tqdm


def main_extract_errors(args):
    
    hashes_counter = HashesCounter(
        kSize=args.ksize,
        sig_paths = args.sig_paths,
        num_threads = args.num_threads
    )
    
    hashes_counter.start_errors_extraction()

    error_kmers = hashes_counter.get_error_hashes()
    
    print("Number of error k-mers:", len(error_kmers))

    random_sig = sourmash.load_one_signature(args.sig_paths[0], ksize=args.ksize)
    final_mh = random_sig.minhash.copy_and_clear().flatten()
    
    
    if len(error_kmers) < 1_000_000:
        final_mh.add_many(error_kmers)
        finalSig = sourmash.SourmashSignature(final_mh, name="error_kmers", filename=args.out) 
        with sourmash.sourmash_args.FileOutput(args.out, 'wt') as fp:
            sourmash.save_signatures([finalSig], fp=fp)    
    else:
        finalSig = sourmash.SourmashSignature(final_mh, name="error_kmers", filename=args.out)
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


def main_filter_sigs(args):
    error_sig_path = args.errors_sig
    kSize = args.ksize
    out_dir = args.out_dir
    sig_paths = args.sig_paths
    scale = args.scale
    
    # create output directory if it doesn't exist
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # initialize HashesCounter
    hashes_counter = HashesCounter(
        kSize=kSize,
        sig_paths = sig_paths,
        num_threads = 1
    )

    print("initializing HashesCounter")
    hashes_counter.initialize_sigs_filtration(error_sig_path)


    for sig_path in tqdm(sig_paths):
        sig_basename = os.path.basename(sig_path).replace(".sig", '')
        output_file_path = os.path.join(out_dir, f"{sig_basename}.sig")
        good_hashes = hashes_counter.filter_sig_return_kmers(sig_path)
        minhash = sourmash.MinHash(n=0, ksize=kSize, scaled=scale, track_abundance=False)
        minhash.add_many(good_hashes)
        final_sig = sourmash.SourmashSignature(minhash, name=sig_basename, filename=output_file_path)
        with sourmash.sourmash_args.FileOutput(output_file_path, 'wt') as fp:
            sourmash.save_signatures([final_sig], fp=fp)
            
        
    
    
    


def cli():
    parser = argparse.ArgumentParser(description="Extract illumina kmers from a large pool of signatures")

    # Add a subparsers object to the parser
    subparsers = parser.add_subparsers(dest="command")

    # Create a parser for the "extract_kmers" subcommand
    parser_extract = subparsers.add_parser("extract", help="Extract erronous k-mers from signature")

    parser_extract.add_argument("sig_paths", nargs="+", help="signatures to process")
    parser_extract.add_argument('-k', "--ksize", type=int, required=True, help="k-mer size")
    parser_extract.add_argument('-t', "--num-threads", type=int, required=True, help="number of threads")
    parser_extract.add_argument('-o', "--out", type=str, help="output signature path", required=True)


    parser_filter = subparsers.add_parser("filter_sigs", help="Filtering out errors from signatures")
    parser_filter.add_argument("sig_paths", nargs="+", help="signatures to process")
    parser_filter.add_argument("--errors-sig", type=str, help="error signature path")
    parser_filter.add_argument('-k', "--ksize", type=int, required=True, help="k-mer size")
    parser_filter.add_argument('-s', "--scale", type=int, required=False, default=1000, help="number of threads")
    parser_filter.add_argument('-o', "--out-dir", type=str, help="output signatures directory", required=True)
    
    # Parse the arguments
    args = parser.parse_args()

    # Take actions based on the subcommand
    if args.command == "extract":
        print(f"Extracting errors")
        main_extract_errors(args)

    elif args.command == "filter_sigs":
        print(f"Filtering signatures")
        main_filter_sigs(args)
    
    
    
    
    






class Command_ExtractErrors(CommandLinePlugin):
    
    command = "extract_errors"
    description = "Extract error k-mers from a set of signatures."
    formatter_class = argparse.RawTextHelpFormatter

    
    def __init__(self, subparser):
        super().__init__(subparser)

        subparser.add_argument("sig_paths", nargs="+", help="signatures to process")
        subparser.add_argument('-k', "--ksize", type=int, required=True, help="k-mer size")
        subparser.add_argument('-t', "--num-threads", type=int, required=True, help="number of threads")
        subparser.add_argument('-o', "--out", type=str, help="output signature path", required=True)
    
    
    def main(self, args):
        
        
        hashes_counter = HashesCounter(
            kSize=args.ksize,
            sig_paths = args.sig_paths,
            num_threads = args.num_threads
        )
    
        hashes_counter.start_errors_extraction()
        print("Done processing signatures.")

        error_kmers = hashes_counter.get_error_hashes()
        print("Number of error k-mers:", len(error_kmers))

        random_sig = sourmash.load_one_signature(args.sig_paths[0], ksize=args.ksize)
        final_mh = random_sig.minhash.copy_and_clear().flatten()
        
        
        if len(error_kmers) < 1_000_000:
            final_mh.add_many(error_kmers)
            finalSig = sourmash.SourmashSignature(final_mh, name="error_kmers", filename=args.out) 
            with sourmash.sourmash_args.FileOutput(args.out, 'wt') as fp:
                sourmash.save_signatures([finalSig], fp=fp)    
        else:
            finalSig = sourmash.SourmashSignature(final_mh, name="error_kmers", filename=args.out)
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
        