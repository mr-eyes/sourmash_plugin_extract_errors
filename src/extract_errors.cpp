#include <nanobind/nanobind.h>
#include <iostream>
#include "simdjson.h"
#include <string>
#include <nanobind/stl/string.h>
#include <vector>
#include <parallel_hashmap/phmap.h>
#include <nanobind/stl/vector.h>
#include <stdexcept>
#include <omp.h>
#include <nlohmann/json.hpp>
#include <fstream>

using json = nlohmann::json;

using phmap::flat_hash_map;
using phmap::parallel_flat_hash_map;

using MAP_KMERCOUNTER = phmap::parallel_flat_hash_map<
    uint64_t, uint32_t,
    std::hash<uint64_t>,
    std::equal_to<uint64_t>,
    std::allocator<std::pair<uint64_t, uint32_t>>,
    1, // sub-maps
    std::mutex>;

namespace nb = nanobind;
using namespace std;
using namespace simdjson;
using namespace nb::literals;

bool valid_file(std::string path)
{
    FILE *fp = fopen(path.c_str(), "r");
    if (fp == NULL)
    {
        return false;
    }
    fclose(fp);
    return true;
}

class HashesCounter
{

private:
    int kSize;
    vector<string> sig_paths;
    MAP_KMERCOUNTER kmerToCount;
    vector<uint64_t> error_hashes;
    int num_threads;
    flat_hash_map<uint64_t, bool> errors_map;

    void count_hashes()
    {
        for (string sig_path : sig_paths)
        {
            vector<uint64_t> hashes;
            get_hashes_from_sig(sig_path, hashes);
            for (uint64_t hash_val : hashes)
            {

                this->kmerToCount.try_emplace_l(
                    hash_val,
                    [](MAP_KMERCOUNTER::value_type &v)
                    { v.second += 1; },
                    1);
            }
        }
    }

    void count_hashes_parallel()
    {

        int thread_num, num_threads, start, end, vec_i;
        int n = this->sig_paths.size();
        omp_set_num_threads(this->num_threads);

#pragma omp parallel private(vec_i, thread_num, num_threads, start, end)
        {
            thread_num = omp_get_thread_num();
            num_threads = omp_get_num_threads();
            start = thread_num * n / num_threads;
            end = (thread_num + 1) * n / num_threads;
            for (vec_i = start; vec_i < end; vec_i++)
            {
                string sig_path = this->sig_paths[vec_i];
                vector<uint64_t> hashes;
                get_hashes_from_sig(sig_path, hashes);
                for (uint64_t hash_val : hashes)
                {

                    this->kmerToCount.try_emplace_l(
                        hash_val,
                        [](MAP_KMERCOUNTER::value_type &v)
                        { v.second += 1; },
                        1);
                }
            }
        }
    }

    void get_hashes_from_large_sig(std::string sigpath, vector<uint64_t> &hashes)
    {

        std::ifstream sig_stream(sigpath);
        json data = json::parse(sig_stream);

        for (auto const &sketch : data)
        {
            for (auto const &signature : sketch["signatures"])
            {
                int _ksize = (int)signature["ksize"];
                if (this->kSize == _ksize)
                {
                    for (auto const &hash_val : signature["mins"])
                    {
                        hashes.emplace_back(hash_val);
                    }
                }
            }
        }
    }

    void get_hashes_from_sig(std::string sig_path, vector<uint64_t> &hashes)
    {
        ondemand::parser parser;
        padded_string json = padded_string::load(sig_path);
        ondemand::document sig = parser.iterate(json);

        ondemand::array sketches = sig.get_array();

        for (auto sketch : sketches)
        {
            ondemand::array signatures = sketch.find_field("signatures");

            for (auto signature : signatures)
            {
                int _ksize = (int)signature["ksize"].get_int64();
                if (this->kSize == _ksize)
                {
                    ondemand::array mins = signature["mins"].get_array();
                    for (uint64_t hash_val : mins)
                    {
                        hashes.emplace_back(hash_val);
                    }
                }
            }
        }
    }

    void extract_errors()
    {
        for (auto const &pair : kmerToCount)
        {
            if (pair.second == 1)
            {
                error_hashes.emplace_back(pair.first);
            }
        }
    }

    void error_kmers_to_hashmap()
    {
        this->errors_map.reserve(this->error_hashes.size());
        for (uint64_t hash_val : this->error_hashes)
        {
            this->errors_map[hash_val] = true;
        }
        cout << "Error hashes size: " << this->errors_map.size() << endl;
    }

    void load_errors_sig(string sig_path)
    {
        if (valid_file(sig_path))
        {
            // First load the errors sig from file
            cerr << "Loading error hashes from sig file" << endl;
            get_hashes_from_large_sig(sig_path, this->error_hashes);
            // Populate it to a hashmap and remove the this->error_hashes vector
            cerr << "Loading error hashes to hashmap" << endl;
            this->error_kmers_to_hashmap();
        }
    }

public:
    HashesCounter(int kSize, vector<string> sig_paths, int num_threads = 1)
    {
        this->kSize = kSize;
        for (string sig_path : sig_paths)
        {
            if (valid_file(sig_path))
            {
                this->sig_paths.push_back(sig_path);
            }
            else
            {
                throw invalid_argument("Invalid file path: " + sig_path);
            }
        }

        this->num_threads = num_threads;
    }

    void start_errors_extraction()
    {
        if (this->num_threads > 1)
        {
            cerr << "Using " << this->num_threads << " threads" << endl;
            count_hashes_parallel();
        }
        else
        {
            count_hashes();
        }
        extract_errors();
    }

    vector<uint64_t> get_error_hashes()
    {
        return error_hashes;
    }

    void initialize_sigs_filtration(string sig_path)
    {

        if (!valid_file(sig_path))
        {
            cerr << "initializing sigs filtration" << endl;
            if (this->error_hashes.size())
            {
                this->error_kmers_to_hashmap();
            }
            else
            {
                throw invalid_argument("Error hashes not loaded, and invalid file path: " + sig_path);
            }
        }
        else
        {
            cerr << "loading errors sig" << endl;
            this->load_errors_sig(sig_path);
        }
    }

    vector<uint64_t> filter_sig_return_kmers(string sig_file_path)
    {
        // load kmers from sig file
        vector<uint64_t> hashes;
        get_hashes_from_sig(sig_file_path, hashes);
        // remove error hashes
        vector<uint64_t> filtered_hashes;
        for (uint64_t hash_val : hashes)
        {
            if (this->errors_map[hash_val])
            {
                filtered_hashes.emplace_back(hash_val);
            }
        }
        // clean some space
        hashes.clear();
        return filtered_hashes;
    }

    void dump_kmers_to_file(string file_path)
    {
        ofstream out_file(file_path);
        cerr << "dumping kmers with number: " << this->kmerToCount.size() << endl;
        for (auto const &pair : kmerToCount)
        {
            out_file << pair.first << '\t' << pair.second << endl;
        }
        out_file.close();
    }
};

NB_MODULE(_extract_errors_impl, m)
{
    nb::class_<HashesCounter>(m, "HashesCounter")
        .def(nb::init<int, vector<string>, int>(), "kSize"_a, "sig_paths"_a, "num_threads"_a = 1)
        .def("start_errors_extraction", &HashesCounter::start_errors_extraction)
        .def("get_error_hashes", &HashesCounter::get_error_hashes)
        .def("initialize_sigs_filtration", &HashesCounter::initialize_sigs_filtration, "sig_path"_a)
        .def("filter_sig_return_kmers", &HashesCounter::filter_sig_return_kmers, "sig_file_path"_a)
        .def("dump_kmers_to_file", &HashesCounter::dump_kmers_to_file, "file_path"_a);
}
