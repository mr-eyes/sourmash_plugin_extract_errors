#include <nanobind/nanobind.h>
#include <iostream>
#include "simdjson.h"
#include <string>
#include <nanobind/stl/string.h>
#include <vector>
#include <parallel_hashmap/phmap.h>
#include <nanobind/stl/vector.h>
#include <stdexcept>

using phmap::flat_hash_map;
using phmap::parallel_flat_hash_map;


using MAP_KMERCOUNTER = phmap::parallel_flat_hash_map<
    uint64_t, uint32_t,
    std::hash<uint64_t>,
    std::equal_to<uint64_t>,
    std::allocator<std::pair<uint64_t, uint32_t>>,
    4, // sub-maps
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

    void get_hashes_from_sig(std::string sig_path, vector<uint64_t> &hashes)
    {
        ondemand::parser parser;
        padded_string json = padded_string::load(sig_path);
        ondemand::document sig = parser.iterate(json);

        ondemand::array sketches = sig.get_array();

        cout << "array length: " << sketches.count_elements() << endl;

        for (auto sketch : sketches)
        {
            ondemand::array signatures = sketch.find_field("signatures");

            for (auto signature : signatures)
            {
                int _ksize = (int)signature["ksize"].get_int64();
                if (this->kSize == _ksize)
                {
                    cout << "Found a signature with ksize: " << _ksize << endl;
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
                error_hashes.push_back(pair.first);
            }
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

    void process()
    {
        count_hashes();
        extract_errors();
    }

    vector<uint64_t> get_error_hashes()
    {
        return error_hashes;
    }
};

NB_MODULE(_extract_errors_impl, m)
{
    nb::class_<HashesCounter>(m, "HashesCounter")
        .def(nb::init<int, vector<string>, int>(), "kSize"_a, "sig_paths"_a, "num_threads"_a = 1)
        .def("process", &HashesCounter::process)
        .def("get_error_hashes", &HashesCounter::get_error_hashes, nb::return_stl_vector());
}
