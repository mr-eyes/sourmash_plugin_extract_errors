#include <nanobind/nanobind.h>
#include <iostream>
#include "simdjson.h"
#include <string>
#include <nanobind/stl/string.h>

using namespace std;
using namespace simdjson;
using namespace nb::literals;

namespace nb = nanobind;

void get_hashes_from_sig(std::string sig_path, int kSize)
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
            if (kSize == _ksize)
            {
                cout << "Found a signature with ksize: " << _ksize << endl;
                ondemand::array mins = signature["mins"].get_array();
                for (uint64_t hash_val : mins)
                {
                    cout << hash_val << endl;
                }
            }
        }
    }
}

NB_MODULE(_extract_errors_impl, m)
{
    m.def("get_hashes_from_sig", &get_hashes_from_sig);
}

