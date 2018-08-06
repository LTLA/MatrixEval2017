#include "H5Cpp.h"
#include <vector>
#include <iostream>
#include <algorithm>

#define MULT 5

/* To compile with the Rhdf5lib libraries, open R and run:

 cat(sprintf("g++ -std=c++11 -I%s access_test.cpp -o HDF5AccessTester %s -ldl\n",
     system.file(package="Rhdf5lib", "include"),
     capture.output(Rhdf5lib::pkgconfig())))

 * Then, assuming 'beach_check.R' has already been run, we can test speed with:

 system.time(system("./HDF5AccessTester bycol.h5 yyy 0"))
 system.time(system("./HDF5AccessTester bycol.h5 yyy 1"))
 */

int main (int argc, const char** argv) {
    if (argc!=4) {
        std::cout << argv[0] << " [FILE] [DATASET] [ONE_CALL]" << std::endl;
        return 1;
    }
    const char* fname=argv[1];
    const char* dname=argv[2];
    const bool one_call=(argv[3][0]=='1');

    H5::H5File hfile(fname, H5F_ACC_RDONLY);
    H5::DataSet hdata=hfile.openDataSet(dname);
    H5::DataSpace hspace=hdata.getSpace();

    hsize_t dims_out[2];
    hspace.getSimpleExtentDims(dims_out, NULL);
    const size_t total_nrows=dims_out[1];
    const size_t total_ncols=dims_out[0];

    /* Setting up the chunk cache specification, ASSUMING COLUMN CHUNKS. */
    H5::DSetCreatPropList cparms = hdata.getCreatePlist();
    hsize_t chunk_dims[2];
    cparms.getChunk(2, chunk_dims);

    const size_t chunk_nrows=chunk_dims[1];
    const size_t chunk_ncols=chunk_dims[0];
    H5::FileAccPropList fapl(H5::FileAccPropList::DEFAULT.getId());
    fapl.setCache(0, 1, hdata.getDataType().getSize() * chunk_nrows, 0);

    // Reopening the file with a new HDF5 cache.
    hfile.close();
    hdata.close();
    hfile.openFile(argv[1], H5F_ACC_RDONLY, fapl);
    hdata=hfile.openDataSet(argv[2]);

    // Assorted odds and ends.
    double total=0;
    hsize_t h5_start[2], col_count[2];
    h5_start[0]=0;
    h5_start[1]=0;
    col_count[0]=1;
    col_count[1]=total_nrows;
    const size_t N=(total_ncols/MULT);

    // Let's say we want to extract every 10th column, either via one "read()" call or multiple.
    if (one_call) { 
        hsize_t output_dims[2];
        output_dims[0]=total_nrows;
        output_dims[1]=N;
        H5::DataSpace outspace(2, output_dims);
        outspace.selectAll();

        hspace.selectNone();
        for (size_t c=0; c<N; ++c) {
            h5_start[0] = c*MULT;
            hspace.selectHyperslab(H5S_SELECT_OR, col_count, h5_start);
        }

        std::vector<double> storage(total_nrows*N);
        hdata.read(storage.data(), H5::PredType::NATIVE_DOUBLE, outspace, hspace);
        total = std::accumulate(storage.begin(), storage.end(), 0.0);
    } else {
        H5::DataSpace outspace(1, col_count+1);
        outspace.selectAll();
        std::vector<double> storage(total_nrows);

        for (size_t c=0; c<N; ++c) {
            h5_start[0] = c*MULT;
            hspace.selectHyperslab(H5S_SELECT_SET, col_count, h5_start);
            hdata.read(storage.data(), H5::PredType::NATIVE_DOUBLE, outspace, hspace);
            total += std::accumulate(storage.begin(), storage.end(), 0.0);
        }
    }

    std::cout << total << std::endl;
    return 0;
}
