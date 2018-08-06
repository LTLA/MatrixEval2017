#include "H5Cpp.h"
#include <vector>
#include <iostream>
#include <algorithm>

#ifndef DECACHE
#define DECACHE 0
#endif

// g++  -std=c++11 tester.cpp -lhdf5_cpp -lhdf5 -o HDF5ChunkTester

int main (int argc, const char** argv) {
    if (argc!=4) {
        std::cout << argv[0] << " [FILE] [DATASET] [BYROW]" << std::endl;
        return 1;
    }
    const char* fname=argv[1];
    const char* dname=argv[2];
    const bool userow=(argv[3][0]=='1');

    H5::H5File hfile(fname, H5F_ACC_RDONLY);
    H5::DataSet hdata=hfile.openDataSet(dname);
    H5::DataSpace hspace=hdata.getSpace();

    hsize_t dims_out[2];
    hspace.getSimpleExtentDims(dims_out, NULL);
    const size_t total_nrows=dims_out[1];
    const size_t total_ncols=dims_out[0];

    hsize_t h5_start[2], col_count[2], row_count[2];
    h5_start[0]=0;
    h5_start[1]=0;
    col_count[0]=1;
    col_count[1]=total_nrows;
    row_count[0]=total_ncols;
    row_count[1]=1;

    const H5::DataType default_type=hdata.getDataType();

    // Resetting chunk cache parameters.
    H5::DSetCreatPropList cparms = hdata.getCreatePlist();

    /* Setting up the chunk cache specification. */
    hsize_t chunk_dims[2];
    cparms.getChunk(2, chunk_dims);
    const size_t chunk_nrows=chunk_dims[1];
    const size_t chunk_ncols=chunk_dims[0];
    const size_t num_chunks_per_row=std::ceil(double(total_ncols)/chunk_ncols); // per row needs to divide by column dimensions.
    const size_t num_chunks_per_col=std::ceil(double(total_nrows)/chunk_nrows); 
#ifdef VERBOSE    
    std::cerr << "Number of row chunks is " << num_chunks_per_row << ", number of column chunks is " << num_chunks_per_col  << std::endl;
#endif

    /* Everything is transposed, so hash indices are filled column-major. 
     * Here, we computing the lowest multiple of # chunks-per-col that is greater than # chunks-per-row, plus 1.
     * This ensures that two chunks in the same row/column do not have the same hash index.
     */
#ifndef BADSLOT
    const size_t nslots = std::ceil(double(num_chunks_per_row)/num_chunks_per_col) * num_chunks_per_col + 1; 
#else
    const size_t nslots=num_chunks_per_col;
#endif    

    /* Computing the size of the cache required to store all chunks in each row or column.
     * The approach used below avoids overflow from computing eachchunk*num_Xchunks.
     */
    const size_t eachchunk=default_type.getSize() * chunk_nrows * chunk_ncols;
    const size_t eachrow=std::max<size_t>(0, eachchunk * num_chunks_per_row - DECACHE); 
    const size_t eachcol=std::max<size_t>(0, eachchunk * num_chunks_per_col - DECACHE);
#ifdef VERBOSE    
    std::cerr << "Nslots is " << nslots <<", raw data cache is " << 
        eachcol << " (column) or " << eachrow << " (row)" << std::endl;
#endif 
    
    // The first argument is ignored, according to https://support.hdfgroup.org/HDF5/doc/RM/RM_H5P.html.
    // Setting w0 to 0 to evict the last used chunk; no need to worry about full vs partial reads here.
    H5::FileAccPropList fapl(H5::FileAccPropList::DEFAULT.getId());
    fapl.setCache(0, nslots, (userow ? eachrow : eachcol), 0);

    // Reopening the file with a new HDF5 cache.
    hfile.close();
    hdata.close();
    hfile.openFile(argv[1], H5F_ACC_RDONLY, fapl);
    hdata=hfile.openDataSet(argv[2]);

#ifdef VERBOSE    
    {
        auto fapl2=hfile.getAccessPlist();
        size_t s, b;
        int m;
        double w;
        fapl2.getCache(m, s, b, w);
        std::cerr << "FileAccPropList parameters are: " <<  m << ", " << s << ", " << b << ", " << w << std::endl;
    }
#endif    

    double total=0;
    H5::DataSpace outspace;
    std::vector<double> storage(total_nrows);

    if (!userow) { 
        for (size_t c=0; c<total_ncols; ++c) {
            outspace.setExtentSimple(1, col_count+1);
            outspace.selectAll();
            h5_start[0] = c;
            hspace.selectHyperslab(H5S_SELECT_SET, col_count, h5_start);
            hdata.read(storage.data(), default_type, outspace, hspace);
            total += std::accumulate(storage.begin(), storage.end(), 0.0);
        }
    } else {
        for (size_t r=0; r<total_nrows; ++r) {
            outspace.setExtentSimple(1, row_count);
            outspace.selectAll();
            h5_start[1] = r;
            hspace.selectHyperslab(H5S_SELECT_SET, row_count, h5_start);
            hdata.read(storage.data(), default_type, outspace, hspace);
            total += std::accumulate(storage.begin(), storage.end(), 0.0);
        }
    }

    std::cout << total << std::endl;
    return 0;
}
