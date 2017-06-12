#include "H5Cpp.h"
#include <vector>
#include <iostream>
#include <algorithm>

// g++  -std=c++11 tester.cpp -lhdf5_cpp -lhdf5

int main (int argc, const char** argv) {
    if (argc!=2) {
        std::cout << "URgh!" << std::endl;
        return 1;
    }

    H5::H5File hfile(argv[1], H5F_ACC_RDONLY);
    H5::DataSet hdata=hfile.openDataSet("yyy");
    H5::DataSpace hspace=hdata.getSpace();

    hsize_t dims_out[2];
    hspace.getSimpleExtentDims(dims_out, NULL);
    const size_t NR=dims_out[1];
    const size_t NC=dims_out[0];

    hsize_t h5_start[2], col_count[2], row_count[2];
    h5_start[0]=0;
    h5_start[1]=0;
    col_count[0]=1;
    col_count[1]=NR;
    row_count[0]=NC;
    row_count[1]=1;

    const H5::DataType default_type=hdata.getDataType();
    std::vector<double> storage(NR);

    // Resetting chunk cache parameters.
    H5::DSetCreatPropList cparms = hdata.getCreatePlist();
    hsize_t chunk_dims[2];
    cparms.getChunk(2, chunk_dims);
    const size_t chunk_nrows=chunk_dims[1];
    const size_t chunk_ncols=chunk_dims[0];
    const size_t num_rowchunks=std::ceil(double(NR)/chunk_nrows); // taking the ceiling.
    const size_t num_colchunks=std::ceil(double(NC)/chunk_ncols); 
    std::cout << "Number of row chunks is " << num_rowchunks << ", number of column chunks is " << num_colchunks  << std::endl;

    // Everything is transposed, so hash indices are filled column-major.
    // Computing the lowest multiple of # row-chunks that is greater than # col-chunks, plus 1.
    const size_t nslots = num_rowchunks * num_colchunks; // size_t(std::ceil(double(num_colchunks)/num_rowchunks) * num_rowchunks + 1; 

    // Computing the size of the cache required to store all chunks in each row or column.
    const size_t eachchunk=default_type.getSize() * chunk_nrows * chunk_ncols;
    const size_t eachrow=eachchunk * num_colchunks;
    const size_t eachcol=eachchunk * num_rowchunks;
    std::cout << "Nslots is " << nslots <<", raw data cache is " << eachcol << std::endl;
    
    // The first argument is ignored, according to https://support.hdfgroup.org/HDF5/doc/RM/RM_H5P.html.
    // Setting w0 to 0 to evict the last used chunk; no need to worry about full vs partial reads here.

    H5::FileAccPropList fapl(H5::FileAccPropList::DEFAULT);
    fapl.setCache(0, nslots, eachcol, 0);
    hfile.close();
    hdata.close();
    hfile.openFile(argv[1], H5F_ACC_RDONLY, fapl);
    hdata=hfile.openDataSet("yyy");
    
    double total=0;
    H5::DataSpace colspace;

    for (size_t c=0; c<100; ++c) {
        colspace.setExtentSimple(1, col_count+1);
        colspace.selectAll();
        h5_start[0] = c;
        hspace.selectHyperslab(H5S_SELECT_SET, col_count, h5_start);
        hdata.read(storage.data(), default_type, colspace, hspace);
        total += std::accumulate(storage.begin(), storage.end(), 0.0);
    }

    std::cout << total << std::endl;
    return 0;
}
