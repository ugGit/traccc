/*
 * TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include "traccc/edm/cell.hpp"
#include "traccc/edm/cluster.hpp"
#include "traccc/edm/measurement.hpp"
#include "vecmem/containers/vector.hpp"

#include "traccc/stdpar/clusterization/component_connection_ssv.hpp"
#include "traccc/stdpar/utils/CountingIterator.hpp"
#include <algorithm>
#include <execution>
#include <vector> // TODO: replace with vecmem vector, once code is running
#include <iostream>
#include <atomic>

// define local constant values
namespace {
using index_t = unsigned short;
/*
 * It is very likely, that small partitions will perform better,
 * because this algo does not control very fine grained the execution
 * pattern, and therefore, e.g. no graphs are created.
 */
static constexpr std::size_t MIN_CELLS_PER_PARTITION = 2048;  
}  // namespace

namespace traccc::stdpar {
namespace details {
/*
 * Structure that defines the start point of a partition and its size.
 */
struct ccl_partition {
  public:
    std::size_t start;
    std::size_t size;
    friend std::ostream& operator<< (std::ostream &out, ccl_partition const& p) {
      out << "[start: " << p.start << ", size: " << p.size <<  "]";
      return out;
    }
};

/*
 * Convenience structure to work with flattened data arrays instead of 
 * an array/vector of cells.
 */
struct cell_container {
    std::size_t size = 0;
    channel_id* channel0 = nullptr;
    channel_id* channel1 = nullptr;
    scalar* activation = nullptr;
    scalar* time = nullptr;
    geometry_id* module_id = nullptr;
};

/*
 * Convenience structure to work with flattened data arrays instead of 
 * an array/vector of measures.
 */
struct measurement_container {
    std::atomic<unsigned int> size = 0;
    scalar* channel0 = nullptr;
    scalar* channel1 = nullptr;
    scalar* variance0 = nullptr;
    scalar* variance1 = nullptr;
    geometry_id* module_id = nullptr;
};

namespace helper{
  /*
 * Helper class to print the partitions'
 */
void print_partitions(std::vector<details::ccl_partition> partitions){
  std::cout << "---------------------------" << std::endl;
  std::cout << "The current partitions" << std::endl;
  std::cout << "---------------------------" << std::endl;
  std::for_each_n(counting_iterator(0), partitions.size(), 
    [=](unsigned int i){;
    std::cout << "Partition " << i << std::endl;
    std::cout << partitions.at(i) << std::endl;
    std::cout << "--------" << std::endl;
  });
}

void print_cells_in_partitions(std::vector<details::ccl_partition> partitions, 
                               channel_id *c0, channel_id *c1){
  std::cout << "---------------------------" << std::endl;
  std::cout << "Cells per partition" << std::endl;
  std::cout << "---------------------------" << std::endl;
  std::for_each_n(counting_iterator(0), partitions.size(), 
    [=](unsigned int i){;
    std::cout << "Partition " << i << std::endl;
    std::cout << "[ ";
    std::for_each_n(counting_iterator(partitions[i].start), partitions[i].size, 
      [=](unsigned int j){
      std::cout << "(" << c0[j] << "," << c1[j] << ") " ;
    });
    std::cout << "] " << std::endl;
    std::cout << "--------" << std::endl;
  });
}

/*
 * Normally, something smart and beautiful could be done, but std::string is
 * not available on the device. Therefore, we write this differently.
 *
 * Therefore, we just assume that the array is never bigger than 1000 chars.
 * Otherwise, manual inspection will not serve anybody anyway.
 */
void print_array_on_device(char *array_label, index_t *array, int size, unsigned short partition_id){
  int max_buffer_size = 1000;
  int max_chars_per_element = 5; // an index_t has never more than 5 digits
  int n_elements = std::min(size, max_buffer_size/max_chars_per_element); 
  int spacer_width = 2; // corresponds to the length of ", "
  int buffer_size = n_elements*max_chars_per_element + n_elements*spacer_width; 

  char *buf = new char[buffer_size];
  char *pos = buf;
  for (int i = 0 ; i < n_elements ; i++) {
      if (i) {
          pos += sprintf(pos, ", ");
      }
      pos += sprintf(pos, "%hi", array[i]);
  }

  printf("Partition %d: %s = %s\n", partition_id, array_label, buf);

  delete[] buf;
}
} // end namespace details::helper


/*
 * Check if two cells are considered close enough to be part of the same cluster.
 */
inline bool is_adjacent(channel_id ac0, channel_id ac1, channel_id bc0,
                            channel_id bc1) {
    unsigned int p0 = (ac0 - bc0);
    unsigned int p1 = (ac1 - bc1);

    return p0 * p0 <= 1 && p1 * p1 <= 1;
}

inline void reduce_problem_cell(const cell_container& cells, index_t tid,
                                unsigned char *adjc, auto adjv) {
                                // Using auto for adjv because other definitions did not work.
                                // Open for suggestions to improve this.
  /*
   * The number of adjacent cells for each cell must start at zero, to
   * avoid uninitialized memory. adjv does not need to be zeroed, as
   * we will only access those values if adjc indicates that the value
   * is set.
   */
  adjc[tid] = 0;

  channel_id c0 = cells.channel0[tid];
  channel_id c1 = cells.channel1[tid];
  geometry_id gid = cells.module_id[tid];

  /*
   * First, we traverse the cells backwards, starting from the current
   * cell and working back to the first, collecting adjacent cells
   * along the way.
   */
  for (index_t j = tid - 1; j < tid; --j) {
      /*
        * Since the data is sorted, we can assume that if we see a cell
        * sufficiently far away in both directions, it becomes
        * impossible for that cell to ever be adjacent to this one.
        * This is a small optimisation.
        */
      if (cells.channel1[j] + 1 < c1 || cells.module_id[j] != gid) {
          break;
      }

      /*
        * If the cell examined is adjacent to the current cell, save it
        * in the current cell's adjacency set.
        */
      if (is_adjacent(c0, c1, cells.channel0[j], cells.channel1[j])) {
          adjv[tid][adjc[tid]++] = j;
      }
  }

  /*
    * Now we examine all the cells past the current one, using almost
    * the same logic as in the backwards pass.
    */
  for (index_t j = tid + 1; j < cells.size; ++j) {
      /*
        * Note that this check now looks in the opposite direction! An
        * important difference.
        */
      if (cells.channel1[j] > c1 + 1 || cells.module_id[j] != gid) {
          break;
      }

      if (is_adjacent(c0, c1, cells.channel0[j], cells.channel1[j])) {
          adjv[tid][adjc[tid]++] = j;
      }
  }
}

/*
 * Implementation of a simplified SV algorithm. 
 *
 * The implementation corresponds to Algorithm 1 of the following paper:
 * https://epubs.siam.org/doi/pdf/10.1137/1.9781611976137.5
 */
void simplified_sv(index_t* f, index_t* gf, unsigned char *adjc,
                          auto adjv, unsigned int size) {
  /*
   * The algorithm finishes if an iteration leaves the arrays unchanged.
   * This varible will be set if a change is made, and dictates if another
   * loop is necessary.
   */
  bool *gfc = new bool();

  do {
    /*
     * Reset the end-parameter to false, so we can set it to true if we
     * make a change to our arrays.
     */
    *gfc = false;

    /*
     * The algorithm executes in a loop of three distinct parallel
     * stages. In this first one, we examine adjacent cells and copy
     * their cluster ID if it is lower than our, essentially merging
     * the two together.
     */
    std::for_each_n(std::execution::unseq, counting_iterator(0), size, 
      [=](unsigned int i){
        for (unsigned char k = 0; k < adjc[i]; ++k) {
          index_t q = gf[adjv[i][k]];

          if (gf[i] > q) {
              f[f[i]] = q;
              f[i] = q;
          }
        }
    });

    /*
     * The second stage is shortcutting, which is an optimisation that
     * allows us to look at any shortcuts in the cluster IDs that we
     * can merge without adjacency information.
     */
    std::for_each_n(std::execution::unseq, counting_iterator(0), size, 
      [=](unsigned int i){
        if (f[i] > gf[i]) {
          f[i] = gf[i];
        }
    });

    /*
     * Update the array for the next generation, keeping track of any
     * changes we make.
     */
    std::for_each_n(std::execution::unseq, counting_iterator(0), size, 
      [=](unsigned int i){
        if (gf[i] != f[f[i]]) {
          gf[i] = f[f[i]];
          *gfc = true;
        }
    });

  } while (*gfc);

  /*
   * Mandatory clean up of memory on the heap.
   */
  delete gfc;
}

/*
 * Aggregate the information of all cells within a cluster to a single measurement.
 */
void aggregate_clusters(const cell_container &cells,
                        const measurement_container* out, index_t* f, unsigned int size) {
  /*
   * Keep track of the number of clusters already processed to properly index them in the output container.
   */
  unsigned int *cluster_index = new unsigned int(0);                               

  /*
   * Iterate over every cell in the partition and perform aggregation once per cluster.
   */
  std::for_each_n(counting_iterator(0), size, 
    [=](unsigned int i){
      /*
       * If and only if the value in the work arrays is equal to the index
       * of a cell, that cell is the "parent" of a cluster of cells. If
       * they are not, there is nothing for us to do. Easy!
       */
      if (f[i] == i) {
        /*
         * These variables keep track of the sums of X and Y coordinates
         * for the final coordinates, the total activation weight, as
         * well as the sum of squares of positions, which we use to
         * calculate the variance.
         */
        float sw = 0.0;
        float mx = 0.0, my = 0.0;
        float vx = 0.0, vy = 0.0;

        /*
         * Now, we iterate over all other cells to check if they belong
         * to our cluster. Note that we can start at the current index
         * because every child cell of a cluster is owned by a cell
         * with a higher ID.
         */
        for (index_t j = i; j < size; j++) {
            /*
             * If the value of this cell is equal to our, that means it
             * is part of our cluster. In that case, we take its values
             * for position and add them to our accumulators.
             */
            if (f[j] == i) {
                float w = cells.activation[j];

                sw += w;

                float pmx = mx, pmy = my;
                float dx = cells.channel0[j] - pmx;
                float dy = cells.channel1[j] - pmy;
                float wf = w / sw;

                mx = pmx + wf * dx;
                my = pmy + wf * dy;

                vx += w * dx * (cells.channel0[j] - mx);
                vy += w * dy * (cells.channel1[j] - my);
            }
        }

        /*
         * Write the average weighted x and y coordinates, as well as
         * the weighted average square position, to the output array.
         */
        out->channel0[*cluster_index] = mx;
        out->channel1[*cluster_index] = my;
        out->variance0[*cluster_index] = vx / sw;
        out->variance1[*cluster_index] = vy / sw;
        out->module_id[*cluster_index] = cells.module_id[i];

        /*
         * Increment the cluster index for the next iteration.
         */
        (*cluster_index)++;
      }
  });     
  
  /*
   * Mandatory clean up of memory on the heap.
   */
  delete cluster_index;
}

/*
 * Function that organized the parallel execution of the algotrithm.
 */
void fast_sv_kernel(
    const cell_container container, 
    const std::vector<details::ccl_partition> *partitions,
    measurement_container* _out_ctnr) {

  /*
   * If there is no partition, there is nothing to do for the kernel.
   */
  if(partitions->size() == 0) return;

  /*
   * Create a local data pointer to access the data in the array of the 
   * vector from within the parallel algorithm.
   */
  const details::ccl_partition *p = partitions->data();

  /*
   * Start running the work in different kernels using std par.
   */
  std::for_each_n(std::execution::par, counting_iterator(0), partitions->size(), // TODO: should be std::execution::par
    [=](unsigned int i){
    /*
     * Seek the correct cell region in the input data. Again, this is all a
     * contiguous block of memory for now, and we use the blocks array to
     * define the different ranges per block/module. At the end of this we
     * have the starting address of the block of cells dedicated to this
     * module, and we have its size.
     */
    const ccl_partition& partition = p[i];
    cell_container cells;
    cells.size = partition.size;
    cells.channel0 = &container.channel0[partition.start];
    cells.channel1 = &container.channel1[partition.start];
    cells.activation = &container.activation[partition.start];
    cells.time = &container.time[partition.start];
    cells.module_id = &container.module_id[partition.start];

    /*
     * As an optimisation, we will keep track of which cells are adjacent to
     * each other cell. To do this, we define, in thread-local memory or
     * registers, up to eight adjacent cell indices and we keep track of how
     * many adjacent cells there are (i.e. adjc[i] determines how many of
     * the eight values in adjv[i] are actually meaningful).
     */
    auto adjv = new index_t[cells.size][8];
    auto adjc = new unsigned char[cells.size];

    /*
     * This loop initializes the adjacency cache, which essentially
     * translates the sparse CCL problem into a graph CCL problem which we
     * can tackle with well-studied algorithms. 
     */
    std::for_each_n(counting_iterator(0), cells.size, 
      [=](unsigned int j){
        reduce_problem_cell(cells, j, adjc, adjv);
    });
    
    /*
     * These arrays are the meat of the pudding of this algorithm, and we
     * will constantly be writing and reading from them.
     */
    index_t *f = new index_t[cells.size];
    index_t *gf = new index_t[cells.size];

    /*
     * At the start, the values of f and gf should be equal to the ID of
     * the cell.
     */
    std::for_each_n(std::execution::unseq, counting_iterator(0), cells.size, 
      [=](unsigned int j){
        f[j] = j;
        gf[j] = j;
    });

    /*
     * Now we move onto the actual processing part.
     */
    details::simplified_sv(f, gf, adjc, adjv, cells.size);

    /*
     * Count the number of clusters by checking how many nodes have
     * themself assigned as a parent.
     */
    // TODO: Verify that adding an execution policy here makes actually sense. -> It would, but it is currently not supported by nvc++
    unsigned int number_of_clusters = std::count_if(std::execution::unseq, counting_iterator(0), counting_iterator(cells.size), 
      [=](unsigned int j){
        return f[j] == j;
    });

    /*
     * Add the number of clusters of each thread block to the total
     * number of clusters. At the same time, a cluster id is retrieved
     * for the next data processing step.
     * Note that this might be not the same cluster as has been treated
     * previously. However, since each thread block spawns a the maximum
     * amount of threads per block, this has no sever implications.
     */
    unsigned int outi = _out_ctnr->size.fetch_add(number_of_clusters);

    /*
     * Seek the correct output location for the results of the aggregated clusters.
     */
    measurement_container *out = new measurement_container();
    out->channel0 = (_out_ctnr->channel0)+outi;
    out->channel1 = (_out_ctnr->channel1)+outi;
    out->variance0 = (_out_ctnr->variance0)+outi;
    out->variance1 = (_out_ctnr->variance1)+outi;
    out->module_id = (_out_ctnr->module_id)+outi;

    /*
     * Perform the equivalent to the measurement creation, i.e. aggregating the clusters.
     */
    aggregate_clusters(cells, out, f, cells.size);

    /*
     * Mandatory clean up of temporary variables on the heap.
     */
    /* TODO: this causes an error when using std::exec::par or par_unseq. We accept the memory leak for the moment...
    for(unsigned int k = 0; k < MIN_CELLS_PER_PARTITION; k++){
      delete[] adjv[k];
    }
    */
    delete[] adjv;
    delete[] adjc;
    delete[] f;
    delete[] gf;
  });
}

/*
 * Split the problem into different partitions that will later be executed in parallel.
 *
 * This implementation uses already the flattened data compared to the CUDA version.
 */
std::vector<details::ccl_partition> partition(
    channel_id *channel1, geometry_id *module_id, std::size_t n_cells) {
    std::vector<details::ccl_partition> partitions;
    std::size_t index = 0;
    std::size_t size = 0;

    /*
     * Keep track of geometry id to create partition if modules change (if adequate)
     */
    geometry_id current_geometry = 0;

    /*
     * We start at 0 since this is the origin of the local coordinate 
     * system within a detector module.
     */
    channel_id last_mid = 0;

    /*
     * Iterate over every cell in the current data set.
     */
    for (std::size_t i = 0; i < n_cells; ++i) {
        /*
         * Create a new partition if an "empty" row is detected. A row 
         * is considered "empty" if the channel1 values between two 
         * consecutive cells have a difference > 1. Another criteria is if 
         * the cells of a new partition are starting to be treated.
         * To prevent creating many small partitions, the current partition
         * must have at least twice the size of threads per block. This 
         * guarantees that each thread handles later at least two cells.
         */
        if ((channel1[i] > last_mid + 1 || module_id[i] != current_geometry) && size >= MIN_CELLS_PER_PARTITION) {
            partitions.push_back(
                details::ccl_partition{.start = index, .size = size});

            index += size;
            size = 0;
        }

        current_geometry = module_id[i];
        last_mid = channel1[i];
        size += 1;
    }

    /*
     * Create the very last partition after having iterated over all cells.
     */
    if (size > 0) {
        partitions.push_back(
            details::ccl_partition{.start = index, .size = size});
    }

    return partitions;
}
}  // namespace details

host_measurement_container component_connection_ssv::operator()(
    const host_cell_container& data) const {
    // TODO: replace with call to host_cell_container.size() once code is working
    /*
     * Calculate the total amount of cells to deal with.
     */
    std::size_t total_cells = 0;

    for (std::size_t i = 0; i < data.size(); ++i) {
        total_cells += data.at(i).items.size();
    }

    /*
     * Flatten input data.
     */
    std::vector<channel_id> channel0;
    std::vector<channel_id> channel1;
    std::vector<scalar> activation;
    std::vector<scalar> time;
    std::vector<geometry_id> module_id;

    channel0.reserve(total_cells);
    channel1.reserve(total_cells);
    activation.reserve(total_cells);
    time.reserve(total_cells);
    module_id.reserve(total_cells);

    for (std::size_t i = 0; i < data.size(); ++i) {
        for (std::size_t j = 0; j < data.at(i).items.size(); ++j) {
            channel0.push_back(data.at(i).items.at(j).channel0);
            channel1.push_back(data.at(i).items.at(j).channel1);
            activation.push_back(data.at(i).items.at(j).activation);
            time.push_back(data.at(i).items.at(j).time);
            module_id.push_back(data.at(i).header.module);
        }
    }

    /*
    * Store the flattened arrays in a convenience data container.
    */
    details::cell_container container;
    container.size = total_cells;
    container.channel0 = channel0.data();
    container.channel1 = channel1.data();
    container.activation = activation.data();
    container.time = time.data();
    container.module_id = module_id.data();

    /*  
     * Run the partitioning algorithm sequentially.
     */
    std::vector<details::ccl_partition> partitions = 
      details::partition(container.channel1, container.module_id, container.size);
    
    /*
     * Reserve space for the result of the algorithm. Currently, there is 
     * enough space allocated that (in theory) each cell could be a single
     * cluster, but this should not be the case with real experiment data.
     */
    details::measurement_container* mctnr =
        new details::measurement_container();

    mctnr->channel0 = new scalar[total_cells];
    mctnr->channel1 = new scalar[total_cells];
    mctnr->variance0 = new scalar[total_cells];
    mctnr->variance1 = new scalar[total_cells];
    mctnr->module_id = new geometry_id[total_cells];

    /*
     * Run the connected component labeling algorithm to retrieve the clusters.
     */
    fast_sv_kernel(container, &partitions, mctnr);

    /*
     * Transform flat data structure to expected output format again.
     */
    host_measurement_container out;

    for (std::size_t i = 0; i < data.size(); ++i) {
        vecmem::vector<measurement> v;

        for (std::size_t j = 0; j < mctnr->size; ++j) {
            if (mctnr->module_id[j] == data.at(i).header.module) {
                measurement m;

                m.local = {mctnr->channel0[j], mctnr->channel1[j]};
                m.variance = {mctnr->variance0[j], mctnr->variance1[j]};

                v.push_back(m);
            }
        }

        out.push_back(cell_module(data.at(i).header), std::move(v));
    }

    return out;
}
}  // namespace traccc::stdpar
