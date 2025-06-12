
#include "mfem.hpp"
#include <functional>

using namespace std;
using namespace mfem;



void dist_tag_func(const Vector & x_dist, Array<int> & tags,  std::function<int(const real_t)> dist_tag){
  if(tags.Size() > 0){ //Only apply if these array(s) are populated
    const bool use_dev = tags.UseDevice() || x_dist.UseDevice();
    const int n = x_dist.Size();
    // Use read+write access for X - we only modify some of its entries
    auto d_X = x_dist.Read(use_dev);
    auto d_tags = tags.ReadWrite(use_dev);
    mfem::forall_switch(use_dev, n, [=] MFEM_HOST_DEVICE (int i)
    {
      d_tags[i] = dist_tag(d_y[i]);
    });
  }
};
