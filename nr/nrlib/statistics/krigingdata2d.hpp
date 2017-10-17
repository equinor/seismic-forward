#ifndef NRLIB_STATISTICS_KRIGINGDATA2D_HPP
#define NRLIB_STATISTICS_KRIGINGDATA2D_HPP

#include <vector>

namespace NRLib {
  class KrigingData2D
  {
  public:
    //KrigingData2D();
    virtual ~KrigingData2D(void) {};

    virtual KrigingData2D*  Clone() const = 0;

    virtual KrigingData2D*  NewInstance() const = 0;

    virtual int             GetNumberOfData(void)        const = 0;

    virtual size_t          CountDataInBlock(size_t imin, size_t imax, size_t jmin, size_t jmax) const = 0;
    virtual void            AddDataToBlock(KrigingData2D & data_block, size_t imin, size_t imax, size_t jmin, size_t jmax) const = 0;

    virtual void            WriteToFile(const std::string & name) const = 0;
  };
}
#endif
