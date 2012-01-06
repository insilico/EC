/**
 * \file best_n.h
 *
 * \brief Find the best n keeping original order for ties - stable sort.
 *
 * \author Nate Barney
 * \version 1.0
 *
 * Contact: bill.c.white@gmail.com
 * Created on: 4/7/04
 */

#include <vector>
#include <algorithm>
#include <boost/pointee.hpp>

namespace insilico
{

  /***************************************************************************//**
   * Get the best n values with ties keeping same original order.
   * \param [in] begin iterator of the beginning of a input container
   * \param [in] end iterator of the end of a input container
   * \param [out] out iterator of the beginning of a output container
   * \param [in] size best n value
   * \param [in] comp compare functor
   * \return path/filename without extension
   ******************************************************************************/
  template <typename InputIt, typename OutputIt, typename Comp>
  void best_n(InputIt begin, InputIt end, OutputIt out, size_t n, Comp comp) {
    typedef typename boost::pointee<InputIt>::type T;
    typename std::vector<T> best;
    std::vector<size_t> indices(n);
    int maxindex = 0;
    int i = 0;

    best.reserve(n);

    for(typename std::vector<T>::const_iterator it = begin;
        it != end; ++it, ++i) {
      if(best.size() < n) {
        indices[maxindex] = i;
        best.push_back(*it);

        if(best.size() == n)
          maxindex = std::distance(best.begin(),
                                   std::max_element(best.begin(), best.end(), comp));
        else
          ++maxindex;

        continue;
      }

      if(comp(*it, best[maxindex])) {
        indices[maxindex] = i;
        best[maxindex] = *it;
        maxindex = std::distance(best.begin(),
                                 std::max_element(best.begin(), best.end(), comp));
      }
    }

    for(typename std::vector<T>::iterator i = best.begin();
        i != best.end(); ++i)
      *out++ = *i;
  }

}
