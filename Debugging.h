/**
 * \file Debugging.h
 *
 * \brief Debugging utilities.
 *
 * \author Bill White
 * \version 1.0
 *
 * Contact: bill.c.white@gmail.com
 * Created on: 8/q/11
 */

#include <iostream>
#include <vector>

#ifndef DEBUGGING_H
#define	DEBUGGING_H

using namespace std;

/***************************************************************************//**
 * Print a vector of T values with optional title.
 * \param [in] vec vector of T type values
 * \param [in] title optional title to print before the vector
 ******************************************************************************/
template <class T> void PrintVector(std::vector<T> vec, std::string title="");

template <class T> void PrintVector(vector<T> vec, string title) {
  if(title != "") {
    cout << title << ": ";
  }
  cout << "[ ";
  copy(vec.begin(), vec.end(), ostream_iterator<T>(cout, " "));
  cout << "]" << endl;
}

#endif	/* DEBUGGING_H */

