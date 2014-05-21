#ifndef __MYSET_HPP__
#define __MYSET_HPP__

#include "stdio.h"
#include <deque>
//#include <vector>
#include <cassert>
#include <cstdlib>

//Data structure which maintains elements in increasing order
class ordered_vec{

 private:

  std::deque<int> elements; //THe list of elements
  
 public:
  
  //Function to return a value at given index
  inline int operator[](const int i) { return elements[i];}
  
  //Function to insert the given element
  //Output : Position at which it is inserted. If element already exists , returns -1
  int insert(const int); 
  
  //Function to return the index which contains the highest element less than or equal to the input
  int posn(const int) const;

  //Function to find the given value in the list of elements, -1 if it is not in the list
  inline int find(const int val) const
  {
    if( elements.size() == 0 || val > elements[elements.size()-1] || val < elements[0] )
      return -1;
    int temp = posn(val);
    if( val == elements[temp] )
  	return temp;
    else
    	return -1;
  }

  //Function to return the size of the list
  inline int size() const{ return (int)elements.size(); }

  //Function to return the iterator to the start of the list
  inline std::deque<int>::const_iterator begin() const{ return elements.begin(); }

  //Function to return the iterator to the end of the list
  inline std::deque<int>::const_iterator end() const{ return elements.end(); }

  //Function to return the value of an element at the input position
  inline int get_value(const int k) const
  { 
    if( k < elements.size() )
      return elements[k]; 
    else
      {
	assert(0); 
	return -1;
      }
  }
    
  //Function to erase the given value if it is in the list
  inline void erase(const int val)
  {
    int position = posn(val);
    if( elements[position] == val )
      elements.erase(elements.begin() + position );
  }

  //Function to erase the value at a given position in the list
  inline void erase_position(const int posn) { assert( posn >= 0 ) ; elements.erase(elements.begin() + posn ); }
};

#endif
