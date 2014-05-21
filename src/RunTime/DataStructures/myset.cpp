#include "stdio.h"
#include "DataStructures/myset.hpp"

using namespace std;

int ordered_vec::insert(const int newnum)
{
  int min = 0, max = elements.size()-1;
  while( min <= max )
    {
      // printf("min = %d, max = %d, mid = %d\n",min,max,mid);
      int mid = (min+max)/2;
      if( elements[mid] == newnum ) break;
      else if( newnum > elements[mid] )
	min = mid + 1;
      else
	max = mid - 1;
    }
  //printf("min = %d, max = %d\n",min,max);
  if( min > max )
    {
      if( min == elements.size() )
	elements.push_back(newnum);
      else
	{
	  deque<int>::iterator pos = elements.begin() + min;
	  elements.insert(pos,newnum);
	}
      return min;
    }
  else
    return -1;
}

int ordered_vec::posn(const int value) const
{
  int min = 0 , max= elements.size() - 1, mid ;
  if( elements.size() == 0 )
    return 0;

  assert( value >= elements[min] && value <= elements[max] );
  while( min <= max )
    {
      mid = (min+max)/2;
      // printf("min = %d, max = %d, mid = %d\n",min,max,mid);
      if( elements[mid] == value ) break;
      else if( value  > elements[mid] )
	min = mid + 1;
      else
	max = mid - 1;      
    }

  // printf("min = %d, max = %d, mid = %d\n",min,max,mid);
  if( min > max ) 
    return ( min + max ) / 2;
  else
    return mid;

}

//Dummy main function used only for testing
// int main(int *argc , char **argv)
// {
//   ordered_vec test;
//   deque<int> help;
//   // test.insert(30);
//   // for( int i = 0; i < test.size() ; i++ )
//   //   printf("%d ",test[i]);
//   // printf("\n");
//   for( int i =0 ; i < 20 ; i+=3 )
//     help.push_back(i);
//   test.set_vals(help.begin(),help.end());
//   for( int i = 0 ; i < test.size() ; i++ )
//      printf("%d ",test[i]);
//   printf("\n");
//   test.insert(atoi(argv[1]));
//   test.insert(atoi(argv[1])+4);
//   for( int i = 0 ; i < test.size() ; i++ )
//     printf("%d ",test[i]);
//   printf("\n");
  
//   int posn = test.posn(atoi(argv[2]));
//   printf("posn = %d\n",posn);
//   printf("Found %d at posn = %d\n",atoi(argv[2]),test.find(atoi(argv[2])));
//   test.erase_position(posn - 3 );
//   test.erase(atoi(argv[2]));
//   posn = test.posn(atoi(argv[2]));
//   printf("Found %d at posn = %d\n",atoi(argv[2]),test.find(atoi(argv[2])));
  
//   for( int i = 0 ; i < test.size() ; i++ )
//     printf("%d ",test[i]);
//   printf("\n");
  
  

//   return 0;
// }
  
	  
      
