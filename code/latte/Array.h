#ifndef ARRAY__H
#define ARRAY__H

template <class T>
class BigArray
{
 public:
  int length;
  T * elements;
  BigArray(int Length){
    length = Length;
    elements = new T[length];
  }

  BigArray( const BigArray & old ){
    delete [] elements;
    length = old.length;
    elements = new T[length];
    for(int i = 0; i < length; i++)
      elements[i] = old.elements[i];
  }
  
  ~BigArray()
    { delete [] elements; }
  
  T & operator[](int index)
    {
      return elements[index];
    }
  
};

#endif



