#ifndef MATRIX_H
#define MATRIX_H

class iMatrix {
 public:
   iMatrix(); //empty constructor
   iMatrix(unsigned rows, unsigned cols);
   int& operator() (unsigned row, unsigned col);
   int operator() (unsigned row, unsigned col) const;  

   void resize(unsigned rows, unsigned cols);

  ~iMatrix();                              // Destructor
//   iMatrix(const iMatrix& m);               // Copy constructor
   iMatrix& operator= (const iMatrix& m);   // Assignment operator

 private:
   unsigned rows_, cols_;
   int* data_;
 };

//---------------------------------------------------------
// this is matrix.cpp. add: #include "matrix.h" if seperate file

 inline
 iMatrix::iMatrix(){
   rows_=0, cols_=0; //empty
 }

 inline
 iMatrix::iMatrix(unsigned rows, unsigned cols)
   : rows_ (rows), cols_ (cols)
   //data_ <--initialized below (after the 'if/throw' statement)
 {
   if (rows == 0 || cols == 0)
     cout<<"iMatrix constructor has 0 size"<<endl;
   data_ = new int[rows * cols];
 }

 inline
 void iMatrix::resize(unsigned rows, unsigned cols){
   if (rows == 0 || cols == 0)
     cout<<"iMatrix constructor has 0 size"<<endl;
  if (rows_ != 0 || cols_ != 0)
    delete[] data_;
  rows_=rows; 
  cols_=cols; 
  data_ = new int[rows * cols];
 }

 inline
 iMatrix::~iMatrix()
 {
   if (rows_ != 0 || cols_ != 0)
    delete[] data_;
 }
 
 inline
 int& iMatrix::operator() (unsigned row, unsigned col)
 {
   if (row >= rows_ || col >= cols_)
     cout<<"iMatrix subscript out of bounds"<<endl;
   return data_[cols_*row + col];
 }
 
 inline
 int iMatrix::operator() (unsigned row, unsigned col) const
 {
   if (row >= rows_ || col >= cols_)
     cout<<"const iMatrix subscript out of bounds"<<endl;
   return data_[cols_*row + col];
 }

#endif
