#ifndef __GSTL_M_SYMETRIC_MATRIX_H__
#define __GSTL_M_SYMETRIC_MATRIX_H__


#include <GsTL/utils/gstl_error_messages.h>

#include <iostream>


namespace GSTL_TNT_LZL {


template <class T>
class Symmetric_matrix {

 public:

  typedef T element_type;
  typedef T value_type;
  typedef T& reference;
  typedef const T& const_reference;

  Symmetric_matrix() {
    N_ =0; 
    size_ =0;
    val_ = 0;
    mat_=0;
  }

  Symmetric_matrix(unsigned long long N);
  Symmetric_matrix(unsigned long long N, reference val);

  ~Symmetric_matrix();

  Symmetric_matrix& resize(unsigned long long N, unsigned long long M);

  inline unsigned long long type() const {return 1;}

  inline unsigned long long num_rows() const;
  inline unsigned long long num_cols() const;
  inline unsigned long long dim() const {return N_;}
  inline reference operator()(unsigned long long i, unsigned long long j);
  inline const_reference operator()(unsigned long long i, unsigned long long j) const;

  inline const	unsigned long long& pt_inx(unsigned long long i, unsigned long long j) const { unsigned long long a = 0;  return a;};
  inline			unsigned long long& pt_inx(unsigned long long i, unsigned long long j)  {unsigned long long a = 0;return a;};;
	  inline const	T& distance(unsigned long long i, unsigned long long j) const {T a = 0;return a;};;
		  inline			T& distance(unsigned long long i, unsigned long long j)  {T a = 0;return a;};;
  
  // This is an awful function, but it can be useful
  inline T* raw_access() {return val_;}


  inline unsigned long long acutal_size() const { return size_;}

  void destroy();

  int Read(char *filename);
  int Write(char *filename);

 private:
  unsigned long long N_;
  unsigned long long size_;
  value_type* val_;
  value_type** mat_;

  void initialize(unsigned long long N);
};  // end class gstl_symetric_matrix





template<class T>
void Symmetric_matrix<T>::initialize(unsigned long long N) {
  N_ = N;
  size_= (N*N + N) / 2;

  val_ = new T[size_];
  mat_ = new T*[N];

  for(unsigned long long i=0; i<size_; i++)
    val_[i] = i;

  T* p = val_;
  for(unsigned long long j=1; j<=N; j++) {
    mat_[j-1] = p;
    p += j;
  }
}


template<class T>
void Symmetric_matrix<T>::destroy() {
  if(val_ == 0 ) return ;

  delete [] mat_;
  delete [] val_;

  mat_ = 0;
  val_ = 0;

  N_ =0; 
  size_ =0;
  val_ = 0;
  mat_=0;
}


template<class T>
Symmetric_matrix<T>::Symmetric_matrix(unsigned long long N) {
  initialize(N);
}

template<class T>
Symmetric_matrix<T>::Symmetric_matrix(unsigned long long N, reference val) {
  initialize(N);

  for(unsigned long long i=0; i<size_; i++)
    val_[i] = val;
}


template<class T>
Symmetric_matrix<T>::~Symmetric_matrix() {
  destroy();
}


template<class T>
Symmetric_matrix<T>& Symmetric_matrix<T>::resize(unsigned long long N, unsigned long long) {
  if (N == N_)
    return *this;
  
  destroy();
  initialize(N);
  return *this;
}

template<class T>
inline unsigned long long Symmetric_matrix<T>::num_rows() const {
  return N_;
}

template<class T>
inline unsigned long long Symmetric_matrix<T>::num_cols() const {
  return N_;
}


template<class T>
inline T& Symmetric_matrix<T>::operator()(unsigned long long i, unsigned long long j) {
  gstl_assert(i>0 && i<=N_);
  gstl_assert(j>0 && j<=N_);
  
  if(i>j)
    return mat_[i-1][j-1];
  
  return mat_[j-1][i-1];
}

template<class T>
inline const T& Symmetric_matrix<T>::operator()(unsigned long long i, unsigned long long j) const {
  if(i>j)
    return mat_[i-1][j-1];
  
  return mat_[j-1][i-1];
}
  

template<class T>
int Symmetric_matrix<T>::Write(char *filename)
{
	unsigned long long  nRows = this->num_rows();
	unsigned long long  nCols = this->num_cols();

	FILE *fp = fopen(filename, "wb");

	fwrite(&nRows, sizeof(unsigned long long), 1, fp);	
	fwrite(&nCols, sizeof(unsigned long long), 1, fp);	

	fwrite(this->raw_access(), sizeof(T), this->acutal_size(), fp);	

	return fclose(fp);
}

template<class T>
int Symmetric_matrix<T>::Read(char *filename)
{
	unsigned long long  nRows = 0, nCols = 0;

	FILE *fp = fopen(filename, "rb");

	fread(&nRows, sizeof(unsigned long long), 1, fp);	
	fread(&nCols, sizeof(unsigned long long), 1, fp);	

	this->resize(nRows, nCols);

	fread(this->raw_access(), sizeof(T), this->acutal_size(), fp);

	return fclose(fp);
}

/* ***************************  I/O  ********************************/

template <class T>
std::ostream& operator<<(std::ostream& s, const Symmetric_matrix<T>& A)
{
    unsigned long long M=A.num_rows();

    s << M << " " << M << "\n";

    for (unsigned long long i=1; i<=M; i++)
    {
        for (unsigned long long j=1; j<=M; j++)
        {
            s << A(i,j) << " ";
        }
        s << "\n";
    }

    return s;
}

template <class T>
std::istream& operator>>(std::istream& s, Symmetric_matrix<T>& A)
{
	//unsigned long long  nRows = A.num_rows();
	//unsigned long long  nCols = A.num_cols();

	unsigned long long  nRows = 0;
	unsigned long long  nCols = 0;

	s >> nRows >> nCols ;

	A.resize(nRows, nCols);

	for (unsigned long long i=1; i<= nRows; i++)
	{
		for (unsigned long long j=1; j<= nCols; j++)
		{
			s>>A(i,j);
		}
	}

	return s;
}

}
#endif