template <typename T, bool lazy, bool alwaysCheck , bool assertion, bool eigenAlignedNew>
inline CheckedItem<T,lazy,alwaysCheck,assertion , eigenAlignedNew>::CheckedItem():
  isSet_(false)
{
}

template <typename T, bool lazy, bool alwaysCheck , bool assertion, bool eigenAlignedNew>
inline CheckedItem<T,lazy,alwaysCheck,assertion , eigenAlignedNew>::CheckedItem(const CheckedItem & c):
  v_(c.v_)
{
  if (do_check_)
  {
    isSet_ = c.isSet_;
  }

  if (do_assert_)
  {
    assertMsg_ = c.assertMsg_;
  }

  if (do_exception_)
  {
    exceptionPtr_=c.exceptionPtr_;
  }
}



template <typename T, bool lazy, bool alwaysCheck , bool assertion, bool eigenAlignedNew>
inline CheckedItem<T,lazy,alwaysCheck,assertion, eigenAlignedNew>::CheckedItem(const T& v):
  isSet_(true), v_(v)
{
}

template <typename T, bool lazy, bool alwaysCheck , bool assertion, bool eigenAlignedNew>
inline T& CheckedItem<T,lazy,alwaysCheck , assertion,eigenAlignedNew>::operator=(const T& v)
{
  isSet_.set(true);
  return v_=v;
}

template <typename T, bool lazy, bool alwaysCheck , bool assertion, bool eigenAlignedNew>
inline CheckedItem<T,lazy,alwaysCheck , assertion,eigenAlignedNew> &
    CheckedItem<T,lazy,alwaysCheck , assertion,eigenAlignedNew>::operator=(const CheckedItem & c)
{
  v_=c.v_;

  if (do_check_)
  {
    isSet_ = c.isSet_;
  }

  if (do_assert_)
  {
    assertMsg_ = c.assertMsg_;
  }

  if (do_exception_)
  {
    exceptionPtr_=c.exceptionPtr_;
  }

  return *this;
}

template <typename T, bool lazy, bool alwaysCheck , bool assertion, bool eigenAlignedNew>
inline  CheckedItem<T,lazy,alwaysCheck , assertion,eigenAlignedNew>::operator T() const
{
  return chckitm_getRef();
}


template <typename T, bool lazy, bool alwaysCheck , bool assertion, bool eigenAlignedNew>
inline  T CheckedItem<T,lazy,alwaysCheck , assertion,eigenAlignedNew>::chckitm_getValue() const
{
  return chckitm_getRef();
}

template <typename T, bool lazy, bool alwaysCheck , bool assertion, bool eigenAlignedNew>
inline  const T& CheckedItem<T,lazy,alwaysCheck , assertion,eigenAlignedNew>::chckitm_getRef() const
{
  chckitm_check_();
  return v_;
}

template <typename T, bool lazy, bool alwaysCheck , bool assertion, bool eigenAlignedNew>
inline bool CheckedItem<T,lazy,alwaysCheck , assertion,eigenAlignedNew>::chckitm_check_() const throw(std::exception)
{
  if (assertion)
  {
    BOOST_ASSERT_MSG(chckitm_isSet(),assertMsg_.get());
  }

  if (alwaysCheck || isDebug)
  {
    if (!chckitm_isSet())
    {
      throw (*(exceptionPtr_.get()));
    }
  }
  return (chckitm_isSet());
}

template <typename T, bool lazy, bool alwaysCheck , bool assertion, bool eigenAlignedNew>
inline bool CheckedItem<T,lazy,alwaysCheck , assertion,eigenAlignedNew>::chckitm_isSet() const
{
  return isSet_.get();
}

template <typename T, bool lazy, bool alwaysCheck , bool assertion, bool eigenAlignedNew>
inline void CheckedItem<T,lazy,alwaysCheck , assertion,eigenAlignedNew>::chckitm_reset()
{
  isSet_.set(false);
}

template <typename T, bool lazy, bool alwaysCheck , bool assertion, bool eigenAlignedNew>
inline void CheckedItem<T,lazy,alwaysCheck , assertion,eigenAlignedNew>::chckitm_set(bool value)
{
  isSet_.set(value);
}


template <typename T, bool lazy, bool alwaysCheck , bool assertion, bool eigenAlignedNew>
inline void CheckedItem<T,lazy,alwaysCheck , assertion,eigenAlignedNew>::chckitm_setAssertMessage(std::string s)
{
  assertMsg_.set(s);
}

template <typename T, bool lazy, bool alwaysCheck , bool assertion, bool eigenAlignedNew>
inline void CheckedItem<T,lazy,alwaysCheck , assertion,eigenAlignedNew>::chckitm_setExceptionPtr(std::exception* e)
{
  exceptionPtr_.set(e);
}





template <typename MatrixType, bool lazy>
inline IndexedMatrixT<MatrixType,lazy>::IndexedMatrixT(const MatrixType& v,TimeIndex k):
  IsSet(true),
  k_(k),
  v_(v)
{
}

template <typename MatrixType, bool lazy>
inline IndexedMatrixT<MatrixType,lazy>::IndexedMatrixT():
  IsSet(false),
  k_(0)
{
}

///Says whether the matrix is initialized or not
template <typename MatrixType, bool lazy>
inline bool IndexedMatrixT<MatrixType,lazy>::isSet()const
{
    return (IsSet::get());
}

///Set the value of the matrix and the time sample
template <typename MatrixType, bool lazy>
inline MatrixType IndexedMatrixT<MatrixType,lazy>::set(const MatrixType& v,TimeIndex k)
{
  IsSet::set(true);
  k_=k;
  v_=v;
  return(v);
}

///Set the value of the matrix and the time sample
template <typename MatrixType, bool lazy>
inline void IndexedMatrixT<MatrixType,lazy>::reset()
{
  IsSet::set(false);
}



///Checks whether the matrix is set or not (assert)
///does nothing in release mode
template <typename MatrixType, bool lazy>
inline bool IndexedMatrixT<MatrixType,lazy>::check_() const
{
    BOOST_ASSERT(isSet() && "Error: Matrix not initialized, if you are initializing it, \
                            use set() function.");
    return isSet();
}




template <typename MatrixType, bool lazy>
inline void IndexedMatrixT<MatrixType,lazy>::setIndex(TimeIndex k)
{
  check_();
  k_=k;
}

///Get the matrix value
template <typename MatrixType, bool lazy>
inline  const MatrixType & IndexedMatrixT<MatrixType,lazy>::operator()()const
{
  check_();
  return v_;
}

///Get the matrix value
template <typename MatrixType, bool lazy>
inline MatrixType & IndexedMatrixT<MatrixType,lazy>::operator()()
{
  check_();
  return v_;
}

///Get the time index
template <typename MatrixType, bool lazy>
TimeIndex IndexedMatrixT<MatrixType,lazy>::getTime()const
{
  check_();
  return k_;
}





///Set the value of the matrix and the time sample
template <typename MatrixType>
inline void IndexedMatrixArrayT<MatrixType>::setValue(const MatrixType& v,TimeIndex k)
{
    if (checkIndex(k))
    {
        (*this)[k]=v;
    }
    else
    {
        checkNext_(k);
        if (v_.size()==0)
            k_=k;

        v_.push_back(v);
    }
}

template <typename MatrixType>
inline void IndexedMatrixArrayT<MatrixType>::pushBack(const MatrixType& v)
{
    v_.push_back(v);
}

template <typename MatrixType>
inline void IndexedMatrixArrayT<MatrixType>::popFront()
{
    check_();
    v_.pop_front();
    ++k_;
}

///Get the matrix value
template <typename MatrixType>
inline MatrixType IndexedMatrixArrayT<MatrixType>::operator[](TimeIndex time)const
{
    check_(time);
    return v_[time - k_];
}

///Get the matrix value
template <typename MatrixType>
inline MatrixType & IndexedMatrixArrayT<MatrixType>::operator[](TimeIndex time)
{
    check_(time);
    return v_[time - k_];
}


///gets the first value
template <typename MatrixType>
inline const MatrixType & IndexedMatrixArrayT<MatrixType>::front() const
{
    return v_.front();
}

///gets the first value
template <typename MatrixType>
inline MatrixType& IndexedMatrixArrayT<MatrixType>::front()
{
    return v_.front();
}

///gets the last value
template <typename MatrixType>
inline const MatrixType & IndexedMatrixArrayT<MatrixType>::back() const
{
    return v_.back();
}

///gets the last value
template <typename MatrixType>
inline MatrixType & IndexedMatrixArrayT<MatrixType>::back()
{
    return v_.back();
}

///Get the time index
template <typename MatrixType>
inline long int IndexedMatrixArrayT<MatrixType>::getLastIndex()const
{
  return long(k_+v_.size())-1;
}

///Get the time index
template <typename MatrixType>
inline TimeIndex IndexedMatrixArrayT<MatrixType>::getNextIndex()const
{
  return k_+v_.size();
}


///Get the time index
template <typename MatrixType>
inline TimeIndex IndexedMatrixArrayT<MatrixType>::getFirstIndex()const
{
  return k_;
}

template <typename MatrixType>
inline TimeIndex IndexedMatrixArrayT<MatrixType>::setLastIndex(int index)
{
  return k_=index-(v_.size()+1);
}

template <typename MatrixType>
inline TimeIndex IndexedMatrixArrayT<MatrixType>::setFirstIndex(int index)
{
  return k_=index;
}

template <typename MatrixType>
inline TimeSize IndexedMatrixArrayT<MatrixType>::size() const
{
    return v_.size();
}

///Switch off the initialization flag, the value is no longer accessible
template <typename MatrixType>
inline void IndexedMatrixArrayT<MatrixType>::reset()
{
    k_=0;
    v_.clear();
}

template <typename MatrixType>
void IndexedMatrixArrayT<MatrixType>::clear()
{
  k_=k_+v_.size();
  v_.clear();
}

template <typename MatrixType>
inline bool IndexedMatrixArrayT<MatrixType>::checkIndex(TimeIndex time) const
{
    return (v_.size()>0 &&
            k_<=time &&
            TimeIndex(k_+v_.size()) > time);
}

///Checks whether the matrix is set or not (assert)
///does nothing in release mode
template <typename MatrixType>
inline void IndexedMatrixArrayT<MatrixType>::check_(TimeIndex time)const
{
    (void)time;//avoid warning in release mode
    BOOST_ASSERT(checkIndex(time) && "Error: Time out of range");
}

///Checks whether the matrix is set or not (assert)
///does nothing in release mode
template <typename MatrixType>
inline void IndexedMatrixArrayT<MatrixType>::check_()const
{
    BOOST_ASSERT(v_.size() && "Error: Matrix array is empty");
}

template <typename MatrixType>
inline void IndexedMatrixArrayT<MatrixType>::checkNext_(TimeIndex time)const
{
    (void)time;//avoid warning
    BOOST_ASSERT( (v_.size()==0 || TimeIndex( k_+v_.size()) == time )&&
                  "Error: New time instants must be consecutive to existing ones");
}

///resizes the array
template <typename MatrixType>
inline void IndexedMatrixArrayT<MatrixType>::resize(TimeSize i, const MatrixType & m )
{
    v_.resize(i,m);
}

///Default constructor
template <typename MatrixType>
IndexedMatrixArrayT<MatrixType>::IndexedMatrixArrayT():
  k_(0)
{
}

template <typename MatrixType>
typename IndexedMatrixArrayT<MatrixType>::Array IndexedMatrixArrayT<MatrixType>::getArray() const
{
  Array v;

  for (TimeSize i=0; i<v_.size(); ++i)
  {
    v.push_back(v_[i]);
  }

  return v;
}

template <typename MatrixType>
void IndexedMatrixArrayT<MatrixType>::truncateAfter(TimeIndex time)
{
  if (v_.size()>0)
  {
    if (time >= getFirstIndex())
    {
      if (time < getLastIndex())
      {
        resize (time-getFirstIndex()+1);
      }
    }
    else
    {
      v_.clear();
    }
  }
}

template <typename MatrixType>
void IndexedMatrixArrayT<MatrixType>::truncateBefore(TimeIndex time)
{
  if (v_.size()>0)
  {
    if (time < getLastIndex())
    {
      for (TimeIndex i=getFirstIndex(); i<time ; ++i)
      {
        v_.pop_front();
      }

      setFirstIndex(time);
    }
    else
    {
      v_.clear();
    }
  }
}

template <typename MatrixType>
void IndexedMatrixArrayT<MatrixType>::readFromFile(const std::string & filename, size_t rows, size_t cols, bool withTimeStamp)
{
  readFromFile(filename.c_str(),rows, cols, withTimeStamp);
}

template <typename MatrixType>
void IndexedMatrixArrayT<MatrixType>::readFromFile(const char * filename, size_t rows, size_t cols, bool withTimeStamp)
{
  reset();

  std::ifstream f;

  f.open(filename);

  if (f.is_open())
  {

    Matrix m(Matrix::Zero(rows,cols));

    bool continuation=true;
    int k=0;

    while (continuation)
    {


      if (withTimeStamp)
      {
        f >> k;
      }


      if (f.fail())
        continuation=false;
      else
      {
        for (size_t i = 0 ; i<rows; ++i)
        {
          for (size_t j = 0 ; j<cols; ++j)
          {
            f >> m(i,j);
          }
        }

        setValue(m,k);
        ++k;
      }
    }
  }
}

template <typename MatrixType>
void IndexedMatrixArrayT<MatrixType>::readVectorsFromFile(const std::string & filename, bool withTimeStamp )
{
  readVectorsFromFile(filename.c_str(),withTimeStamp);
}

template <typename MatrixType>
void IndexedMatrixArrayT<MatrixType>::readVectorsFromFile(const char * filename, bool withTimeStamp )
{
  reset();

  std::ifstream f;

  f.open(filename);

  if (f.is_open())
  {
    std::string s;
    Vector v;
    int k=0;

    bool continuation=true;

    while (continuation)
    {
      std::getline(f,s);

      std::stringstream ss(s);

      if (withTimeStamp)
      {
        ss >> k;
      }



      if (f.fail())
        continuation=false;
      else
      {
        std::vector<double> doublecontainer;
        double component;
        bool readingVector = true;
        while (readingVector)
        {
          ss>> component;

          if (ss.fail())
          {
            readingVector=false;
          }
          else
          {
            doublecontainer.push_back(component);
          }
        }
        v.resize(doublecontainer.size());
        for (unsigned i=0 ; i<doublecontainer.size() ; ++i)
        {
          v(i)=doublecontainer[i];
        }
        setValue(v,k);
        ++k;

      }
    }
  }
}

template <typename MatrixType>
void IndexedMatrixArrayT<MatrixType>::writeInFile(const std::string & filename,  bool clearLog, bool append)
{
  writeInFile(filename.c_str(),  clearLog, append);
}


template <typename MatrixType>
void IndexedMatrixArrayT<MatrixType>::writeInFile(const char * filename, bool clearLog, bool append)
{
  std::ofstream f;
  if (!append)
  {
    f.open(filename);
  }
  else
  {
    f.open(filename,std::ofstream::app);
  }

  if (f.is_open())
  {
    if (size()>0)
    {

      for (TimeIndex k=getFirstIndex(); k<getNextIndex(); ++k)
      {

        f << k;

        MatrixType & m = operator[](k);

        for (int i = 0 ; i< m.rows(); ++i)
        {
          for (int j = 0 ; j< m.cols(); ++j)
          {
            f << " "<< m(i,j);
          }
        }
        f << std::endl;
      }
    }

    if (clearLog)
    {
      clear();
    }
  }
  else
  {
    std::stringstream ss;
    ss<< "Logger: File " <<filename<<" could not be created/opened.";
    std::runtime_error e(ss.str().c_str());
    throw e;
  }

}

namespace tools
{

  inline void SimplestStopwatch::start()
  {
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
  }

  inline double SimplestStopwatch::stop()
  {
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time3);

    return diff(time2,time3)-diff(time1,time2);
  }

  double SimplestStopwatch::diff(const timespec & start, const timespec & end)
  {
    return 1e9*double(end.tv_sec-start.tv_sec) + double(end.tv_nsec - start.tv_nsec);
  }
}



