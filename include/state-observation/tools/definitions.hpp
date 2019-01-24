/**
 * \file     definitions.hpp
 * \author   Mehdi Benallegue
 * \date     2013
 * \brief    Definitions of types and some structures.
 *
 * \details
 *
 *
 */

#ifndef STATEOBSERVATIONDEFINITIONSHPP
#define STATEOBSERVATIONDEFINITIONSHPP

//#define STATEOBSERVATION_VERBOUS_CONSTRUCTORS

#include <vector>
#include <deque>
#include <stdexcept>

#ifdef STATEOBSERVATION_VERBOUS_CONSTRUCTORS
#   include <iostream>
#endif

#include <boost/assert.hpp>
#include <boost/timer/timer.hpp>

#include <Eigen/Core>
#include <Eigen/Geometry>



// basic file operations
#include <fstream>

namespace stateObservation
{
  ///Dynamic sized scalar vector
  typedef Eigen::VectorXd Vector;

  ///1D Vector
  typedef Eigen::Matrix<double,1,1> Vector1;

  ///3D vector
  typedef Eigen::Vector3d Vector3;

  ///3D vector unaligned
  typedef Eigen::Matrix<double, 3, 1, Eigen::DontAlign> Vector3Unaligned;

  ///4D vector
  typedef Eigen::Vector4d Vector4;

  /// 6D vector
  typedef Eigen::Matrix<double,6,1> Vector6;

  ///Dynamic sized Matrix
  typedef Eigen::MatrixXd Matrix;

  ///3x3 Scalar Matrix
  typedef Eigen::Matrix3d Matrix3;

  ///3x3 Scalar Matrix Unaligned
  typedef Eigen::Matrix<double, 3, 3, Eigen::DontAlign> Matrix3Unaligned;

  ///4x4 Scalar Matrix
  typedef Eigen::Matrix4d Matrix4;

  ///6x6 Scalar Matrix
  typedef Eigen::Matrix<double, 6, 6> Matrix6;

  ///12x12 scalar Matrix
  typedef Eigen::Matrix<double, 12, 12> Matrix12;

  ///Quaternion
  typedef Eigen::Quaterniond Quaternion;

  ///Quaternion Unaligned
  typedef Eigen::Quaternion<double, Eigen::DontAlign> QuaternionUnaligned;

  ///Euler Axis/Angle representation of orientation
  typedef Eigen::AngleAxis<double> AngleAxis;

  typedef long int TimeIndex;
  typedef size_t TimeSize;


#ifndef NDEBUG
  static const bool isDebug=true;
#else
  static const bool isDebug=false;
#endif // NDEBUG

  template <typename T, const T& defaultValue=T(), bool debug=true>
  class DebugItem
  {
  public:
    DebugItem():b_(defaultValue) {}
    explicit DebugItem(const T& v):b_(v) {}
    inline T& operator=(T v)
    {
      return b_=v;
    }
    inline operator T()const
    {
      return b_;
    }
    inline T set(const T& v)
    {
      return b_=v;
    }
    T get()const
    {
      return b_;
    }
  private:
    T b_;

  };

  template <typename T, const T& defaultValue>
  class DebugItem<T,defaultValue,false>
  {
  public:
    DebugItem() {}
    explicit DebugItem(T v) {}
    inline T& operator=(T v)
    {
      return defaultValue;
    }
    inline operator T() const
    {
      return defaultValue;
    }
    inline T set(const T& v)
    {
      return defaultValue;
    }
    T get()const
    {
      return defaultValue;
    }
  private:
    T b_;
  };

  namespace detail
  {
    extern const bool defaultTrue;
    extern const char* defaultErrorMSG;
    extern const std::runtime_error defaultException;
    extern const std::exception* defaultExcepionAddr;

    void defaultSum(const  Vector& stateVector, const Vector& tangentVector, Vector& sum);
    void defaultDifference(const  Vector& stateVector1, const Vector& stateVector2, Vector& difference);

  }

  ///this is simply a structure allowing for automatically verifying that
  /// the item has been initialized or not. The chckitm_reset() function allows to
  /// set it back to "not initialized" state.
  /// -lazy means that the "set" value is true all the time if NDEBUG is defined
  /// -alwaysCheck means that the check is always performed and throws exception
  /// if it fails. Otherwise, the check is performed only for debug.
  /// warning, this has no effect if lazy is set to true
  /// - assertion means that an assertion will be introduced for the check.
  /// - eigenAlignedNew should be set to true if any alignment is required for the
  /// new operator (see eigen documentation)
  template <typename T, bool lazy=false, bool alwaysCheck = false,
                  bool assertion=true, bool eigenAlignedNew=false>
  class CheckedItem
  {
  public:

    /// The parameter initialize sets whether the isSet() parameter is initialized to false
    CheckedItem( bool initialize = true);
    explicit CheckedItem(const T&);
    CheckedItem( const CheckedItem &);
    virtual ~CheckedItem() {}

    inline CheckedItem& operator=(const CheckedItem & );
    inline T& operator=(const T&);
    inline operator T() const ;
    inline operator const T&() const ;

    inline T chckitm_getValue() const;

    inline T& operator()();
    inline const T& operator()() const;

    inline bool isSet() const;
    inline void reset();

    /// set the value of the initialization check boolean
    inline void set(bool value=true);

    void setAssertMessage(std::string s);
    void setExceptionPtr(std::exception* e);

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW_IF(eigenAlignedNew)
  protected:

    static const bool do_check_ = !lazy || isDebug;
    static const bool do_assert_ = do_check_ && assertion;
    static const bool do_exception_ = do_check_ && !assertion;

    typedef DebugItem<bool,detail::defaultTrue,do_check_> IsSet;
    typedef DebugItem<const char*,detail::defaultErrorMSG, do_assert_> AssertMsg;
    typedef DebugItem<const std::exception*,detail::defaultExcepionAddr,
                                                    do_exception_> ExceptionPtr;



    IsSet isSet_;
    AssertMsg assertMsg_;
    ExceptionPtr exceptionPtr_;

    bool chckitm_check_() const throw(std::exception);
    T v_;
  };

  typedef CheckedItem<Matrix3,false,false,true,true> CheckedMatrix3;
  typedef CheckedItem<Matrix6,false,false,true,true> CheckedMatrix6;
  typedef CheckedItem<Matrix12,false,false,true,true> CheckedMatrix12;
  typedef CheckedItem<Vector3,false,false,true,true> CheckedVector3;
  typedef CheckedItem<Vector6,false,false,true,true> CheckedVector6;
  typedef CheckedItem<Quaternion,false,false,true,true> CheckedQuaternion;


  /**
   * \class    IndexedMatrixT
   * \brief    This class describes a structure composed by a matrix
   *           of a given size and a time-index parameter. It can tell also if
   *           it initialized or not.
   *
   *
   */
  template <typename MatrixType=Matrix, bool lazy = false>
  class IndexedMatrixT:
    protected DebugItem<bool,detail::defaultTrue,!lazy || isDebug>
  {
    typedef DebugItem<bool,detail::defaultTrue,!lazy || isDebug> IsSet;
  public:
    ///Default constructor
    IndexedMatrixT();

    ///A constructor with a given matrix value and a time index
    IndexedMatrixT(const MatrixType& v, TimeIndex k);

    ///Set the value of the matrix and the time sample
    inline MatrixType set(const MatrixType& v,TimeIndex k);

    ///Switch the vector to "initialized" state
    inline void set(bool value=true);

    ///set the index of the matrix
    inline void setIndex(TimeIndex index);

    ///Get the matrix value
    inline const MatrixType & operator()() const;

    ///Get the matrix value
    inline MatrixType & operator()();

    ///Get the time index
    inline TimeIndex getTime() const;

    ///Says whether the matrix is initialized or not
    inline bool isSet() const;

    ///Switch off the initalization flag, the value is no longer accessible
    inline void reset();

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  protected:
    ///Checks whether the matrix is set or not (assert)
    ///does nothing in release mode
    inline bool check_() const;
    TimeIndex k_;
    MatrixType v_;
  };

  typedef IndexedMatrixT<Matrix> IndexedMatrix;
  typedef IndexedMatrixT<Vector> IndexedVector;
  typedef IndexedMatrixT<Vector3> IndexedVector3;
  typedef IndexedMatrixT<Matrix3> IndexedMatrix3;

  /**
   * \class    IndexedMatrixArray
   * \brief    This class describes a structure that enables to store array of matrices
   *           with time indexation.
   *
   */

  template <typename MatrixType=Matrix>
  class IndexedMatrixArrayT
  {
  public:
    ///Default constructor
    IndexedMatrixArrayT();



    typedef std::vector< MatrixType > Array;

    ///Sets the vector v at the time index k
    ///It checks the time index, the array must have contiguous indexes
    ///It can be used to push a value into the back of the array
    inline void setValue(const MatrixType& v, TimeIndex  k);

    ///Pushes back the matrix to the array, the new value will take the next time
    inline void pushBack(const MatrixType& v);

    ///removes the first (oldest) element of the array
    inline void popFront();

    ///gets the value with the given time index
    inline MatrixType operator[](TimeIndex index) const;

    ///gets the value with the given time index, non const version
    inline MatrixType  & operator[](TimeIndex index);

    ///gets the first value
    inline const MatrixType & front() const;

    ///gets the first value
    inline MatrixType& front();

    ///gets the last value
    inline const MatrixType & back() const;

    ///gets the last value
    inline MatrixType & back();

    ///removes all the elements with larger indexes than timeIndex
    void truncateAfter(TimeIndex timeIndex);

    ///removes all the elements with smaller indexes than timeIndex
    void truncateBefore(TimeIndex timeIndex);

    ///resizes the array
    inline void resize(TimeSize i, const MatrixType & m= MatrixType());

    ///Get the time index
    inline long int getLastIndex() const;

    ///Get the time index of the next value that will be pushed back
    /// Can be used in for loops
    inline TimeIndex getNextIndex() const;

    ///Set the time index of the last element
    inline TimeIndex setLastIndex(int index);

    ///Get the time index
    inline TimeIndex getFirstIndex() const;

    ///set the time index of the first element
    inline TimeIndex setFirstIndex(int index);

    inline TimeSize size() const;

    ///Resets the array to initial state
    ///the value is no longer accessible
    inline void reset();

    ///Clears the vector but keeps the last index
    inline void clear();

    ///converts the array into a standard vector
    Array getArray() const;

    ///checks whether the index is present in the array
    inline bool checkIndex(TimeIndex k) const;

    ///gets the array from a file
    ///the line starts with the time index and then the matrix is read
    ///row by row
    ///WARNING: this resets the array
    void readFromFile(const char * filename, size_t rows, size_t cols=1, bool withTimeStamp = true);
    void readFromFile(const std::string &  filename, size_t rows, size_t cols=1, bool withTimeStamp = true);

    ///gets the array from a file
    ///the line starts with the time index and then every line of the file
    /// is converted into a vector
    ///WARNING: this resets the array
    void readVectorsFromFile(const std::string &  filename, bool withTimeStamp = true );
    void readVectorsFromFile(const char * filename, bool withTimeStamp = true );

    ///write the array in a a file
    ///the line starts with the time index and then the matrix is described
    ///row by row
    /// When clear is set, the array is cleared but the time index is conserved
    /// When append is set to true, the output is appended to file
    void writeInFile(const char * filename, bool clear=false, bool append =false);

    ///write the array in a a file
    ///the line starts with the time index and then the matrix is described
    ///row by row
    /// When clear is set, the array is cleared but the time index is conserved
    /// When append is set to true, the output is appended to file
    void writeInFile(const std::string & filename, bool clear=false, bool append =false);

  protected:
    typedef std::deque< MatrixType > Deque;

    ///Asserts that the index is present in the array
    ///does nothing in release mode
    inline void check_(TimeIndex time) const;

    ///Asserts that the array is not empty
    ///does nothing in release mode
    inline void check_() const;

    ///Asserts that the given time can be pushed at the back of the vector
    ///does nothing in release mode
    inline void checkNext_(TimeIndex time) const;

    TimeIndex k_;

    Deque v_;

  };

  typedef IndexedMatrixArrayT<Matrix> IndexedMatrixArray;
  typedef IndexedMatrixArrayT<Vector> IndexedVectorArray;


  namespace cst
  {
    const double gravityConstant = 9.8;

    ///Gravity Vector along Z
    const Vector gravity= gravityConstant * Vector3::UnitZ();

    ///angles considered Zero
    const double epsilonAngle=1e-16;

  }

  typedef boost::timer::auto_cpu_timer auto_cpu_timer;
  typedef boost::timer::cpu_timer cpu_timer;
  typedef boost::timer::cpu_timer cpu_times;

  namespace tools
  {
    struct SimplestStopwatch
    {
      inline void start();

      ///provides the time since the start
      ///the value is in nanoseconds
      inline double stop();

      inline double diff(const timespec & start, const timespec & end);

      timespec time1, time2, time3;
    };

    std::string matrixToString(const Matrix& mat);

    std::string vectorToString(const Vector& v);

    Matrix stringToMatrix(const std::string& str, unsigned rows, unsigned cols);

    Vector stringToVector(const std::string& str, unsigned length);

    Vector stringToVector(const std::string& str);
  }


#include <state-observation/tools/definitions.hxx>
}

#endif //STATEOBSERVATIONDEFINITIONSHPP
