#include <fstream>
#include <state-observation/tools/logger.hpp>


namespace stateObservation
{
  namespace tools
  {
    Logger::Logger()
    {
      scalar_.resize(1,1);
    }

    Logger::~Logger()
    {
      //dtor
    }

     void Logger::setPath(const std::string& path)
    {
      path_=path;

    }

    void Logger::save(bool clear, bool append)
    {
      for (std::map<const void *, log_s >::iterator i=logs_.begin();i!=logs_.end();++i)
      {
        if (i->second.filename!=std::string(""))
        {
          i->second.array.writeInFile((path_+std::string("/")+i->second.filename).c_str(),
                                      clear,append);
        }
      }
    }

    void Logger::clearTracking()
    {
      logs_.clear();
    }

    void Logger::clearLogs()
    {
      for (std::map<const void *, log_s >::iterator i=logs_.begin();i!=logs_.end();++i)
      {
        i->second.array.clear();
      }
    }

    const IndexedMatrixArray & Logger::getRecord(const void* address) const
    {
      Tmap::const_iterator i= logs_.find(address);
      if (i==logs_.end())
      {
        BOOST_ASSERT(false && "The logger cannot find the data, please use record function");
        throw std::invalid_argument
                ("The logger cannot find the data, please use record function");
      }
      else
        return i->second.array;
    }

    IndexedMatrixArray & Logger::getRecord(const void* address)
    {
      Tmap::iterator i= logs_.find(address);
      if (i==logs_.end())
      {
        BOOST_ASSERT(false && "The logger cannot find the data, please use record function");
        throw std::invalid_argument
                ("The logger cannot find the data, please use record function");
      }
      else
        return i->second.array;
    }

    void Logger::push()
    {
      for (std::map<const void *, log_s >::iterator i=logs_.begin();i!=logs_.end();++i)
      {
        update_(i);
      }
    }

    void Logger::update_(const Tmap::iterator & i)
    {

      ///If it is a Matrix or a Vector type
      if (*i->second.type == typeid(Matrix))
      {
        i->second.array.pushBack(*static_cast<const Matrix *>(i->first));
        return;
      }
      if (*i->second.type == typeid(Vector))
      {
        i->second.array.pushBack(*static_cast<const Vector *>(i->first));
        return;
      }
      if (*i->second.type == typeid(Matrix3))
      {
        i->second.array.pushBack(*static_cast<const Matrix3 *>(i->first));
        return;
      }
      if (*i->second.type == typeid(Vector3))
      {
        i->second.array.pushBack(*static_cast<const Vector3 *>(i->first));
        return;
      }

      ///if it is a scalar type
      if (*i->second.type == typeid(int))
      {
        scalar_(0,0)=double(*static_cast<const int *>(i->first));
      }
      else if (*i->second.type == typeid(unsigned))
      {
        scalar_(0,0)=double(*static_cast<const unsigned *>(i->first));
      }
      else if (*i->second.type == typeid(long))
      {
        scalar_(0,0)=double(*static_cast<const long *>(i->first));
      }
      else if (*i->second.type == typeid(size_t))
      {
        scalar_(0,0)=double(*static_cast<const size_t *>(i->first));
      }
      else if (*i->second.type == typeid(double))
      {
        scalar_(0,0)=*static_cast<const double *>(i->first);
      }
      else if (*i->second.type == typeid(float))
      {
        scalar_(0,0)=double(*static_cast<const float *>(i->first));
      }

      i->second.array.pushBack(scalar_);
    }

  }
}
