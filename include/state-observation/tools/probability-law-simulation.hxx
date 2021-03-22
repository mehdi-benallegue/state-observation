template<typename ReturnType, typename StdType, typename BiasType>
typename MatrixType<ReturnType>::type ProbabilityLawSimulation::getGaussianVector(StdType std,
                                                                                  BiasType bias,
                                                                                  Index rows,
                                                                                  Index cols)
{
  static_assert(isEigen<StdType>::value && isEigen<BiasType>::value,
                "Standard deviation and bias need to be eigen matrices");

  std::normal_distribution<double> g(0, 1);

  ReturnType ret;

  if(ReturnType::RowsAtCompileTime == -1 || ReturnType::ColsAtCompileTime == -1) // the matrix is dynamic size
  {
    ret = Matrix::Zero(rows, cols);
  }
  else
  {
    ret.setZero();
  }

  for(Index i = 0; i < rows; ++i)
  {
    for(Index j = 0; j < cols; ++j)
    {
      ret(i, j) = g(gen_);
    }
  }
  ret = std * ret + bias;

  return ret;
}