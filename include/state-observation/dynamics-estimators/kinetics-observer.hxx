inline unsigned KineticsObserver::kineIndex() const
{
  return 0;
}
inline unsigned KineticsObserver::posIndex() const
{
  return 0;
}
inline unsigned KineticsObserver::oriIndex() const
{
  return posIndex()+sizePos;
}
inline unsigned KineticsObserver::linVelIndex() const
{
  return oriIndex()+sizeOri;
}
inline unsigned KineticsObserver::angVelIndex() const
{
  return linVelIndex()+sizeLinVel;
}
inline unsigned KineticsObserver::gyroBiasIndex() const
{
  return angVelIndex() + sizeAngVel;
}
inline unsigned KineticsObserver::unmodeledWrenchIndex() const
{
  return gyroBiasIndex()+sizeGyroBias;
}
inline unsigned KineticsObserver::unmodeledForceIndex() const
{
  return unmodeledWrenchIndex();
}
inline unsigned KineticsObserver::unmodeledTorqueIndex() const
{
  return unmodeledForceIndex()+sizeForce;
}
inline unsigned KineticsObserver::contactKineIndex(int contactNbr) const
{
  return contacts_.find(contactNbr)->second.stateIndex;
}
inline unsigned KineticsObserver::contactPosIndex(int contactNbr) const
{
  return contactKineIndex(contactNbr);
}
inline unsigned KineticsObserver::contactOriIndex(int contactNbr) const
{
  return contactPosIndex(contactNbr)+sizePos;
}
inline unsigned KineticsObserver::contactWrenchIndex(int contactNbr) const
{
  return contactKineIndex(contactNbr)+sizeContactKine;
}
inline unsigned KineticsObserver::contactForceIndex(int contactNbr) const
{
  return contactWrenchIndex(contactNbr);
}
inline unsigned KineticsObserver::contactTorqueIndex(int contactNbr) const
{
  return contactForceIndex(contactNbr)+sizeForce;
}

template <int blockSize>
void KineticsObserver::setBlockStateCovariance(Matrix & covMat, const Matrix & covBlock, int blockIndex,int matrixSize)
{
  covMat.block<blockSize,blockSize>(blockIndex,blockIndex)=covBlock;
  covMat.block(blockIndex,0,blockSize,blockIndex).setZero();
  covMat.block(0,blockIndex,blockIndex,blockSize).setZero();
  covMat.block(blockIndex+blockSize,blockIndex,matrixSize-blockIndex-blockSize,blockSize).setZero();
  covMat.block(blockIndex,blockIndex+blockSize,blockSize,matrixSize-blockIndex-blockSize).setZero();
}







