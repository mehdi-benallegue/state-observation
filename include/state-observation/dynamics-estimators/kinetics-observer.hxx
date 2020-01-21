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
inline unsigned KineticsObserver::gyroBiasIndex( unsigned numberOfIMU) const
{
  BOOST_ASSERT( numberOfIMU < maxImuNumber_ && \
                  "The requested IMU number is higher than the maximum");
  return angVelIndex() + sizeAngVel+ sizeGyroBias * numberOfIMU;
}
inline unsigned KineticsObserver::unmodeledWrenchIndex() const
{
  return gyroBiasIndex(0)+ sizeGyroBias*maxImuNumber_;
}
inline unsigned KineticsObserver::unmodeledForceIndex() const
{
  return unmodeledWrenchIndex();
}
inline unsigned KineticsObserver::contactsIndex() const
{
  return unmodeledWrenchIndex()+sizeWrench;
}
inline unsigned KineticsObserver::unmodeledTorqueIndex() const
{
  return unmodeledForceIndex()+sizeForce;
}
inline unsigned KineticsObserver::contactIndex(unsigned contactNbr) const
{
  BOOST_ASSERT(contactNbr < maxContacts_ && \
                  "The requested contact number is higher than the maximum");
  BOOST_ASSERT(contacts_[contactNbr].isSet && \
                  "The requested contact is not set yet, please add it before");
  return contacts_[contactNbr].stateIndex;
}
inline unsigned KineticsObserver::contactKineIndex(unsigned contactNbr) const
{
  return contactIndex(contactNbr);
}
inline unsigned KineticsObserver::contactPosIndex(unsigned contactNbr) const
{
  return contactKineIndex(contactNbr);
}
inline unsigned KineticsObserver::contactOriIndex(unsigned contactNbr) const
{
  return contactPosIndex(contactNbr)+sizePos;
}
inline unsigned KineticsObserver::contactWrenchIndex(unsigned contactNbr) const
{
  return contactKineIndex(contactNbr)+sizeContactKine;
}
inline unsigned KineticsObserver::contactForceIndex(unsigned contactNbr) const
{
  return contactWrenchIndex(contactNbr);
}
inline unsigned KineticsObserver::contactTorqueIndex(unsigned contactNbr) const
{
  return contactForceIndex(contactNbr)+sizeForce;
}

inline unsigned KineticsObserver::contactIndex(VectorContactConstIterator i) const
{  
  return i->stateIndex;
}
inline unsigned KineticsObserver::contactKineIndex(VectorContactConstIterator i) const
{
  return contactIndex(i);
}
inline unsigned KineticsObserver::contactPosIndex(VectorContactConstIterator i) const
{
  return contactKineIndex(i);
}
inline unsigned KineticsObserver::contactOriIndex(VectorContactConstIterator i) const
{
  return contactPosIndex(i)+sizePos;
}
inline unsigned KineticsObserver::contactWrenchIndex(VectorContactConstIterator i) const
{
  return contactKineIndex(i)+sizeContactKine;
}
inline unsigned KineticsObserver::contactForceIndex(VectorContactConstIterator i) const
{
  return contactWrenchIndex(i);
}
inline unsigned KineticsObserver::contactTorqueIndex(VectorContactConstIterator i) const
{
  return contactForceIndex(i)+sizeForce;
}


inline unsigned KineticsObserver::kineIndexTangent() const
{
  return 0;
}
inline unsigned KineticsObserver::posIndexTangent() const
{
  return 0;
}
inline unsigned KineticsObserver::oriIndexTangent() const
{
  return posIndexTangent()+sizePos;
}
inline unsigned KineticsObserver::linVelIndexTangent() const
{
  return oriIndexTangent()+sizeOriTangent;
}
inline unsigned KineticsObserver::angVelIndexTangent() const
{
  return linVelIndexTangent()+sizeLinVel;
}
inline unsigned KineticsObserver::gyroBiasIndexTangent( unsigned numberOfIMU) const
{
  BOOST_ASSERT( numberOfIMU < maxImuNumber_ && \
                  "The requested IMU number is higher than the maximum");
  return angVelIndexTangent() + sizeAngVel+ sizeGyroBias * numberOfIMU;
}
inline unsigned KineticsObserver::unmodeledWrenchIndexTangent() const
{
  return gyroBiasIndexTangent(0)+sizeGyroBias*maxImuNumber_;
}
inline unsigned KineticsObserver::unmodeledForceIndexTangent() const
{
  return unmodeledWrenchIndexTangent();
}
inline unsigned KineticsObserver::contactsIndexTangent() const
{
  return unmodeledWrenchIndexTangent()+sizeWrench;
}
inline unsigned KineticsObserver::unmodeledTorqueIndexTangent() const
{
  return unmodeledForceIndexTangent()+sizeForce;
}
inline unsigned KineticsObserver::contactIndexTangent(unsigned contactNbr) const
{
  BOOST_ASSERT(contactNbr < maxContacts_ && \
                  "The requested contact number is higher than the maximum");
  BOOST_ASSERT(contacts_[contactNbr].isSet && \
                  "The requested contact is not set yet, please add it before");
  return contacts_[contactNbr].stateIndexTangent;
}
inline unsigned KineticsObserver::contactKineIndexTangent(unsigned contactNbr) const
{
  return contactIndexTangent(contactNbr);
}
inline unsigned KineticsObserver::contactPosIndexTangent(unsigned contactNbr) const
{
  return contactKineIndexTangent(contactNbr);
}
inline unsigned KineticsObserver::contactOriIndexTangent(unsigned contactNbr) const
{
  return contactPosIndexTangent(contactNbr)+sizePos;
}
inline unsigned KineticsObserver::contactWrenchIndexTangent(unsigned contactNbr) const
{
  return contactKineIndexTangent(contactNbr)+sizeContactKineTangent;
}
inline unsigned KineticsObserver::contactForceIndexTangent(unsigned contactNbr) const
{
  return contactWrenchIndexTangent(contactNbr);
}
inline unsigned KineticsObserver::contactTorqueIndexTangent(unsigned contactNbr) const
{
  return contactForceIndexTangent(contactNbr)+sizeForce;
}

inline unsigned KineticsObserver::contactIndexTangent(VectorContactConstIterator i) const
{  
  return i->stateIndexTangent;
}
inline unsigned KineticsObserver::contactKineIndexTangent(VectorContactConstIterator i) const
{
  return contactIndexTangent(i);
}
inline unsigned KineticsObserver::contactPosIndexTangent(VectorContactConstIterator i) const
{
  return contactKineIndexTangent(i);
}
inline unsigned KineticsObserver::contactOriIndexTangent(VectorContactConstIterator i) const
{
  return contactPosIndexTangent(i)+sizePos;
}
inline unsigned KineticsObserver::contactWrenchIndexTangent(VectorContactConstIterator i) const
{
  return contactKineIndexTangent(i)+sizeContactKineTangent;
}
inline unsigned KineticsObserver::contactForceIndexTangent(VectorContactConstIterator i) const
{
  return contactWrenchIndexTangent(i);
}
inline unsigned KineticsObserver::contactTorqueIndexTangent(VectorContactConstIterator i) const
{
  return contactForceIndexTangent(i)+sizeForce;
}







