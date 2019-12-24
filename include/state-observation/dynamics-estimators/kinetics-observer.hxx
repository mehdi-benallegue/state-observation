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
inline unsigned KineticsObserver::contactsIndex() const
{
  return unmodeledWrenchIndex()+sizeWrench;
}
inline unsigned KineticsObserver::unmodeledTorqueIndex() const
{
  return unmodeledForceIndex()+sizeForce;
}
inline unsigned KineticsObserver::contactIndex(int contactNbr) const
{
    BOOST_ASSERT(contacts_.find(contactNbr)!=contacts_.end() && \
                  "The requested contact is not set yet, please add it before");
  return contacts_.find(contactNbr)->second.stateIndex;
}
inline unsigned KineticsObserver::contactKineIndex(int contactNbr) const
{
  return contactIndex(contactNbr);
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

inline unsigned KineticsObserver::contactIndex(MapContactConstIterator i) const
{  
  return i->second.stateIndex;
}
inline unsigned KineticsObserver::contactKineIndex(MapContactConstIterator i) const
{
  return contactIndex(i);
}
inline unsigned KineticsObserver::contactPosIndex(MapContactConstIterator i) const
{
  return contactKineIndex(i);
}
inline unsigned KineticsObserver::contactOriIndex(MapContactConstIterator i) const
{
  return contactPosIndex(i)+sizePos;
}
inline unsigned KineticsObserver::contactWrenchIndex(MapContactConstIterator i) const
{
  return contactKineIndex(i)+sizeContactKine;
}
inline unsigned KineticsObserver::contactForceIndex(MapContactConstIterator i) const
{
  return contactWrenchIndex(i);
}
inline unsigned KineticsObserver::contactTorqueIndex(MapContactConstIterator i) const
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
inline unsigned KineticsObserver::gyroBiasIndexTangent() const
{
  return angVelIndexTangent() + sizeAngVel;
}
inline unsigned KineticsObserver::unmodeledWrenchIndexTangent() const
{
  return gyroBiasIndexTangent()+sizeGyroBias;
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
inline unsigned KineticsObserver::contactIndexTangent(int contactNbr) const
{
    BOOST_ASSERT(contacts_.find(contactNbr)!=contacts_.end() && \
                  "The requested contact is not set yet, please add it before");
  return contacts_.find(contactNbr)->second.stateIndexTangent;
}
inline unsigned KineticsObserver::contactKineIndexTangent(int contactNbr) const
{
  return contactIndexTangent(contactNbr);
}
inline unsigned KineticsObserver::contactPosIndexTangent(int contactNbr) const
{
  return contactKineIndexTangent(contactNbr);
}
inline unsigned KineticsObserver::contactOriIndexTangent(int contactNbr) const
{
  return contactPosIndexTangent(contactNbr)+sizePos;
}
inline unsigned KineticsObserver::contactWrenchIndexTangent(int contactNbr) const
{
  return contactKineIndexTangent(contactNbr)+sizeContactKineTangent;
}
inline unsigned KineticsObserver::contactForceIndexTangent(int contactNbr) const
{
  return contactWrenchIndexTangent(contactNbr);
}
inline unsigned KineticsObserver::contactTorqueIndexTangent(int contactNbr) const
{
  return contactForceIndexTangent(contactNbr)+sizeForce;
}

inline unsigned KineticsObserver::contactIndexTangent(MapContactConstIterator i) const
{  
  return i->second.stateIndexTangent;
}
inline unsigned KineticsObserver::contactKineIndexTangent(MapContactConstIterator i) const
{
  return contactIndexTangent(i);
}
inline unsigned KineticsObserver::contactPosIndexTangent(MapContactConstIterator i) const
{
  return contactKineIndexTangent(i);
}
inline unsigned KineticsObserver::contactOriIndexTangent(MapContactConstIterator i) const
{
  return contactPosIndexTangent(i)+sizePos;
}
inline unsigned KineticsObserver::contactWrenchIndexTangent(MapContactConstIterator i) const
{
  return contactKineIndexTangent(i)+sizeContactKineTangent;
}
inline unsigned KineticsObserver::contactForceIndexTangent(MapContactConstIterator i) const
{
  return contactWrenchIndexTangent(i);
}
inline unsigned KineticsObserver::contactTorqueIndexTangent(MapContactConstIterator i) const
{
  return contactForceIndexTangent(i)+sizeForce;
}







