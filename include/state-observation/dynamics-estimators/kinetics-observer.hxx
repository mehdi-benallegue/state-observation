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

inline unsigned KineticsObserver::contactIndex(MapContactIterator i) const
{  
  return i->second.stateIndex;
}
inline unsigned KineticsObserver::contactKineIndex(MapContactIterator i) const
{
  return contactIndex(i);
}
inline unsigned KineticsObserver::contactPosIndex(MapContactIterator i) const
{
  return contactKineIndex(i);
}
inline unsigned KineticsObserver::contactOriIndex(MapContactIterator i) const
{
  return contactPosIndex(i)+sizePos;
}
inline unsigned KineticsObserver::contactWrenchIndex(MapContactIterator i) const
{
  return contactKineIndex(i)+sizeContactKine;
}
inline unsigned KineticsObserver::contactForceIndex(MapContactIterator i) const
{
  return contactWrenchIndex(i);
}
inline unsigned KineticsObserver::contactTorqueIndex(MapContactIterator i) const
{
  return contactForceIndex(i)+sizeForce;
}








