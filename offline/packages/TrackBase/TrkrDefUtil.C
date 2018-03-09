#include <bitset>
#include "TrkrDefUtil.h"

void TrkrDefUtil::PrintBits(const TrkrDefs::hitsetkey key, std::ostream& os)
{
  os << "key: " << std::bitset<32>(key) << std::endl;
}

void TrkrDefUtil::PrintBits(const TrkrDefs::cluskey key, std::ostream& os)
{
  os << "key: " << std::bitset<64>(key) << std::endl;
}

uint8_t
TrkrDefUtil::GetTrkrId(const TrkrDefs::hitsetkey key)
{
  TrkrDefs::hitsetkey tmp = (key >> kBitShiftTrkrId);
  return tmp;
}

uint8_t
TrkrDefUtil::GetTrkrId(const TrkrDefs::cluskey key)
{
  TrkrDefs::hitsetkey tmp = (key >> kBitShiftClusId);
  return GetTrkrId(tmp);
}

uint8_t
TrkrDefUtil::GetLayer(const TrkrDefs::hitsetkey key)
{
  TrkrDefs::hitsetkey tmp = (key >> kBitShiftLayer);
  return tmp;
}

uint8_t
TrkrDefUtil::GetLayer(const TrkrDefs::cluskey key)
{
  TrkrDefs::hitsetkey tmp = (key >> kBitShiftClusId);
  return GetLayer(tmp);
}

uint32_t
TrkrDefUtil::GetClusIndex(const TrkrDefs::cluskey key)
{
  return key;
}

uint32_t
TrkrDefUtil::GetHitSetKeyFromClusKey(const TrkrDefs::cluskey key)
{
  return (key >> kBitShiftClusId);
}

TrkrDefs::hitsetkey
TrkrDefUtil::GetHitSetKeyLo(const TrkrDefs::TRKRID trkr_id)
{
  return GenHitSetKey(trkr_id, 0);
}

TrkrDefs::hitsetkey
TrkrDefUtil::GetHitSetKeyHi(const TrkrDefs::TRKRID trkr_id)
{
  return GenHitSetKey(static_cast<TrkrDefs::TRKRID>(trkr_id + 1), 0) - 1;
}

TrkrDefs::hitsetkey
TrkrDefUtil::GetHitSetKeyLo(const TrkrDefs::TRKRID trkr_id, const char lyr)
{
  return GenHitSetKey(trkr_id, lyr);
}

TrkrDefs::hitsetkey
TrkrDefUtil::GetHitSetKeyHi(const TrkrDefs::TRKRID trkr_id, const char lyr)
{
  return GenHitSetKey(trkr_id, lyr + 1) - 1;
}

TrkrDefs::cluskey
TrkrDefUtil::GetClusKeyLo(const TrkrDefs::TRKRID trkr_id)
{
  TrkrDefs::cluskey tmp = GenHitSetKey(trkr_id, 0);
  return (tmp << kBitShiftClusId);
}

TrkrDefs::cluskey
TrkrDefUtil::GetClusKeyHi(const TrkrDefs::TRKRID trkr_id)
{
  TrkrDefs::cluskey tmp = GenHitSetKey(static_cast<TrkrDefs::TRKRID>(trkr_id + 1), 0);
  return (tmp << kBitShiftClusId) - 1;
}

TrkrDefs::cluskey
TrkrDefUtil::GetClusKeyLo(const TrkrDefs::TRKRID trkr_id, const char lyr)
{
  TrkrDefs::cluskey tmp = GenHitSetKey(trkr_id, lyr);
  return (tmp << kBitShiftClusId);
}

TrkrDefs::cluskey
TrkrDefUtil::GetClusKeyHi(const TrkrDefs::TRKRID trkr_id, const char lyr)
{
  TrkrDefs::cluskey tmp = GenHitSetKey(trkr_id, lyr + 1);
  return (tmp << kBitShiftClusId) - 1;
}

TrkrDefs::hitsetkey
TrkrDefUtil::GenHitSetKey(const TrkrDefs::TRKRID trkr_id, const char lyr)
{
  TrkrDefs::hitsetkey tmp = trkr_id;
  TrkrDefs::hitsetkey key = tmp << kBitShiftTrkrId;  // detector id
  tmp = lyr;
  key |= (tmp << kBitShiftLayer);  // layer
  return key;
}
