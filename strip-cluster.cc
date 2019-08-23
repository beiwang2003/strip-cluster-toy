#include <fstream>
#include <iostream>
#include <algorithm>

#include "Clusterizer.h"
#include "FEDChannel.h"
#include "FEDZSChannelUnpacker.h"

/*
class StripByStripAdder {
public:
  typedef std::output_iterator_tag iterator_category;
  typedef void value_type;
  typedef void difference_type;
  typedef void pointer;
  typedef void reference;

  StripByStripAdder(Clusterizer& clusterizer,
                    Clusterizer::State& state,
                    std::vector<SiStripCluster>& record)
    : clusterizer_(clusterizer), state_(state), record_(record) {}

  StripByStripAdder& operator= ( SiStripDigi digi )
  {
    clusterizer_.stripByStripAdd(state_, digi.strip(), digi.adc(), record_);
    return *this;
  }

  StripByStripAdder& operator*  ()    { return *this; }
  StripByStripAdder& operator++ ()    { return *this; }
  StripByStripAdder& operator++ (int) { return *this; }
private:
  Clusterizer& clusterizer_;
  Clusterizer::State& state_;
  std::vector<SiStripCluster>& record_;
};

template<typename OUT>
OUT unpackZS(const FEDChannel& chan, uint16_t stripOffset, OUT out, detId_t idet)
{
  auto unpacker = FEDZSChannelUnpacker::zeroSuppressedModeUnpacker(chan);
  while (unpacker.hasData()) {
    auto digi = SiStripDigi(stripOffset+unpacker.sampleNumber(), unpacker.adc());
    //    std::cout << "unpackZS det " << idet << " digi " << digi.strip() << " sample " << (unsigned int) unpacker.sampleNumber() << " adc " << (unsigned int) unpacker.adc() << std::endl;
    if (digi.strip() != 0) {
      *out++ = digi;
    }
    unpacker++;
  }
  return out;
}
*/

template<typename OUT>
void unpackZS2(const FEDChannel& chan, uint16_t stripOffset, OUT& out, detId_t idet)
{
  auto unpacker = FEDZSChannelUnpacker::zeroSuppressedModeUnpacker(chan);
  while (unpacker.hasData()) {
    auto digi = SiStripDigi(stripOffset+unpacker.sampleNumber(), unpacker.adc());
    //    std::cout << "unpackZS det " << idet << " digi " << digi.strip() << " sample " << (unsigned int) unpacker.sampleNumber() << " adc " << (unsigned int) unpacker.adc() << std::endl;
    if (digi.strip() != 0) {
      out.push_back(digi);
    }
    unpacker++;
  }
}

FEDSet fillFeds()
{
  std::ifstream fedfile("stripdata.bin", std::ios::in | std::ios::binary);

  FEDSet feds;
  detId_t detid;  

  while (fedfile.read((char*)&detid, sizeof(detid)).gcount() == sizeof(detid)) {
    FEDChannel fed(fedfile);
    //    std::cout << "Det " << detid << " fed " << fed.fedId() << " channel " << (int) fed.fedCh() << " length " << fed.length() << " offset " << fed.offset() << " ipair " << fed.iPair() << std::endl;
    feds[detid].push_back(std::move(fed));
  }
  return feds;
}

std::vector<SiStripCluster>
fillClusters(detId_t idet, Clusterizer& clusterizer, Clusterizer::State& state, const std::vector<FEDChannel>& channels)
{
  static bool first = true;
  std::vector<SiStripCluster> out;
  std::vector<SiStripDigi> input; 

  // start clusterizing for idet
  auto const & det = clusterizer.stripByStripBegin(idet);
  state.reset(det);

  // unpack RAW data to SiStripDigi
  for (auto const& chan : channels) {
    //        std::cout << "Processing channel for detid " << idet << " fed " << chan.fedId() << " channel " << (int) chan.fedCh() << " len:off " << chan.length() << ":" << chan.offset() << " ipair " << chan.iPair() << std::endl;
    //auto perStripAdder = StripByStripAdder(clusterizer, state, out);
    //unpackZS(chan, chan.iPair()*256, perStripAdder, idet);
    unpackZS2(chan, chan.iPair()*256, input, idet);
  }

  // process each digi 
  for (auto const& digi : input) {
    clusterizer.stripByStripAdd(state, digi.strip(), digi.adc(), out);
  }

  // end clusterizing for idet
  clusterizer.stripByStripEnd(state, out);

  if (first) {
    first = false;
    std::cout << "Printing clusters for detid " << idet << std::endl;
    for (const auto& cluster : out) {
      std::cout << "Cluster " << cluster.firstStrip() << ": ";
      for (const auto& ampl : cluster.amplitudes()) {
        std::cout << (int) ampl << " ";
      }
      std::cout << std::endl;
    }
  }

  return out;
}

int main()
{
  //std::cout << "start initializing clusterizer" << std::endl;
  Clusterizer clusterizer;
  Clusterizer::State state;

  //std::cout << "start fillFeds" <<std::endl;
  FEDSet feds(fillFeds());
  for (auto idet : clusterizer.allDetIds()) {
    if (feds.find(idet) != feds.end()) {
      auto out = fillClusters(idet, clusterizer, state, feds[idet]);
    }
  }
}
