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
fillClusters(detId_t idet, Clusterizer& clusterizer, const std::vector<FEDChannel>& channels, std::ofstream& digidata_out)
{
  static bool first = true;
  std::vector<SiStripCluster> out;
  std::vector<SiStripDigi> detDigis;
  std::vector<uint16_t> seedStrips;
  std::vector<uint16_t> seedStrips_noduplicate;

  // start clusterizing for idet
  //auto const & det = clusterizer.stripByStripBegin(idet);
  //state.reset(det);

  // unpack RAW data to SiStripDigi
  for (auto const& chan : channels) {
    //        std::cout << "Processing channel for detid " << idet << " fed " << chan.fedId() << " channel " << (int) chan.fedCh() << " len:off " << chan.length() << ":" << chan.offset() << " ipair " << chan.iPair() << std::endl;
    //auto perStripAdder = StripByStripAdder(clusterizer, state, out);
    //unpackZS(chan, chan.iPair()*256, perStripAdder, idet);
    unpackZS2(chan, chan.iPair()*256, detDigis, idet);
  }

  // start clusterizing for idet
  auto& det = clusterizer.getDet(idet);

  // initialize ADC for strips
  for (auto const& digi : detDigis) {
    det.setADC(digi.strip(), digi.adc());
  }

  for (auto const& chan : channels) {
    for(auto const& digi : detDigis) {
      fedId_t fedId=chan.fedId();
      fedCh_t fedCh=chan.fedCh();
      uint16_t strip=digi.strip();
      uint16_t adc=digi.adc();
      float noise=det.noise(strip);
      float gain=det.gain(strip);
      bool bad=det.bad(strip);
      digidata_out.write((char*)&idet, sizeof(detId_t));
      digidata_out.write((char*)&fedId, sizeof(fedId_t));
      digidata_out.write((char*)&fedCh, sizeof(fedCh_t));
      digidata_out.write((char*)&strip, sizeof(uint16_t));
      digidata_out.write((char*)&adc, sizeof(uint16_t));
      digidata_out.write((char*)&noise, sizeof(float));
      digidata_out.write((char*)&gain, sizeof(float));
      digidata_out.write((char*)&bad, sizeof(bool));
    }
  }

  // create seedStrips
  for (auto const& digi : detDigis) {
    Clusterizer::State state(det);
    if (clusterizer.seedStrip(state, digi.strip())) {
      seedStrips.push_back(digi.strip());
      //std::cout<<"detID="<<idet<<"seed strip="<<digi.strip()<<std::endl;
    }
  }

  // create non-duplicated seedStrips
  if (!seedStrips.empty()) seedStrips_noduplicate.push_back(seedStrips[0]);
  for (int i=1; i<seedStrips.size(); i++) {
    if (seedStrips[i] - seedStrips[i-1]!=1)
      seedStrips_noduplicate.push_back(seedStrips[i]);
  }

  // process each seedStrip digi
  //for (auto const& digi : detDigis) {
  for (auto ss : seedStrips_noduplicate) {
    //std::cout<<"detID="<<idet<<"seed strip="<<ss<<std::endl;
    Clusterizer::State state(det);
    state.reset(det);
    if (clusterizer.seedStrip(state, ss)) {
      SiStripCluster cluster;
      clusterizer.findCluster(state, ss, cluster);
      if (cluster.amplitudes().size() > 0)
	out.push_back(cluster);
    }
  }

  //if (first) {
    first = false;
    //std::cout << "Printing clusters for detid " << idet << std::endl;
    for (const auto& cluster : out) {
      std::cout <<" det id "<<idet<<" strip " << cluster.firstStrip() << ": ";
      for (const auto& ampl : cluster.amplitudes()) {
        std::cout << (int) ampl << " ";
      }
      std::cout << std::endl;
    }
    //}

  return out;
}

int main()
{
  //std::cout << "start initializing clusterizer" << std::endl;
  Clusterizer clusterizer;
  //  Clusterizer::State state;

  std::ofstream digidata_out("digidata.bin", std::ofstream::out | std::ios::binary);
  std::ifstream digidata_in("digidata.bin", std::ofstream::in | std::ios::binary);

  FEDSet feds(fillFeds());
  for (auto idet : clusterizer.allDetIds()) {
    if (feds.find(idet) != feds.end()) {
      auto out = fillClusters(idet, clusterizer, feds[idet], digidata_out);
    }
  }

  // test digidata output]
  detId_t detid;
  while (digidata_in.read((char*)&detid, sizeof(detid)).gcount() == sizeof(detid)) {
    fedId_t fedId;
    fedCh_t fedCh;
    uint16_t strip;
    uint16_t adc;
    float noise;
    float gain;
    bool bad;
    digidata_in.read((char*)&fedId, sizeof(fedId_t));
    digidata_in.read((char*)&fedCh, sizeof(fedCh_t));
    digidata_in.read((char*)&strip, sizeof(uint16_t));
    digidata_in.read((char*)&adc, sizeof(uint16_t));
    digidata_in.read((char*)&noise, sizeof(float));
    digidata_in.read((char*)&gain, sizeof(float));
    digidata_in.read((char*)&bad, sizeof(bool));
    //    std::cout<<" detid "<< detid <<" fedId "<<fedId<<" fedCh "<<(int)fedCh<<" strip "<<strip<<" adc "<<adc<<" noise "<<noise<<" gain "<<gain<<" bad "<<bad<<std::endl;
  }

  return 0;

}
