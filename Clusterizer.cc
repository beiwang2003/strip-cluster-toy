#include <numeric>
#include <cmath>
#include <iostream>

#include "Clusterizer.h"

template<typename Iter>
inline float chargePerCM(Iter a, Iter b) {
  return float(std::accumulate(a,b,int(0)))/0.047f;
}


Clusterizer::Det::Det(detId_t detid, std::ifstream& file)
  : id_(detid)
{
  fedId_t fid;

  std::vector<DetStrip> dets;

  while (file.read((char*)&fid, sizeof(fid)).gcount() == sizeof(fid) && fid != UINT16_MAX) {
    fedIDs_.push_back(fid);
    size_t count;
    file.read((char*)&count, sizeof(count));
    //std::cout << "Fed " << fid << " det " << detid << std::endl;
    if (count > 0) {
      std::vector<DetStrip> tmpDets(count);
      file.read((char*)tmpDets.data(), count*sizeof(DetStrip));
      dets.insert(dets.end(), tmpDets.begin(), tmpDets.end());
    }
  }
  if (dets.size() > 0) {
    offset_ = dets.front().strip_;
    auto last = dets.back().strip_;
    strips_.resize(last - offset_ + 1);
    ADCs_.resize(last - offset_ + 1);
    valid_.resize(last -offset_ + 1);
    for (int i=0; i<ADCs_.size(); i++) {ADCs_[i] = 0;}
    initBad();
    //std::cout <<"strips size"<<strips_.size();
    for (const auto& det : dets) {
      auto ind = det.strip_-offset_;
      //std::cout << "strip #" << det.strip_ << " " << "ind" << ind << std::endl;
      strips_[det.strip_-offset_] = det;
      valid_[det.strip_-offset_] = true;
    }
    //std::cout << std::endl;
  }
  //  std::cout << "det " <<detid << "DetStrip size" << strips_.size() << "offset " << offset_ << std::endl;
}

Clusterizer::Clusterizer()
{
  std::ifstream detfile("stripdets.bin", std::ios::in | std::ios::binary);

  detId_t detid;

  while (detfile.read((char*)&detid, sizeof(detid)).gcount() == sizeof(detid)) {
    detIds_.push_back(detid);
    Det det(detid, detfile);
    dets_[det.id()] = std::move(det);
  }
}

Clusterizer::Det& Clusterizer::getDet(const uint32_t id) {
  auto it = dets_.find(id);
  if (it != dets_.end()) {
    return it->second;
  }
  std::cout << "no detector id found, exit"  << std::endl;
  exit (2);
}

//Clusterizer::Det Clusterizer::findDetId(const uint32_t id) const
//{
//  const auto it = dets_.find(id);
//  if (it != dets_.end()) {
//    return it->second;
//  }
//  return Det();
//}

inline bool
Clusterizer::candidateEndedLeft(State const & state, const uint16_t& testStrip) const {
  uint16_t holes = state.lastStripLeft - testStrip - 1;
  return ( testStrip < state.det().getOffset() ||
           (( (!state.ADCs.empty())  &&                    // a candidate exists, and
             (holes > MaxSequentialHoles )       // too many holes if not all are bad strips, and
           ) &&
           ( holes > MaxSequentialBad ||       // (too many bad strips anyway, or
             !state.det().allBadBetween( testStrip, state.lastStripLeft ) // not all holes are bad strips)
	     ))
           );
}

inline bool
Clusterizer::candidateEndedRight(State const & state, const uint16_t& testStrip) const {
  uint16_t holes = testStrip - state.lastStripRight - 1;
  //  std::cout<<"holes"<<holes<<std::endl;
  return ( testStrip > state.det().getLastStripID() ||
	  ( (!state.ADCs.empty())  &&                    // a candidate exists, and
	    (holes > MaxSequentialHoles )       // too many holes if not all are bad strips, and
	  ) &&
	  ( holes > MaxSequentialBad ||       // (too many bad strips anyway, or
	    !state.det().allBadBetween( state.lastStripRight, testStrip ) // not all holes are bad strips)
	  )
	 );
}

inline void
Clusterizer::addToCandidateLeft(State & state, uint16_t strip) const {
  float Noise = state.det().noise( strip );
  uint8_t adc = state.det().getADC( strip );
  if(  !state.det().valid( strip ) || Noise ==0 || adc < static_cast<uint8_t>( Noise * ChannelThreshold) || state.det().bad(strip) )
    return;

  while( --state.lastStripLeft < strip ) state.ADCs.push_back(0); // pad holes

  state.ADCs.insert( state.ADCs.begin(), adc);
  state.noiseSquared += Noise*Noise;
}

inline void
Clusterizer::addToCandidateRight(State & state, uint16_t strip) const {
  float Noise = state.det().noise( strip );
  uint8_t adc = state.det().getADC( strip);
  if( !state.det().valid( strip ) || Noise ==0 || adc < static_cast<uint8_t>( Noise * ChannelThreshold) || state.det().bad(strip) )
    return;

  //std::cout<<"strip"<<strip<<"noise"<<Noise<<"adc"<<(int)adc<<"lastStripRight"<<state.lastStripRight<<std::endl;
  while( ++state.lastStripRight < strip ) state.ADCs.push_back(0); // pad holes

  state.ADCs.push_back( adc );
  state.noiseSquared += Noise*Noise;
}

template <class T> inline void
Clusterizer::endCandidate(State & state, T& out) const {
  if(candidateAccepted(state)) {
    applyGains(state);
    appendBadNeighbors(state);
    if(chargePerCM(state.ADCs.begin(), state.ADCs.end()) > minGoodCharge)
      return out.initialize(firstStrip(state), state.ADCs.begin(), state.ADCs.end());
  }
  clearCandidate(state);
}

inline bool
Clusterizer::candidateAccepted(State const & state) const {
  return ( state.noiseSquared * ClusterThresholdSquared
	   <=  std::pow( float(std::accumulate(state.ADCs.begin(),state.ADCs.end(), int(0))), 2.f));
}

inline void
Clusterizer::applyGains(State & state) const {
  uint16_t strip = firstStrip(state);
  for( auto &  adc :  state.ADCs) {
    auto charge = int( float(adc)/state.det().gain(strip++) + 0.5f ); //adding 0.5 turns truncation into rounding
    if(adc < 254) adc = ( charge > 1022 ? 255 :
			  ( charge >  253 ? 254 : charge ));
  }
}

inline void
Clusterizer::appendBadNeighbors(State & state) const {
  uint8_t max = MaxAdjacentBad;
  while(0 < max--) {
    if( state.det().bad( firstStrip(state)-1) ) { std::cout<<"bad strip in the left"<<std::endl; state.ADCs.insert( state.ADCs.begin(), 0);  }
    if( state.det().bad(  state.lastStripRight + 1) ) { std::cout<<"bad strip in the right"<<std::endl; state.ADCs.push_back(0); state.lastStripRight++; }
  }
}

//Clusterizer::Det
//Clusterizer::stripByStripBegin(uint32_t id) const {
//return findDetId( id );
//}

//void
//Clusterizer::stripByStripAdd(State & state, uint16_t strip, uint8_t adc, std::vector<SiStripCluster>& out) const {

//  if(candidateEnded_right(state,strip)) endCandidate(state,out);
//  addToCandidate(state, SiStripDigi(strip,adc));
//}

void
Clusterizer::findLeftBoundary(State & state, uint16_t strip) const {
  strip--;
  //std::cout <<"left strip in"<<strip;
  while(!candidateEndedLeft(state,strip)) {
    addToCandidateLeft(state, strip);
    strip--;
  }
  //std::cout <<"left strip out"<<strip << std::endl;
}

void
Clusterizer::findRightBoundary(State & state, uint16_t strip) const{
  strip++;
  //std::cout <<"right strip in"<<strip;
  while(!candidateEndedRight(state,strip)) {
    addToCandidateRight(state, strip);
    strip++;
    //std::cout <<"right strip mid"<<strip;
  }
  //std::cout <<"right strip out"<<strip << std::endl;
}

bool
Clusterizer::seedStrip(State& state, uint16_t strip) const {
  float Noise = state.det().noise( strip );
  uint8_t adc = state.det().getADC( strip );
  return  adc >= static_cast<uint8_t>( Noise * SeedThreshold);
}

void
Clusterizer::findCluster(State & state, uint16_t strip, SiStripCluster& cluster) const {
  state.lastStripLeft = strip;
  state.lastStripRight = strip;
  float Noise = state.det().noise( strip );
  uint8_t adc = state.det().getADC(strip);
  state.ADCs.push_back( adc );
  state.noiseSquared += Noise*Noise;
  findLeftBoundary(state, strip);
  findRightBoundary(state, strip);
  endCandidate(state, cluster);
}

//void
//Clusterizer::stripByStripEnd(State & state, std::vector<SiStripCluster>& out) const {
//  endCandidate(state, out);
//}
