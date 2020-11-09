#include <cstdio>
#include <vector>
#include <complex>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cstring>
#include <assert.h>
#include "simplelin.hh"
#include "ifloat.hh"
template <typename T> using complex = Complex<T>;
#include "p0.hh"

typedef SimpleFloat<uint32_t, uint64_t, 32, short> sfloat;

void usage() {
  std::cerr << "Usage: psnd [-e | -d]" << std::endl;
  exit(0);
}

inline int countMSB(const int16_t& x) {
  int res(0);
  for(int i = 1; i < 16; i ++)
    if(x & (1 << i)) res = i;
  return res;
}

SimpleVector<uint8_t> blockout(const SimpleVector<int16_t>& buf) {
  assert(buf.size() == 32);
  std::vector<std::pair<std::pair<int16_t, int16_t>, int> > s0;
  s0.reserve(buf.size());
  for(int i = 0; i < buf.size(); i ++)
    s0.emplace_back(std::make_pair(std::make_pair(abs(buf[i]), buf[i]), i));
  std::sort(s0.begin(), s0.end());
  DUInt<uint64_t, 64> permraw;
  permraw ^= permraw;
  std::vector<uint8_t> outbase;
  uint32_t punch;
  punch ^= punch;
  int  bidx(0);
  bool flag(false);
  for(int i = 0; i < buf.size(); i ++) {
    const auto& tb(s0[buf.size() - i - 1]);
    int cnt(0);
    for(int j = 0; j < tb.second; j ++)
      if(! (punch & (j ? (1 << j) : 1))) cnt ++;
    punch |= tb.second ? (1 << tb.second) : 1;
    assert(0 <= cnt && cnt < buf.size() - i);
    permraw += DUInt<uint64_t, 64>(cnt);
    permraw *= DUInt<uint64_t, 64>(buf.size() - i);
    if(flag) continue;
    if(! i) {
      outbase.emplace_back(tb.first.second & 0xff);
      outbase.emplace_back((tb.first.second >> 8) & 0xff);
    } else {
      auto bits(min(countMSB(s0[buf.size() - i].first.first) + 2, 16));
      while(0 < bits) {
        const auto residue(8 - bidx);
        if(! bidx) outbase.emplace_back(0);
        outbase[outbase.size() - 1] |=
          (bits > residue ? tb.first.second >> (bits - residue)
           : (bits == residue ? tb.first.second
                              : tb.first.second << (residue - bits))) &
          ((1 << residue) - 1);
        bidx += min(bits, residue);
        bidx %= 8;
        bits -= residue;
      }
    }
    if(! tb.first.first)
      flag = true;
  }
  assert(!((punch + 1) & 0xffffffff));
  static const auto szperm(sizeof(DUInt<uint64_t, 64>) / sizeof(uint8_t));
  SimpleVector<uint8_t> out(outbase.size() + szperm);
  for(int i = 0; i < szperm; i ++)
    out[i] = uint8_t(int(i ? (permraw >> (i * 8)) : permraw) & 0xff);
  for(int i = 0; i < outbase.size(); i ++)
    out[i + szperm] = outbase[i];
  return out;
}

SimpleVector<int16_t> blockin(std::istream& in) {
  SimpleVector<int16_t> out(32);
  // decode permraw.
  DUInt<uint64_t, 64> permraw;
  permraw ^= permraw;
  static const auto szperm(sizeof(DUInt<uint64_t, 64>) / sizeof(uint8_t));
  for(int i = 0; i < szperm; i ++) {
     uint8_t buf;
     in.read(reinterpret_cast<char*>(&buf), sizeof(uint8_t));
     permraw |= DUInt<uint64_t, 64>(uint64_t(buf)) << (i * 8);
  }
  std::vector<int> isort;
  std::vector<int> idxs;
  uint32_t punch;
  isort.resize(out.size(), 0);
  idxs.resize(out.size(), 0);
  punch ^= punch;
  for(int i = 2; i <= out.size() + 1; i ++) {
    isort[i - 2] = int(permraw % DUInt<uint64_t, 64>(i));
    permraw     /= DUInt<uint64_t, 64>(i);
  }
  assert(! permraw);
  for(int i = 0; i < isort.size(); i ++) {
    int cnt(i < isort.size() - 1 ? isort[isort.size() - 1 - i] : 0);
    int j(0);
    assert(0 <= cnt && cnt < isort.size() - i);
    for(int cnt = 0;
            cnt <= isort[isort.size() - 1 - i] && j < 32;
            j ++)
      if(! (punch & (j ? (1 << j) : 1))) {
        if(++ cnt <= isort[isort.size() - 1 - i])
          continue;
        else
          break;
      }
    punch   |= j ? (1 << j) : 1;
    idxs[i]  = j;
  }
  assert(!((punch + 1) & 0xffffffff));
  // decode after them.
  uint8_t work(0);
  int     bidx(0);
  for(int i = 0; i < out.size(); i ++) {
    auto& outs(out[idxs[i]]);
    if(! i) {
      in.read(reinterpret_cast<char*>(&work), sizeof(uint8_t));
      outs  = work;
      in.read(reinterpret_cast<char*>(&work), sizeof(uint8_t));
      outs |= work << 8;
    } else {
      const auto  bits0(min(countMSB(abs(out[idxs[i - 1]])) + 2, 16));
            auto  bits(bits0);
      outs ^= outs;
      while(0 < bits) {
        const auto residue(8 - bidx);
        if(! bidx) in.read(reinterpret_cast<char*>(&work), sizeof(uint8_t));
        outs |= (bits > residue ? (work << (bits - residue))
                 : (bits == residue ? work : (work >> (residue - bits)))) &
                ((1 << bits) - 1);
        bidx += min(bits, residue);
        bidx %= 8;
        bits -= residue;
      }
      if(outs & (1 << (bits0 - 1)))
        outs |= (- 1) & ~ ((1 << bits0) - 1);
    }
    if(! outs) {
      for(int j = i + 1; j < idxs.size(); j ++)
        out[idxs[j]] = 0;
      break;
    }
  }
  return out;
}

int main(int argc, char* argv[]) {
  if(argc < 2) usage();
  const auto pnext(P0<sfloat>().next(32));
  if(std::string(argv[1]) == std::string("-e")) {
    try {
      SimpleVector<int16_t> work(32);
      SimpleVector<int16_t> pbuf(32);
      for(int i = 0; i < pbuf.size(); i ++)
        std::cin.read(reinterpret_cast<char*>(&pbuf[i]), sizeof(int16_t));
      const auto out(blockout(pbuf));
      for(int i = 0; i < out.size(); i ++)
        std::cout.write(reinterpret_cast<const char*>(&out[i]), sizeof(uint8_t));
      while(! std::cin.eof() && ! std::cin.bad()) {
        for(int i = 0; i < work.size(); i ++)
          std::cin.read(reinterpret_cast<char*>(&work[i]), sizeof(int16_t));
        for(int i = 0; i < work.size(); i ++) {
          const auto bbuf(work[i]);
          auto pwork(pbuf.template cast<sfloat>());
          for(int j = 0; j < pwork.size() - 1; j ++)
            pwork[j] += pwork[j + 1];
          pwork[pwork.size() - 1] += pwork[pwork.size() - 1];
          work[i] -= int(pnext.dot(pwork) / sfloat(2));
          for(int j = 1; j < pbuf.size(); j ++)
            pbuf[j - 1] = pbuf[j];
          pbuf[pbuf.size() - 1] = bbuf;
        }
        const auto out(blockout(work));
        for(int i = 0; i < out.size(); i ++)
          std::cout.write(reinterpret_cast<const char*>(&out[i]), sizeof(uint8_t));
        std::cout.flush();
      }
    } catch(...) {
      std::cerr << "XXX: some error had occur." << std::endl;
    }
  } else if(std::string(argv[1]) == std::string("-d")) {
    try {
      auto pbuf(blockin(std::cin));
      for(int i = 0; i < pbuf.size(); i ++)
        std::cout.write(reinterpret_cast<char*>(&pbuf[i]), sizeof(int16_t));
      while(! std::cin.eof() && ! std::cin.bad()) {
        auto work(blockin(std::cin));
        for(int i = 0; i < work.size(); i ++) {
          auto pwork(pbuf.template cast<sfloat>());
          for(int j = 0; j < pwork.size() - 1; j ++)
            pwork[j] += pwork[j + 1];
          pwork[pwork.size() - 1] += pwork[pwork.size() - 1];
          work[i] += int(pnext.dot(pwork) / sfloat(2));
          for(int j = 1; j < pbuf.size(); j ++)
            pbuf[j - 1] = pbuf[j];
          pbuf[pbuf.size() - 1] = work[i];
          std::cout.write(reinterpret_cast<char*>(&work[i]), sizeof(int16_t));
        }
        std::cout.flush();
      }
    } catch(...) {
      std::cerr << "XXX: some error had occur." << std::endl;
    }
  } else usage();
  return 0;
}

