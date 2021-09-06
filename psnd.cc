#include <cstdio>
#include <vector>
#include <map>
#include <complex>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cstring>
#include <assert.h>
#include "lieonn.hh"
typedef myfloat num_t;
typedef num_t sfloat;
#include "p0.hh"

const auto pblocks(80);

/*
const auto blocks(128);
const auto szperm((718 + 7) / 8);
typedef DUInt<DUInt<DUInt<DUInt<uint64_t, 64>, 128>, 256>, 512> perm_t;
typedef DUInt<uint64_t, 64> punch_t;

const auto blocks(64);
const auto szperm((297 + 7) / 8);
typedef DUInt<DUInt<DUInt<uint64_t, 64>, 128>, 256> perm_t;
typedef uint64_t punch_t;

const auto blocks(16);
const auto szperm((46 + 7) / 8);
typedef uint64_t perm_t;
typedef uint16_t punch_t;
*/
const auto blocks(7);
const auto szperm(2);
typedef uint32_t perm_t;
typedef uint8_t punch_t;

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
  assert(buf.size() == blocks);
  std::vector<std::pair<std::pair<uint16_t, int16_t>, int> > s0;
  s0.reserve(buf.size());
  for(int i = 0; i < buf.size(); i ++)
    s0.emplace_back(std::make_pair(std::make_pair(abs(buf[i]), buf[i]), i));
  std::sort(s0.begin(), s0.end());
  perm_t permraw;
  permraw ^= permraw;
  std::vector<uint8_t> outbase;
  punch_t punch;
  punch ^= punch;
  int  bidx(0);
  bool flag(false);
  for(int i = 0; i < buf.size(); i ++) {
    const auto& tb(s0[buf.size() - i - 1]);
    int cnt(0);
    for(int j = 0; j < tb.second; j ++)
      if(! (punch & (j ? (punch_t(1) << j) : punch_t(1)))) cnt ++;
    punch |= tb.second ? (punch_t(1) << tb.second) : punch_t(1);
    assert(0 <= cnt && cnt < buf.size() - i);
    if(i < buf.size() - 1) {
      permraw += perm_t(int(cnt));
      permraw *= perm_t(int(buf.size() - i - 1));
    }
    if(flag) continue;
    if(! i) {
      outbase.emplace_back(tb.first.second & 0xff);
      outbase.emplace_back((tb.first.second >> 8) & 0xff);
    } else {
      auto bits(countMSB(s0[buf.size() - i].first.first) + 2);
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
  //assert(!(++ punch));
  SimpleVector<uint8_t> out(outbase.size() + szperm);
  for(int i = 0; i < szperm; i ++)
    out[i] = uint8_t(int(i ? (permraw >> (i * 8)) : permraw) & 0xff);
  for(int i = 0; i < outbase.size(); i ++)
    out[i + szperm] = outbase[i];
  return out;
}

std::pair<SimpleVector<int16_t>, bool> blockin(std::istream& in) {
  std::pair<SimpleVector<int16_t>, bool> out;
  out.first.resize(blocks);
  out.second = false;
  // decode permraw.
  perm_t permraw;
  permraw ^= permraw;
  for(int i = 0; i < szperm; i ++) {
     uint8_t buf;
     in.read(reinterpret_cast<char*>(&buf), sizeof(uint8_t));
     if(i == szperm - 1 && (buf & 0x80)) {
       buf = ~ buf;
       out.second = true;
     }
     permraw |= perm_t(uint64_t(buf)) << (i * 8);
  }
  std::vector<int> isort;
  std::vector<int> idxs;
  punch_t punch;
  isort.resize(out.first.size(), 0);
  idxs.resize(out.first.size(), 0);
  punch ^= punch;
  for(int i = 2; i <= out.first.size(); i ++) {
    isort[i - 1] = int(permraw % perm_t(i));
    permraw     /= perm_t(i);
  }
  // XXX: clang is something buggy:
  // assert(! permraw);
  for(int i = 0; i < isort.size(); i ++) {
    int cnt(i < isort.size() - 1 ? isort[isort.size() - 1 - i] : 0);
    int j(0);
    assert(0 <= cnt && cnt < isort.size() - i);
    for(int cnt = 0;
            cnt <= isort[isort.size() - 1 - i] && j < blocks;
            j ++)
      if(! (punch & (j ? (punch_t(1) << j) : punch_t(1)))) {
        if(++ cnt <= isort[isort.size() - 1 - i])
          continue;
        else
          break;
      }
    punch   |= j ? (punch_t(1) << j) : punch_t(1);
    idxs[i]  = j;
  }
  //assert(!(++ punch));
  // decode after them.
  uint8_t work(0);
  int     bidx(0);
  for(int i = 0; i < out.first.size(); i ++) {
    auto& outs(out.first[idxs[i]]);
    if(! i) {
      in.read(reinterpret_cast<char*>(&work), sizeof(uint8_t));
      outs  = work;
      in.read(reinterpret_cast<char*>(&work), sizeof(uint8_t));
      outs |= work << 8;
    } else {
      const auto  bits0(countMSB(abs(out.first[idxs[i - 1]])) + 2);
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
        out.first[idxs[j]] = 0;
      break;
    }
  }
  return out;
}

int main(int argc, char* argv[]) {
  if(argc < 2) usage();
  P0<sfloat, sumFeeder<sfloat, idFeeder<sfloat> > > p(pblocks);
  if(std::string(argv[1]) == std::string("-e")) {
    try {
      SimpleVector<int16_t> work(blocks);
      int nbuf(0);
      for(int i = 0; i < work.size(); i ++) {
        std::cin.read(reinterpret_cast<char*>(&work[i]), sizeof(int16_t));
        nbuf = int(p.next(sfloat(int(work[i]))));
      }
      const auto out(blockout(work));
      for(int i = 0; i < out.size(); i ++)
        std::cout.write(reinterpret_cast<const char*>(&out[i]), sizeof(uint8_t));
      while(! std::cin.eof() && ! std::cin.bad()) {
        for(int i = 0; i < work.size(); i ++)
          std::cin.read(reinterpret_cast<char*>(&work[i]), sizeof(int16_t));
        auto out0(blockout(work));
        out0[szperm - 1] = ~ out0[szperm - 1];
        for(int i = 0; i < work.size(); i ++) {
          const auto bbuf(work[i]);
          work[i] -= nbuf;
          nbuf     = int(p.next(sfloat(int(bbuf))));
        }
        const auto out(blockout(work));
        if(out.size() <= out0.size())
          for(int i = 0; i < out.size(); i ++)
            std::cout.write(reinterpret_cast<const char*>(&out[i]), sizeof(uint8_t));
        else
          for(int i = 0; i < out0.size(); i ++)
            std::cout.write(reinterpret_cast<const char*>(&out0[i]), sizeof(uint8_t));
        std::cout.flush();
      }
    } catch(...) {
      std::cerr << "XXX: some error had occur." << std::endl;
    }
  } else if(std::string(argv[1]) == std::string("-d")) {
    try {
      const auto w(blockin(std::cin));
      auto work(w.first);
      int  nbuf(0);
      for(int i = 0; i < work.size(); i ++) {
        std::cout.write(reinterpret_cast<char*>(&work[i]), sizeof(int16_t));
        nbuf = int(p.next(sfloat(int(work[i]))));
      }
      while(! std::cin.eof() && ! std::cin.bad()) {
        const auto w(blockin(std::cin));
        work = w.first;
        if(! w.second)
          for(int i = 0; i < work.size(); i ++) {
            work[i] += nbuf;
            nbuf     = int(p.next(sfloat(int(work[i]))));
            std::cout.write(reinterpret_cast<char*>(&work[i]), sizeof(int16_t));
          }
        else
          for(int i = 0; i < work.size(); i ++) {
            nbuf = int(p.next(sfloat(int(work[i]))));
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

