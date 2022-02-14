/*
 BSD 3-Clause License

Copyright (c) 2019-2021, bitsofcotton (kazunobu watatsu)
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
#if !defined(_P0_)

using std::ifstream;
using std::ofstream;
using std::cerr;
using std::flush;
using std::endl;

template <typename T> SimpleVector<T> pnext(const int& size, const int& step = 1) {
  SimpleVector<T> p;
  if(size <= 1) {
    p.resize(1);
    p[0] = T(1);
  } else if(size <= 2) {
    p.resize(2);
    p[0] = T(1) / T(2);
    p[1] = T(1) / T(2) + T(1);
    p   /= T(2);
  } else {
    const auto file(std::string("./.cache/lieonn/next-") +
      std::to_string(size) + std::string("-") + std::to_string(step) +
#if defined(_FLOAT_BITS_)
      std::string("-") + std::to_string(_FLOAT_BITS_)
#else
      std::string("-ld")
#endif
    );
    ifstream cache(file.c_str());
    if(cache.is_open()) {
      cache >> p;
      cache.close();
    } else {
      p = taylor(size, T(size + step - 1));
      cerr << "." << flush;
      const auto pp(pnext<T>(size - 1, step));
      for(int j = 0; j < pp.size(); j ++)
        p[j - pp.size() + p.size()] += pp[j] * T(size - 1);
      p /= T(size);
      ofstream ocache(file.c_str());
      ocache << p << endl;
      ocache.close();
    }
  }
  assert(p.size() == size);
  return p;
}

template <typename T, typename feeder> class P0 {
public:
  typedef SimpleVector<T> Vec;
  inline P0() { ; }
  inline P0(const int& size, const int& step = 1) {
    f = feeder(size);
    p = pnext<T>(size, step);
  }
  inline ~P0() { ; };
  inline T next(const T& in) {
    const auto& ff(f.next(in));
    static const T zero(int(0));
    return f.full ? p.dot(ff) : zero;
  }
  Vec p;
  feeder f;
};

template <typename T, typename P> class northPole {
public:
  inline northPole() { ; }
  inline northPole(P&& p, const T& r = T(int(0))) {
    this->p = p;
    this->r = r;
    M0 = M1 = T(int(0));
  }
  inline ~northPole() { ; }
  inline T next(const T& in) {
    static const T zero(int(0));
    static const T one(int(1));
    if(! isfinite(in)) return in;
    if(M0 < abs(in)) M0 = abs(in) * T(int(2));
    if(in == zero || M0 == zero) return in;
    auto s(r == zero ? one / atan(in) : one / atan(in * r / M0));
    if(! isfinite(s)) return in;
    if(M1 < abs(s)) M1 = abs(s) * T(int(2));
    if(s  == zero || M1 == zero) return in;
    const auto pn(max(- atan(r), min(atan(r), p.next(atan(r == zero ? s : s * r / M1)))));
    if(pn == zero || ! isfinite(pn)) return in;
    auto res(tan(max(- atan(r), min(atan(r), one / (tan(pn) * (r == zero ? one : M1 / r))))) * (r == zero ? one : M0 / r));
    if(isfinite(res)) return res;
    return in;
  }
  P p;
  T r;
  T M0;
  T M1;
};

template <typename T, typename P> class sumChain {
public:
  inline sumChain() { ; }
  inline sumChain(P&& p) { S = T(int(0)); this->p = p; }
  inline ~sumChain() { ; }
  inline T next(const T& in) {
    auto res(- S); return res += p.next(S += in);
  }
  T S;
  P p;
};

#define _P0_
#endif

