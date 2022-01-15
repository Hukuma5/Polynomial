// Copyright 2018 Your Name <your_email>

#ifndef INCLUDE_POLYNOMIAL_HPP_
#define INCLUDE_POLYNOMIAL_HPP_

#include <math.h>
#include <limits>
#include <type_traits>
#include <vector>

using namespace std;
template <class T>
class polinom {
  static_assert(std::is_arithmetic<T>::value, "Not arithmetic type");
  int N;  //размер вектора и степень полинома + 1
  vector<T> c;

 public:
  polinom() {  //конструктор по уммолчанию
    N = 0;
    c = {};
  }
  polinom(int new_N) {
    N = new_N;
    c.resize(N);
    for (int i = 0; i < N; i++) {
      c[i] = 0;
    }
  }
  polinom(vector<T> k) {  //конструктор по вектору коеф.
    N = k.size();
    c.resize(N);
    for (int i = 0; i < static_cast<int> (N); i++) {
      c[i] = k[i];
    }
  }
  polinom operator+(polinom& rhs) {
    if (rhs.c.size() > c.size()) {
      N = rhs.c.size();
      c.resize(N);
    } else {
      N = c.size();
      rhs.c.resize(N);
    }
    polinom sum(N);
    for (int i = 0; i < static_cast<int> (N); i++) {
      sum.c[i] = c[i] + rhs.c[i];
    }
    return sum;
  }
  polinom& operator+=(polinom& rhs) {
    if (rhs.c.size() > c.size()) {
      N = rhs.c.size();
      c.resize(N);
    } else {
      N = c.size();
      rhs.c.resize(N);
    }
    for (int i = 0; i < static_cast<int> (N); i++) {
      c[i] += rhs.c[i];
    }
    return *this;
  }
  polinom operator-(polinom& rhs) {
    if (rhs.c.size() > c.size()) {
      N = rhs.c.size();
      c.resize(N);
    } else {
      N = c.size();
      rhs.c.resize(N);
    }
    polinom raz(N);
    for (int i = 0; i < static_cast<int> (N); i++) {
      raz.c[i] = c[i] - rhs.c[i];
    }
    return raz;
  }
  polinom operator-=(polinom& rhs) {
    if (rhs.c.size() > c.size()) {
      N = rhs.c.size();
      c.resize(N);
    } else {
      N = c.size();
      rhs.c.resize(N);
    }
    for (int i = 0; i < static_cast<int> (N); i++) {
      c[i] -= rhs.c[i];
    }
    return *this;
  }
  polinom operator*(const polinom& rhs) const{
    polinom temp(N + rhs.N);
    for (int i = 0; i < static_cast<int> (N); i++) {
      for (int j = 0; j < static_cast<int> (rhs.N); j++) {
        temp.c[i + j] += c[i] * rhs.c[j];
      }
    }
    temp.N = N + rhs.N - 1;
    return temp;
  }
  polinom& operator*=(const polinom& rhs) {
    polinom temp(N + rhs.N);
    for (int i = 0; i < static_cast<int> (N); i++) {
      for (int j = 0; j < static_cast<int> (rhs.N); j++) {
        temp.c[i + j] += c[i] * rhs.c[j];
      }
    }
    temp.N = N + rhs.N - 1;
    N = temp.N;
    c = temp.c;
    return *this;
  }
  polinom operator*(T k) const{
    polinom r(c);
    for (int i = 0; i < c.size(); i++) {
      r.c[i] = c[i] * k;
    }
    return r;
  }
  polinom& operator=(const polinom& a) {
    N = a.N;
    c = a.c;
    return *this;
  }
  polinom operator/(const polinom& rhs) const {
    if (c.size() < rhs.c.size()) {
      polinom<T> R(1);
      R.c[0] = 0;
      return R;
    } else {
      if (c.size() == rhs.c.size()) {
        double f = c[0] / rhs.c[0];
        polinom<T> R(1);
        R.c[0] = f;
        return R;
      } else {
        polinom<T> ost(static_cast<int>(c.size()));
        int h = static_cast<int>(c.size() - rhs.c.size()) + 1;
        polinom<T> chast(h);
        for (int i = 0; i < static_cast<int>(c.size()); i++) {
          ost.c[i] = c[i];
        }
        for (int i = 0; i < static_cast<int>(chast.c.size()); i++) {
          double m = ost.c[0] / rhs.c[0];
          for (int j = 0; j < static_cast<int>(rhs.c.size()); j++) {
            ost.c[j] -= m * rhs.c[j];
          }
          chast.c[i] = m;
          auto iter = ost.c.cbegin();
          ost.c.erase(iter);
        }
        return chast;
      }
    }
  }
  polinom operator/=(const polinom& rhs) {
    if (c.size() < rhs.c.size()) {
      c.resize(1);
      c[0] = 0;
      return *this;
    } else {
      if (c.size() == rhs.c.size()) {
        double f = c[0] / rhs.c[0];
        polinom<T> R(1);
        R.c[0] = f;
        return R;
      } else {
        polinom<T> ost(static_cast<int>(c.size()));
        int h = static_cast<int>(c.size() - rhs.c.size()) + 1;
        polinom<T> chast(h);
        for (int i = 0; i < static_cast<int>(c.size()); i++) {
          ost.c[i] = c[i];
        }
        for (int i = 0; i < static_cast<int>(chast.c.size()); i++) {
          double m = ost.c[0] / rhs.c[0];
          for (int j = 0; j < static_cast<int>(rhs.c.size()); j++) {
            ost.c[j] -= m * rhs.c[j];
          }
          chast.c[i] = m;
          auto iter = ost.c.cbegin();
          ost.c.erase(iter);
        }
        c.clear();
        c.resize(chast.c.size());
        for (int i = 0; i < static_cast<int>(c.size()); i++) {
          c[i] = chast.c[i];
        }
        return *this;
      }
    }
  }
  polinom operator%(const polinom& rhs) const {
    if (c.size() < rhs.c.size()) {
      return *this;
    } else {
      if (c.size() == rhs.c.size()) {
        double f = c[0] / rhs.c[0];
        polinom<T> R(static_cast<int>(c.size()));
        for (int i = 0; i < static_cast<int>(c.size()); i++) {
          R.c[i] = c[i] - f * rhs.c[i];
        }
        while (R.c[0] < std::numeric_limits<double>::epsilon()) {
          auto iter = R.c.cbegin();
          R.c.erase(iter);
        }
        return R;
      } else {
        polinom<T> ost(static_cast<int>(c.size()));
        int h = static_cast<int>(c.size() - rhs.c.size()) + 1;
        polinom<T> chast(h);
        for (int i = 0; i < static_cast<int>(c.size()); i++) {
          ost.c[i] = c[i];
        }
        for (int i = 0; i < static_cast<int>(chast.c.size()); i++) {
          double m = ost.c[0] / rhs.c[0];
          for (int j = 0; j < static_cast<int>(rhs.c.size()); j++) {
            ost.c[j] -= m * rhs.c[j];
          }
          chast.c[i] = m;
          auto iter = ost.c.cbegin();
          ost.c.erase(iter);
        }
        return ost;
      }
    }
  }
  polinom operator%=(const polinom& rhs) {
    if (c.size() < rhs.c.size()) {
      return *this;
    } else {
      if (c.size() == rhs.c.size()) {
        double f = c[0] / rhs.c[0];
        polinom<T> R(static_cast<int>(c.size()));
        for (int i = 0; i < static_cast<int>(c.size()); i++) {
          R.c[i] = c[i] - f * rhs.c[i];
        }
        while (R.c[0] < std::numeric_limits<double>::epsilon()) {
          auto iter = R.c.cbegin();
          R.c.erase(iter);
        }
        return R;
      } else {
        polinom<T> ost(static_cast<int>(c.size()));
        int h = static_cast<int>(c.size() - rhs.c.size()) + 1;
        polinom<T> chast(h);
        for (int i = 0; i < static_cast<int>(c.size()); i++) {
          ost.c[i] = c[i];
        }
        for (int i = 0; i < static_cast<int>(chast.c.size()); i++) {
          double m = ost.c[0] / rhs.c[0];
          for (int j = 0; j < static_cast<int>(rhs.c.size()); j++) {
            ost.c[j] -= m * rhs.c[j];
          }
          chast.c[i] = m;
          auto iter = ost.c.cbegin();
          ost.c.erase(iter);
        }
        c.clear();
        c.resize(ost.c.size());
        for (int i = 0; i < static_cast<int>(c.size()); i++) {
          c[i] = ost.c[i];
        }
        return *this;
      }
    }
  }
  int get_N() const { return N; }
  T get_ci(int i) const { return c[i]; }
  void set_c(int j, T new_c) {
    for (int i = 0; i < static_cast<int> (N); i++) {
      if (i == j) c[i] = new_c;
    }
  }

  T get_x(T x) {
    T rez = 0;
    for (int i = 0; i < static_cast<int> (N); i++) {
      rez += c[i] * pow(x, i);
    }
    return rez;
  }
  ~polinom() { 
	  c.clear();
  }
};

#endif  // INCLUDE_POLYNOMIAL_HPP_
