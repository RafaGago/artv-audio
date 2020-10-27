#pragma once

#include <utility>
#include <vector>

namespace artv { namespace dragonfly {

struct compile_firewall {
public:
  compile_firewall() {};
  compile_firewall (compile_firewall const&) = delete;
  compile_firewall (compile_firewall&& other)
  {
    *this = std::forward<compile_firewall> (other);
  }

  compile_firewall& operator= (compile_firewall const&) = delete;
  compile_firewall& operator= (compile_firewall&& other)
  {
    if (this == &other) {
      return *this;
    }
    _dsp       = other._dsp;
    other._dsp = nullptr;
    return *this;
  }
  virtual ~compile_firewall() {};
  virtual void set_parameter (unsigned id, float v)   = 0;
  virtual void reset (unsigned samplerate)            = 0;
  virtual void process (float** io, unsigned samples) = 0;

protected:
  void* _dsp = nullptr;
};

class er_compile_firewall : public compile_firewall {
public:
  er_compile_firewall();
  ~er_compile_firewall();

  er_compile_firewall (er_compile_firewall const&) = delete;
  er_compile_firewall (er_compile_firewall&&)      = default;
  er_compile_firewall& operator= (er_compile_firewall const&) = delete;
  er_compile_firewall& operator= (er_compile_firewall&&) = default;

  void set_parameter (unsigned id, float v) override;
  void reset (unsigned samplerate) override;
  void process (float** io, unsigned samples) override;
};

class plate_compile_firewall : public compile_firewall {
public:
  plate_compile_firewall();
  ~plate_compile_firewall();

  plate_compile_firewall (plate_compile_firewall const&) = delete;
  plate_compile_firewall (plate_compile_firewall&&)      = default;
  plate_compile_firewall& operator= (plate_compile_firewall const&) = delete;
  plate_compile_firewall& operator= (plate_compile_firewall&&) = default;

  void set_parameter (unsigned id, float v) override;
  void reset (unsigned samplerate) override;
  void process (float** io, unsigned samples) override;
};

class hall_compile_firewall : public compile_firewall {
public:
  hall_compile_firewall();
  ~hall_compile_firewall();

  hall_compile_firewall (hall_compile_firewall const&) = delete;
  hall_compile_firewall (hall_compile_firewall&&)      = default;
  hall_compile_firewall& operator= (hall_compile_firewall const&) = delete;
  hall_compile_firewall& operator= (hall_compile_firewall&&) = default;

  void set_parameter (unsigned id, float v) override;
  void reset (unsigned samplerate) override;
  void process (float** io, unsigned samples) override;
};

struct room_compile_firewall : public compile_firewall {
public:
  room_compile_firewall();
  ~room_compile_firewall();

  room_compile_firewall (room_compile_firewall const&) = delete;
  room_compile_firewall (room_compile_firewall&&)      = default;
  room_compile_firewall& operator= (room_compile_firewall const&) = delete;
  room_compile_firewall& operator= (room_compile_firewall&&) = default;

  void set_parameter (unsigned id, float v) override;
  void reset (unsigned samplerate) override;
  void process (float** io, unsigned samples) override;
};

}}; // namespace artv::dragonfly
