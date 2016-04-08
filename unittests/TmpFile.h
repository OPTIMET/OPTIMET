#ifndef OPTIMET_TMPFILE_H
#include <cstdio>
#include <fstream>

namespace optimet {
//! Creates a file that deletes itself when it goes out of scope
class TmpFile {
public:
  TmpFile() : name_(std::tmpnam(nullptr)) {}
  ~TmpFile() { std::remove(filename().c_str()); }
  void write(std::string const &string) {
    std::ofstream stream(filename().c_str());
    stream << string;
    stream.close();
  }
  std::string filename() const { return name_; }

protected:
  std::string name_;
};
}
#endif
