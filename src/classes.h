#include <string>
#include <map>

#include "DataFormats/Common/interface/Wrapper.h"

namespace {
  struct dictionary {
    std::map<std::string, float> dummy_i1;
    edm::Wrapper<std::map<std::string,float>> dummy_wi1;
  };
}
