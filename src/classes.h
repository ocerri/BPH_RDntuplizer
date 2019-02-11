#include <string>
#include <map>
#include <utility>

#include "DataFormats/Common/interface/Wrapper.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"


namespace {
  struct dictionary {
    std::map<std::string, float> dummy_i1;
    edm::Wrapper<std::map<std::string,float>> dummy_wi1;

    std::vector<string> dummy_s2;
    edm::Wrapper<std::vector<string>> dummy_w2;

    // pat::PackedCandidate dummy_ci2;
    // std::pair<const std::basic_string<char,std::char_traits<char> >,pat::PackedCandidate> dummy_pi2;
    // std::pair<std::string, pat::PackedCandidate> dummy_pi2;
    // std::map<std::string, pat::PackedCandidate> dummy_mi2;
    // edm::Wrapper<std::map<std::string, pat::PackedCandidate>> dummy_wi2;
  };
}
