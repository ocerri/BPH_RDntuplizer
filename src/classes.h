#include <string>
#include <map>
#include <utility>

#include "DataFormats/Common/interface/Wrapper.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

template<>
inline
bool edm::Wrapper<std::pair<const std::basic_string<char,std::char_traits<char> >,std::vector<float> >>::hasSwap_() const {
  return true;
}

template<>
inline
void edm::Wrapper<std::pair<const std::basic_string<char,std::char_traits<char> >,std::vector<float> >>::swapProduct_(WrapperBase* newProduct) {}


namespace {
  struct dictionary {
    std::map<std::string, float> dummy_i1;
    edm::Wrapper<std::map<std::string,float>> dummy_wi1;

    std::vector<string> dummy_s2;
    edm::Wrapper<std::vector<string>> dummy_w2;

    std::map<std::string, std::vector<float>> dummy_i3;
    edm::Wrapper<std::map<std::string,std::vector<float>>> dummy_wi3;

    std::char_traits<char> d4;
    std::basic_string<char,std::char_traits<char> > d5;
    std::pair<const std::basic_string<char,std::char_traits<char> >,std::vector<float> > d6;
    edm::Wrapper<std::pair<const std::basic_string<char,std::char_traits<char> >,std::vector<float> >> wd6;

    // pat::PackedCandidate dummy_ci2;
    // std::pair<const std::basic_string<char,std::char_traits<char> >,pat::PackedCandidate> dummy_pi2;
    // std::pair<std::string, pat::PackedCandidate> dummy_pi2;
    // std::map<std::string, pat::PackedCandidate> dummy_mi2;
    // edm::Wrapper<std::map<std::string, pat::PackedCandidate>> dummy_wi2;
  };
}
