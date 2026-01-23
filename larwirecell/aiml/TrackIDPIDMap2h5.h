#ifndef WIRECELL_AIML_TRACKIDPIDMAP2H5
#define WIRECELL_AIML_TRACKIDPIDMAP2H5

#include "WireCellAux/Logger.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IFrameFilter.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larwirecell/Interfaces/IArtEventVisitor.h"

#include <hdf5.h>
#include <string>
#include <unordered_map>
#include <vector>

namespace WireCell::AIML {
  class TrackIDPIDMap2h5 : public Aux::Logger,
                           public wcls::IArtEventVisitor,
                           public IFrameFilter,
                           public IConfigurable {
  public:
    TrackIDPIDMap2h5();
    ~TrackIDPIDMap2h5() override;

    // IFrameFilter
    bool operator()(const input_pointer& in, output_pointer& out) override;

    // IConfigurable
    void configure(const WireCell::Configuration& config) override;
    WireCell::Configuration default_configuration() const override;

    // IArtEventVisitor
    void visit(art::Event& event) override;

  private:
    void cache_simchannels(const std::vector<sim::SimChannel>& simchs);
    void populate_trackid_pid_map();
    void clear_cache();
    void ensure_file();
    void write_mapping(int frame_ident);

    std::string m_simchannel_label;
    std::string m_output_file;

    std::vector<sim::SimChannel> m_simchannels;
    std::unordered_map<int, int> m_trackid_to_pid;

    hid_t m_file;
  };
}

#endif
