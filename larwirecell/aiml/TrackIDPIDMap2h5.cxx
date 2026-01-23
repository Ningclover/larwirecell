#include "TrackIDPIDMap2h5.h"
#include "WireCellUtil/Configuration.h"
#include "WireCellUtil/NamedFactory.h"

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "cetlib_except/exception.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "nusimdata/SimulationBase/MCParticle.h"

#include <algorithm>
#include <hdf5.h>

WIRECELL_FACTORY(TrackIDPIDMap2h5,
                 WireCell::AIML::TrackIDPIDMap2h5,
                 WireCell::INamed,
                 WireCell::IFrameFilter,
                 WireCell::IConfigurable)

using namespace WireCell;

AIML::TrackIDPIDMap2h5::TrackIDPIDMap2h5()
  : Aux::Logger("TrackIDPIDMap2h5", "aiml")
  , m_simchannel_label("tpcrawdecoder:simpleSC")
  , m_output_file("trackid_pid_map.h5")
  , m_file(-1)
{}

AIML::TrackIDPIDMap2h5::~TrackIDPIDMap2h5()
{
  if (m_file >= 0) {
    H5Fclose(m_file);
    m_file = -1;
  }
}

Configuration AIML::TrackIDPIDMap2h5::default_configuration() const
{
  Configuration cfg;
  cfg["simchannel_label"] = m_simchannel_label;
  cfg["output_file"] = m_output_file;
  return cfg;
}

void AIML::TrackIDPIDMap2h5::configure(const Configuration& cfg)
{
  m_simchannel_label = get(cfg, "simchannel_label", m_simchannel_label);
  m_output_file = get(cfg, "output_file", m_output_file);

  clear_cache();
  if (m_file >= 0) {
    H5Fclose(m_file);
    m_file = -1;
  }
}

void AIML::TrackIDPIDMap2h5::visit(art::Event& event)
{
  log->info("TrackIDPIDMap2h5::visit called for event: {}", event.event());
  clear_cache();

  if (m_simchannel_label.empty()) {
    log->warn("TrackIDPIDMap2h5: SimChannel label not configured; skipping");
    return;
  }

  log->debug("TrackIDPIDMap2h5: Fetching SimChannel with label '{}'", m_simchannel_label);
  art::Handle<std::vector<sim::SimChannel>> handle;
  if (!event.getByLabel(art::InputTag{m_simchannel_label}, handle)) {
    log->warn("TrackIDPIDMap2h5 failed to fetch SimChannel with label '{}'", m_simchannel_label);
    return;
  }

  log->debug("TrackIDPIDMap2h5: Got {} SimChannels", handle->size());
  cache_simchannels(*handle);
  populate_trackid_pid_map();
  log->info("TrackIDPIDMap2h5 cached {} track->pid entries", m_trackid_to_pid.size());
}

bool AIML::TrackIDPIDMap2h5::operator()(const input_pointer& in, output_pointer& out)
{
  out = in;
  if (!in) {
    log->debug("TrackIDPIDMap2h5: null input frame");
    return true;
  }

  log->info("TrackIDPIDMap2h5::operator() called for frame {}, have {} track->pid entries",
            in->ident(),
            m_trackid_to_pid.size());

  if (!m_trackid_to_pid.empty()) {
    ensure_file();
    write_mapping(in->ident());
  }
  else {
    log->warn("TrackIDPIDMap2h5: no track ID to PID mapping to write for frame {} - was visit() called?",
              in->ident());
  }

  return true;
}

void AIML::TrackIDPIDMap2h5::cache_simchannels(const std::vector<sim::SimChannel>& simchs)
{
  m_simchannels = simchs;
  m_trackid_to_pid.clear();
}

void AIML::TrackIDPIDMap2h5::populate_trackid_pid_map()
{
  m_trackid_to_pid.clear();
  if (m_simchannels.empty()) { return; }

  try {
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    for (auto const& sc : m_simchannels) {
      for (auto const& tdc_entry : sc.TDCIDEMap()) {
        for (auto const& ide : tdc_entry.second) {
          const int track_id = ide.trackID;
          if (track_id == 0) { continue; }
          if (m_trackid_to_pid.find(track_id) != m_trackid_to_pid.end()) { continue; }

          int pid = 0;
          try {
            auto const particle = pi_serv->TrackIdToParticle_P(track_id);
            if (particle) { pid = particle->PdgCode(); }
          }
          catch (const cet::exception& ex) {
            log->debug("TrackIDPIDMap2h5: failed to fetch MCParticle for track {}: {}",
                       track_id,
                       ex.what());
          }
          m_trackid_to_pid.emplace(track_id, pid);
        }
      }
    }
  }
  catch (const cet::exception& ex) {
    log->warn("TrackIDPIDMap2h5: ParticleInventoryService unavailable: {}", ex.what());
  }
}

void AIML::TrackIDPIDMap2h5::clear_cache()
{
  m_simchannels.clear();
  m_trackid_to_pid.clear();
}

void AIML::TrackIDPIDMap2h5::ensure_file()
{
  if (m_file >= 0) { return; }

  m_file = H5Fopen(m_output_file.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  if (m_file < 0) {
    m_file = H5Fcreate(m_output_file.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (m_file < 0) {
      log->error("TrackIDPIDMap2h5 failed to open output file {}", m_output_file);
    }
  }
}

void AIML::TrackIDPIDMap2h5::write_mapping(int frame_ident)
{
  if (m_file < 0) { return; }
  if (m_trackid_to_pid.empty()) { return; }

  // Extract and sort track IDs for consistent ordering
  std::vector<int> track_ids;
  std::vector<int> pids;
  track_ids.reserve(m_trackid_to_pid.size());
  pids.reserve(m_trackid_to_pid.size());

  for (auto const& entry : m_trackid_to_pid) {
    track_ids.push_back(entry.first);
  }
  std::sort(track_ids.begin(), track_ids.end());

  for (int track_id : track_ids) {
    pids.push_back(m_trackid_to_pid.at(track_id));
  }

  const hsize_t dims[1] = {track_ids.size()};
  const std::string group_name = "/" + std::to_string(frame_ident);
  const std::string trackid_dset_name = group_name + "/track_ids";
  const std::string pid_dset_name = group_name + "/pids";

  // Create group if it doesn't exist
  hid_t lcpl = H5Pcreate(H5P_LINK_CREATE);
  H5Pset_create_intermediate_group(lcpl, 1);

  // Create and write track IDs dataset
  hid_t dataspace = H5Screate_simple(1, dims, nullptr);
  if (dataspace < 0) {
    log->warn("TrackIDPIDMap2h5 failed to create dataspace for track_ids");
    H5Pclose(lcpl);
    return;
  }

  hid_t trackid_dset = H5Dcreate2(
    m_file, trackid_dset_name.c_str(), H5T_NATIVE_INT, dataspace, lcpl, H5P_DEFAULT, H5P_DEFAULT);

  if (trackid_dset >= 0) {
    herr_t status = H5Dwrite(trackid_dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, track_ids.data());
    if (status < 0) {
      log->warn("TrackIDPIDMap2h5 failed to write track_ids dataset");
    }
    H5Dclose(trackid_dset);
  }
  else {
    log->warn("TrackIDPIDMap2h5 failed to create track_ids dataset {}", trackid_dset_name);
  }

  // Create and write PIDs dataset
  hid_t pid_dset = H5Dcreate2(
    m_file, pid_dset_name.c_str(), H5T_NATIVE_INT, dataspace, lcpl, H5P_DEFAULT, H5P_DEFAULT);

  if (pid_dset >= 0) {
    herr_t status = H5Dwrite(pid_dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, pids.data());
    if (status < 0) { log->warn("TrackIDPIDMap2h5 failed to write pids dataset"); }
    H5Dclose(pid_dset);
  }
  else {
    log->warn("TrackIDPIDMap2h5 failed to create pids dataset {}", pid_dset_name);
  }

  H5Sclose(dataspace);
  H5Pclose(lcpl);

  log->debug("TrackIDPIDMap2h5 wrote {} track ID to PID mappings to frame {}", track_ids.size(), frame_ident);
}
